!!****if* source/Driver/DriverMain/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***



subroutine Driver_sourceTerms(blockCount, blockList, dt, pass) 
    use Polytrope_interface, ONLY : Polytrope
    use Driver_data, ONLY: dr_simTime, dr_meshComm, dr_dt, dr_meshMe
    use Flame_interface, ONLY : Flame_step
    use Stir_interface, ONLY : Stir
    use Heat_interface, ONLY : Heat
    use Heatexchange_interface, ONLY : Heatexchange
    use Burn_interface, ONLY : Burn
    use Cool_interface, ONLY : Cool
    use Ionize_interface, ONLY : Ionize
    use EnergyDeposition_interface, ONLY : EnergyDeposition
    use Deleptonize_interface, ONLY : Deleptonize
    use Simulation_data, ONLY: sim_smallX, &
        sim_tRelax, sim_relaxRate, sim_fluffDampCoeff, sim_fluffDampCutoff, sim_accRadius, sim_accCoeff, &
        sim_fluidGamma, sim_softenRadius, sim_rotFac, sim_rotAngle, sim_tSpinup, obj_ipos, obj_radius, &
        sim_objMass, sim_msun, sim_xCenter, sim_yCenter, sim_zCenter, sim_cylinderScale, &
        sim_cylinderDensity, sim_cylinderTemperature, obj_mu, obj_gamc, stvec, obj_xn, &
        sim_cylinderMDot, sim_cylinderNCells, sim_ptMass, sim_cylinderRadius, &
        sim_kind, sim_windVelocity, sim_windMdot, sim_windTemperature, sim_windNCells, &
        sim_fixedPartTag, sim_windKernel, sim_cylinderType, sim_orbEcc, sim_periDist, &
        sim_tDelay, sim_mpoleVX, sim_mpoleVY, sim_mpoleVZ
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize, Grid_fillGuardCells
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
    !use Gravity_data, ONLY: grv_ptvec, grv_obvec, grv_ptmass, grv_exactvec, grv_optmass, grv_momacc, &
    !    grv_angmomacc, grv_eneracc, grv_massacc, grv_comPeakCut, grv_totmass, grv_totmass
    use gr_mpoleData, ONLY: gr_mpoleXcenterOfMass, gr_mpoleYcenterOfMass, gr_mpoleZcenterOfMass
    use Grid_data, ONLY: gr_smalle, gr_meshMe
    use Eos_interface, ONLY: Eos_wrapped
    use Eos_data, ONLY : eos_gasConstant
    use pt_sinkInterface, ONLY : pt_sinkGatherGlobal
    use Particles_sinkData, ONLY : localnpf, particles_global
    use Multitidal_interface, ONLY : Multitidal_findExtrema
    implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    real, intent(IN)    :: dt
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    integer, OPTIONAL, intent(IN) :: pass
    
    integer :: i, j, k, lb, ierr, istat
    integer :: sizeX, sizeY, sizeZ
    real    :: xx, yy, zz, dist, y2, z2, tot_mass, gtot_mass, peak_mass, gpeak_mass
    real    :: relax_rate, mass_acc, tot_mass_acc, gtot_mass_acc, tot_ener_acc, gtot_ener_acc
    real    :: distxy, vpara, vspin, x, y, z, vx, vy, vz, velprojy, velprojz, rotang
    real    :: tinitial, vol, ldenscut, denscut, extrema, newton, mcs, new_dens
    real    :: flow_dist, flow_vel, polyk, rho0, kb, mp, yr, mdot, perp_dist, flow_ecc
    real    :: damp_dens
  
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    real, dimension(6) :: pvec
    real, dimension(3) :: vvec
    real, dimension(:,:,:,:),pointer :: solnData
    real, dimension(MDIM) :: tot_avg_vel, peak_avg_vel, tot_mom_acc, gtot_mom_acc, tot_com_acc, gtot_com_acc, &
        tot_angmom_acc, gtot_angmom_acc, tot_mom, gtot_mom, peak_mom, gpeak_mom
    real, dimension(2*MDIM) :: pt_pos
  
    logical :: gcell = .true.
  
    tot_mass_acc = 0.d0
    tot_com_acc = 0.d0
    tot_mom_acc = 0.d0
    tot_angmom_acc = 0.d0
    tot_ener_acc = 0.d0
    gtot_mass_acc = 0.d0
    gtot_com_acc = 0.d0
    gtot_mom_acc = 0.d0
    gtot_angmom_acc = 0.d0
    gtot_ener_acc = 0.d0
    tot_mass = 0.d0
    tot_mom = 0.d0
    gtot_mass = 0.d0
    gtot_mom = 0.d0
    peak_mass = 0.d0
    peak_mom = 0.d0
    gpeak_mass = 0.d0
    gpeak_mom = 0.d0

    yr = 3.15569252E7
    call PhysicalConstants_get("Newton", newton)
    call PhysicalConstants_get("Boltzmann", kb)
    call PhysicalConstants_get("proton mass", mp)
    call RuntimeParameters_get('tinitial',tinitial)
    call Grid_getMinCellSize(mcs)

    rotang = PI/180.*sim_rotAngle

    if (sim_kind .eq. 'wind') then
        call pt_sinkGatherGlobal()
        do i = 1, localnpf
            if (idnint(particles_global(TAG_PART_PROP,i)) .ne. sim_fixedPartTag) then
                pvec(1:3) = particles_global(POSX_PART_PROP:POSZ_PART_PROP,i)
                pvec(4:6) = particles_global(VELX_PART_PROP:VELZ_PART_PROP,i)
            endif
        enddo
    endif

    if (sim_kind .eq. 'polytrope') then
        call Multitidal_findExtrema(DENS_VAR, 1, damp_dens)
        damp_dens = damp_dens*sim_fluffDampCutoff
    else
        damp_dens = sim_fluffDampCutoff
    endif

    if (dr_simTime .lt. tinitial + sim_tRelax) then
        relax_rate = max(0.0d0, dr_simTime - sim_tSpinup)/(sim_tRelax - sim_tSpinup)*(1.0 - sim_relaxRate) + sim_relaxRate 
        if (sim_kind .eq. "polytrope") then
            do lb = 1, blockCount
                call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
                sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
                allocate(xCoord(sizeX),stat=istat)
                sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
                allocate(yCoord(sizeY),stat=istat)
                sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
                allocate(zCoord(sizeZ),stat=istat)
  
                if (NDIM == 3) call Grid_getCellCoords&
                                    (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
                if (NDIM >= 2) call Grid_getCellCoords&
                                    (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
                call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
                call Grid_getBlkPtr(blockList(lb),solnData)

                do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                            ! subtract out COM motion
                            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) - sim_mpoleVX
                            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) - sim_mpoleVY
                            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) - sim_mpoleVZ

                            ! Need to work out math.
                            if (sim_rotFac .gt. 0.d0) then
                                x = xCoord(i) - gr_mpoleXcenterOfMass
                                y = (yCoord(j) - gr_mpoleYcenterOfMass)*dcos(rotang) - (zCoord(k) - gr_mpoleZcenterOfMass)*dsin(rotang)
                                !z = (yCoord(j) - gr_mpoleYcenterOfMass)*dsin(rotang) + (zCoord(k) - gr_mpoleZcenterOfMass)*dcos(rotang)
                                vx = solnData(VELX_VAR,i,j,k)
                                vy = solnData(VELY_VAR,i,j,k)*dcos(rotang) - solnData(VELZ_VAR,i,j,k)*dsin(rotang)
                                vz = solnData(VELY_VAR,i,j,k)*dsin(rotang) + solnData(VELZ_VAR,i,j,k)*dcos(rotang)

                                distxy = dsqrt(x**2 + y**2)
                                vpara = (x*vx + y*vy)/distxy
                                if (dr_simTime .lt. sim_tSpinup) then
                                    vspin = dr_simTime/sim_tSpinup*sim_rotFac*&
                                        min(dsqrt(newton*sim_objMass/obj_radius(obj_ipos)**3.d0)*distxy, &
                                            dsqrt(newton*sim_objMass/obj_radius(obj_ipos)))
                                    solnData(VELX_VAR,i,j,k) = -vspin*y/distxy - x/distxy*vpara*(1.d0 - sim_relaxRate)
                                    velprojy =                  vspin*x/distxy - y/distxy*vpara*(1.d0 - sim_relaxRate)
                                else
                                    solnData(VELX_VAR,i,j,k) = vx - x/distxy*vpara*(1.d0 - relax_rate)
                                    velprojy =                 vy - y/distxy*vpara*(1.d0 - relax_rate)
                                endif
                                velprojz = vz*relax_rate
                                solnData(VELY_VAR,i,j,k) =  dcos(rotang)*velprojy + dsin(rotang)*velprojz
                                solnData(VELZ_VAR,i,j,k) = -dsin(rotang)*velprojy + dcos(rotang)*velprojz
                            else
                                solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k)*relax_rate
                                solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)*relax_rate
                                solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k)*relax_rate
                            endif

                            solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                        enddo
                    enddo
                enddo
  
                call Grid_releaseBlkPtr(blockList(lb), solnData)
                deallocate(xCoord)
                deallocate(yCoord)
                deallocate(zCoord)
            enddo
        endif
    else
        !pt_pos = grv_exactvec - grv_obvec + grv_ptvec !Used for calculating properties of mass accreted onto point mass
        relax_rate = sim_fluffDampCoeff

        do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
            sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
            allocate(xCoord(sizeX),stat=istat)
            sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
            allocate(yCoord(sizeY),stat=istat)
            sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
            allocate(zCoord(sizeZ),stat=istat)
  
            if (NDIM == 3) call Grid_getCellCoords&
                                (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
            if (NDIM >= 2) call Grid_getCellCoords&
                                (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
            call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
            call Grid_getBlkPtr(blockList(lb),solnData)

            if (sim_kind .eq. 'polytrope') then
                do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) - sim_mpoleVX
                            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) - sim_mpoleVY
                            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) - sim_mpoleVZ

                            solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                        enddo
                    enddo
                enddo
            endif

            do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                        if (solnData(DENS_VAR,i,j,k) .lt. damp_dens) then
                            !reduce by constant factor
                            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k)*relax_rate
                            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)*relax_rate
                            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k)*relax_rate
                            solnData(EINT_VAR,i,j,k) = max(solnData(EINT_VAR,i,j,k)*relax_rate, gr_smalle)
                            solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                        endif
                    enddo
                enddo
            enddo

            if (sim_kind .eq. 'cylinder') then
                !mdot = sim_cylinderMdot*8.726109204062576e-21*(-5.648776878000001e15 + &
                !       2.2327568400807356e22*(1.d0/(5.993628868330768e9 + dr_simTime))**0.6666666666666666)**3* &
                !       dlog(3.738175969361188e22/(-5.648776878000001e15 + &
                !       2.2327568400807356e22*(1.d0/(5.993628868330768e9 + dr_simTime))**0.6666666666666666)**1.5d0) 
                mdot = sim_cylinderMdot*1.456480385842967d40/(dexp(9.5414615205053d12/(dr_simTime+sim_tDelay)**(4.d0/3.d0))*&
                       (dr_simTime+sim_tDelay)**(5.d0/3.d0))
                !if (dr_meshMe .eq. MASTER_PE) print *, 'mdot', mdot

                if (mdot .gt. 0.d0) then
                    if (sim_cylinderType .eq. 1) then
                        flow_dist = dsqrt(stvec(1)**2 + stvec(2)**2 + stvec(3)**2)
                        flow_vel = dsqrt(stvec(4)**2 + stvec(5)**2 + stvec(6)**2)
                    else
                        flow_dist = sim_periDist*((dr_simTime+sim_tDelay)/(2.d0*PI*dsqrt(sim_periDist**3.d0/(newton*sim_ptMass))))**(2.d0/3.d0)
                        flow_ecc = 1.d0 - 2.d0*sim_periDist/flow_dist
                        flow_vel = dsqrt((1.d0-flow_ecc)/(1.d0+flow_ecc)*newton*sim_ptMass/(0.5d0*flow_dist))
                    endif

                    polyk = kb*sim_cylinderTemperature/(mp*obj_mu)

                    rho0 = dr_dt/mcs/sim_cylinderNCells*mdot/&
                        (2.d0*PI*flow_dist**3*polyk/newton/sim_ptMass)
                    !rho0 = dr_dt/mcs/sim_cylinderNCells*&
                    !    sim_cylinderMDot*sim_msun/yr/&
                    !    (2.d0*PI*flow_dist**3*polyk/newton/sim_ptMass)

                    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                        if (sim_cylinderType .eq. 1) then
                            zz = zCoord(k) - (sim_zCenter + stvec(3))
                        else
                            zz = zCoord(k) - sim_zCenter
                        endif
                        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                            if (sim_cylinderType .eq. 1) then
                                yy = yCoord(j) - (sim_yCenter + stvec(2))
                            else
                                yy = yCoord(j) - (sim_yCenter + flow_dist)
                            endif
                            do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                                if (sim_cylinderType .eq. 1) then
                                    xx = xCoord(i) - (sim_xCenter + stvec(1))
                                else
                                    xx = xCoord(i) - sim_xCenter
                                endif
                               
                                if (sim_cylinderType .eq. 1) then
                                    dist = dsqrt(xx**2 + zz**2)
                                    perp_dist = yy
                                else
                                    dist = dsqrt(yy**2 + zz**2)
                                    perp_dist = xx
                                endif

                                if (dist .le. sim_cylinderRadius .and. dabs(perp_dist) .le. 0.5d0*sim_cylinderNCells*mcs) then
                                    new_dens = rho0*dexp(-0.5d0*newton*sim_ptMass*dist**2/(flow_dist**3*polyk))

                                    !new_dens = dr_dt/(mcs/flow_vel)*sim_cylinderDensity*&
                                    !    dexp(-(dist/sim_cylinderScale)**2)
                                    solnData(EINT_VAR,i,j,k) = &
                                        (solnData(EINT_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) + &
                                        (eos_gasConstant*sim_cylinderTemperature/((obj_gamc-1.e0)*obj_mu))*&
                                        new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                    if (sim_cylinderType .eq. 1) then
                                        solnData(VELX_VAR,i,j,k) = (solnData(VELX_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) &
                                            + stvec(4)*new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                        solnData(VELY_VAR,i,j,k) = (solnData(VELY_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) &
                                            + stvec(5)*new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                        solnData(VELZ_VAR,i,j,k) = (solnData(VELZ_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) &
                                            + stvec(6)*new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                    else
                                        solnData(VELX_VAR,i,j,k) = (solnData(VELX_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) &
                                            - flow_vel*new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                    endif
                                    solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) + new_dens
                                    ! Not obvious why this needs to be done, but for
                                    ! some reason xn changes in this part of the code.
                                    solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = obj_xn

                                    solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                        0.5d0*(solnData(VELX_VAR,i,j,k)**2.d0 + &
                                               solnData(VELY_VAR,i,j,k)**2.d0 + &
                                               solnData(VELZ_VAR,i,j,k)**2.d0)
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif

            if (sim_kind .eq. 'wind') then
                rho0 = dr_dt*sim_windMdot*sim_msun/yr/((2.d0*PI)**1.5d0*(sim_windKernel*mcs)**3.)

                do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                    zz = zCoord(k) - pvec(3)
                    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                        yy = yCoord(j) - pvec(2)
                        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                            xx = xCoord(i) - pvec(1)
                           
                            dist = dsqrt( xx**2 + yy**2 + zz**2 )

                            if (dist .le. sim_windNCells*mcs) then
                                new_dens = rho0*dexp(-0.5d0*(dist/(sim_windKernel*mcs))**2)

                                !new_dens = dr_dt/(mcs/flow_vel)*sim_cylinderDensity*&
                                !    dexp(-(dist/sim_cylinderScale)**2)
                                solnData(EINT_VAR,i,j,k) = &
                                    (solnData(EINT_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) + &
                                    (eos_gasConstant*sim_windTemperature/((obj_gamc-1.e0)*obj_mu))*&
                                    new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)

                                vvec = (/ xx, yy, zz /)/dist*sim_windVelocity + pvec(4:6)
                                solnData(VELX_VAR,i,j,k) = (solnData(VELX_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) &
                                    + vvec(1)*new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                solnData(VELY_VAR,i,j,k) = (solnData(VELY_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) &
                                    + vvec(2)*new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                solnData(VELZ_VAR,i,j,k) = (solnData(VELZ_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) &
                                    + vvec(3)*new_dens)/(solnData(DENS_VAR,i,j,k) + new_dens)
                                solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) + new_dens
                                ! Not obvious why this needs to be done, but for
                                ! some reason xn changes in this part of the code.
                                solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = obj_xn

                                solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                    0.5d0*(solnData(VELX_VAR,i,j,k)**2.d0 + &
                                           solnData(VELY_VAR,i,j,k)**2.d0 + &
                                           solnData(VELZ_VAR,i,j,k)**2.d0)
                            endif
                        enddo
                    enddo
                enddo
            endif

            call Eos_wrapped(MODE_DENS_EI, blkLimits, lb)
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo

        call Polytrope(blockCount, blockList, dt)
        call Stir(blockCount, blockList, dt) 
        call Flame_step(blockCount, blockList, dt)
        call Burn(blockCount, blockList, dt) 
        call Heat(blockCount, blockList, dt, dr_simTime) 
        call Heatexchange(blockCount, blockList, dt) 
        call Cool(blockCount, blockList, dt, dr_simTime) 
        call Ionize(blockCount, blockList, dt, dr_simTime)
        call EnergyDeposition(blockCount, blockList, dt, dr_simTime)
        call Deleptonize(blockCount, blockList, dt, dr_simTime)
    endif

    return
end subroutine Driver_sourceTerms
