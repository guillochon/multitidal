!!****if* source/Simulation/SimulationMain/WDAccrection/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!! ARGUMENTS
!!
!!   myPE -   current processor number
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***

subroutine Simulation_init()

    use Simulation_data 
    use Particles_sinkData, ONLY : particles_local, ipm, ipvx, ipvy, ipvz, iptag, localnp
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Multispecies_interface, ONLY : Multispecies_getSumFrac, Multispecies_getSumInv, Multispecies_getAvg
    use Grid_data, ONLY : gr_globalMe, gr_globalNumProcs, gr_globalComm
    use Logfile_interface, ONLY : Logfile_stampMessage
    use tree, ONLY : lrefine_max
    use Eos_data, ONLY : eos_gasConstant
    use pt_sinkInterface, ONLY : pt_sinkCreateParticle
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use Driver_data, ONLY : dr_simTime, dr_restart
    use IO_data, ONLY : io_scalar

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Multispecies.h"
#include "Eos.h"
#include "Starprof.h"
#include "Particles.h"

    integer             :: ierr, i, nsend, reqr
    double precision  polyk, &
                    x(np),y(np),yp(np), &
                    mass(np),ebind(np), &
                    rhom(np),zbeta(np),ztemp(np),exact(np), &
                    xsurf,ypsurf
    double precision :: cfl, rho0, a, period, start_dist, ecc_anom, newton
    integer mode,iend,pno
    character(len=200) :: logstr
    double precision, dimension(6) :: obvec, ptvec
    integer, dimension(gr_globalNumProcs-1) :: reqs
    integer :: stats(MPI_STATUS_SIZE,gr_globalNumProcs-1)
    integer :: statr(MPI_STATUS_SIZE)

    call PhysicalConstants_get("Newton", newton)

    call RuntimeParameters_get('smallx', sim_smallX)
    call RuntimeParameters_get('smlrho', sim_smallRho)
    call RuntimeParameters_get('smallp', sim_smallP)
    call RuntimeParameters_get("smallt", sim_smallT)
    call RuntimeParameters_get('xmin',sim_xMin)
    call RuntimeParameters_get('ymin',sim_yMin)
    call RuntimeParameters_get('zmin',sim_zMin)
    call RuntimeParameters_get('xmax',sim_xMax)
    call RuntimeParameters_get('ymax',sim_yMax)
    call RuntimeParameters_get('zmax',sim_zMax)
    call RuntimeParameters_get('tinitial',sim_tInitial)

    call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
    call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
    call RuntimeParameters_get('sim_xcenter',sim_xCenter)
    call RuntimeParameters_get('sim_ycenter',sim_yCenter)
    call RuntimeParameters_get('sim_zcenter',sim_zCenter)
    call RuntimeParameters_get('sim_tRelax',sim_tRelax)
    call RuntimeParameters_get('sim_tDelay',sim_tDelay)
    call RuntimeParameters_get('sim_relaxRate',sim_relaxRate)
    call RuntimeParameters_get('sim_softenRadius',sim_softenRadius)
    call RuntimeParameters_get('sim_accRadius',sim_accRadius)
    call RuntimeParameters_get('sim_accCoeff',sim_accCoeff)
    call RuntimeParameters_get('sim_fluffDampCoeff',sim_fluffDampCoeff)
    call RuntimeParameters_get('sim_fluffDampCutoff',sim_fluffDampCutoff)
    call RuntimeParameters_get('sim_fluffRefineCutoff',sim_fluffRefineCutoff)
    call RuntimeParameters_get('sim_maxBlocks',sim_maxBlocks)
    call RuntimeParameters_get("sim_periBeta", sim_periBeta)
    call RuntimeParameters_get("sim_startBeta", sim_startBeta)
    call RuntimeParameters_get("sim_periodFac", sim_periodFac)
    call RuntimeParameters_get("sim_orbEcc", sim_orbEcc)
    call RuntimeParameters_get("sim_useInitialPeakDensity", sim_useInitialPeakDensity)
    call RuntimeParameters_get("sim_ptMassRefine", sim_ptMassRefine)
    call RuntimeParameters_get("sim_ptMassRefRad", sim_ptMassRefRad)
    call RuntimeParameters_get("sim_totForceSub", sim_totForceSub)
    call RuntimeParameters_get("sim_rotFac", sim_rotFac)
    call RuntimeParameters_get("sim_rotAngle", sim_rotAngle)
    call RuntimeParameters_get("sim_tSpinup", sim_tSpinup)
    call RuntimeParameters_get("sim_fixedParticle", sim_fixedParticle)
    call RuntimeParameters_get("sim_xRayFraction", sim_xRayFraction)

    call RuntimeParameters_get("sim_powerLawScale", sim_powerLawScale)
    call RuntimeParameters_get("sim_powerLawExponent", sim_powerLawExponent)
    call RuntimeParameters_get("sim_powerLawExtent", sim_powerLawExtent)
    call RuntimeParameters_get("sim_powerLawMass", sim_powerLawMass)
    call RuntimeParameters_get("sim_powerLawTemperature", sim_powerLawTemperature)

    call RuntimeParameters_get("sim_cylinderScale", sim_cylinderScale)
    call RuntimeParameters_get("sim_cylinderDensity", sim_cylinderDensity)
    call RuntimeParameters_get("sim_cylinderTemperature", sim_cylinderTemperature)
    call RuntimeParameters_get("sim_cylinderRadius", sim_cylinderRadius)
    call RuntimeParameters_get("sim_cylinderMDot", sim_cylinderMDot)
    call RuntimeParameters_get("sim_cylinderNCells", sim_cylinderNCells)
    call RuntimeParameters_get("sim_cylinderType", sim_cylinderType)

    call RuntimeParameters_get("sim_windVelocity", sim_windVelocity)
    call RuntimeParameters_get("sim_windMdot", sim_windMdot)
    call RuntimeParameters_get("sim_windTemperature", sim_windTemperature)
    call RuntimeParameters_get("sim_windNCells", sim_windNCells)
    call RuntimeParameters_get("sim_windLaunchRadius", sim_windLaunchRadius)
    call RuntimeParameters_get("sim_windKernel", sim_windKernel)

    call RuntimeParameters_get("sim_ptMass", sim_ptMass)
    call RuntimeParameters_get("sim_starPtMass", sim_starPtMass)

    if (sim_powerLawExponent .le. -3.d0) then
        call Driver_abortFlash('Error: sim_powerLawExponent must be greater than -3.0')
    endif

    call RuntimeParameters_get('sim_kind',sim_kind)

#ifdef LOADPROFILE
    call RuntimeParameters_get('sim_tAmbient', sim_tAmbient)
    call RuntimeParameters_get('sim_profFile',sim_profFile)

    if (gr_globalMe .eq. MASTER_PE) then
        call read_table_dims(sim_profFile, sim_tableRows, sim_tableCols)
        allocate(sim_table(sim_tableRows, sim_tableCols))
        call read_table(sim_profFile, sim_table, sim_tableRows, sim_tableCols)
        write(*,*) 'r_prof', sim_table(:,R_PROF)
        write(*,*) 'rho_prof', sim_table(:,RHO_PROF)
        write(*,*) 'temp_prof', sim_table(:,TEMP_PROF)
    endif

    call MPI_BCAST(sim_tableRows, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(sim_tableCols, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    if (gr_globalMe .ne. MASTER_PE) allocate(sim_table(sim_tableRows, sim_tableCols))
    call MPI_BCAST(sim_table, sim_tableRows*sim_tableCols, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

    sim_objRadius = sim_table(sim_tableRows,R_PROF)
    sim_objMass = sum(PI*(sim_table(2:sim_tableRows,R_PROF) + sim_table(1:sim_tableRows-1,R_PROF))**2.d0*&
                      (sim_table(2:sim_tableRows,R_PROF) - sim_table(1:sim_tableRows-1,R_PROF))*&
                      sim_table(:,RHO_PROF))
    sim_objCentDens = sim_table(1,RHO_PROF)
#else
    call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
    call RuntimeParameters_get('sim_objPolyN',sim_objPolyN)
    call RuntimeParameters_get('sim_objMass',sim_objMass)
    call RuntimeParameters_get('sim_objCentDens',sim_objCentDens)
    call RuntimeParameters_get('sim_condCoeff',sim_condCoeff)

    obj_xn = 0.d0
    obj_xn(H1_SPEC) = 0.7d0
    obj_xn(HE4_SPEC) = 0.3d0

    call Multispecies_getSumInv(A, obj_mu, obj_xn)
    obj_mu = 1.e0 / obj_mu
    call Multispecies_getSumInv(GAMMA, obj_gamc, obj_xn)
    obj_gamc = 1.e0 / obj_gamc

    if (sim_kind .eq. 'polytrope') then
        mode = 1
        call polytr(sim_objPolyN,sim_objMass/sim_msun,sim_objCentDens,polyk,obj_mu,mode, &
            x,y,yp,obj_radius,obj_rhop,mass,obj_prss,ebind, &
            rhom,ztemp,zbeta,exact,xsurf,ypsurf,np,iend,obj_ipos)
        sim_objRadius = obj_radius(obj_ipos)
    elseif (sim_kind .eq. 'powerlaw') then
        obj_ipos = np
        rho0 = sim_powerLawMass*sim_powerLawExtent**(-3.d0 - sim_powerLawExponent)*&
            sim_powerLawScale**sim_powerLawExponent*(3.d0 + sim_powerLawExponent)/(2.d0*PI)
        do i = 1, np
            obj_radius(i) = sim_powerLawExtent*dble(i)/np
            obj_rhop(i) = rho0*(obj_radius(i)/sim_powerLawScale)**sim_powerLawExponent
            obj_prss(i) = eos_gasConstant*obj_rhop(i) * &
                             sim_powerLawTemperature / obj_mu
        enddo 
        sim_objMass = sim_powerLawMass
        sim_objCentDens = obj_rhop(1)
        sim_objRadius = obj_radius(obj_ipos)
    elseif (sim_kind .eq. 'cylinder') then
        obj_ipos = np
        do i = 1, np
            obj_radius(i) = 5.d0*sim_cylinderScale*dble(i)/np
            obj_rhop(i) = sim_cylinderDensity*dexp((obj_radius(i)/sim_cylinderScale)**2)
            obj_prss(i) = eos_gasConstant*obj_rhop(i) * &
                             sim_cylinderTemperature / obj_mu
        enddo 
        sim_objMass = sim_ptMass*1.d-10
        sim_objCentDens = obj_rhop(1)
        sim_objRadius = obj_radius(obj_ipos)
    elseif (sim_kind .eq. 'wind') then
        sim_objMass = sim_starPtMass
        sim_objRadius = sim_windLaunchRadius
    endif

    !call MPI_BCAST(obj_xn, NSPECIES, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    !call MPI_BCAST(obj_mu, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    !call MPI_BCAST(obj_radius, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    !call MPI_BCAST(obj_rhop, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    !call MPI_BCAST(obj_prss, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    !call MPI_BCAST(obj_ipos, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)              

#endif

    if (gr_globalMe .eq. MASTER_PE) then
        write(logstr, fmt='(A30, ES15.8)') 'Object mass (m_sun):', sim_objMass/sim_msun
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, ES15.8)') 'Object radius:', sim_objRadius
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, ES15.8)') 'Object central density:', sim_objCentDens
        call Logfile_stampMessage(logstr)

        call RuntimeParameters_get("cfl", cfl)

        write(logstr, fmt='(A30, ES15.8)') 'Ambient CFL timestep:', cfl*(sim_xMax - sim_xMin)/NXB/2.d0**(lrefine_max-1)/dsqrt(sim_fluidGamma*sim_pAmbient/sim_rhoAmbient)
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, ES15.8)') 'Fluff CFL timestep:', cfl*(sim_xMax - sim_xMin)/NXB/2.d0**(lrefine_max-1)/dsqrt(sim_fluidGamma*sim_smallP/sim_smallRho)
        call Logfile_stampMessage(logstr)
    endif

    !Sink radius is scaled relative to the original pericenter distance
    
    sim_softenRadius = sim_softenRadius*sim_objRadius/sim_periBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
    sim_accRadius = sim_accRadius*sim_objRadius/sim_periBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
    sim_startDistance = sim_objRadius/sim_startBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)

    write(logstr, fmt='(A30, ES15.8)') 'Sink radius:', sim_softenRadius
    call Logfile_stampMessage(logstr)

    sim_inSubInv = 1.d0/2.d0**sim_nSubZones
    sim_nSubZones = 2**(sim_nSubZones - 1)
    sim_inszd      = (1.d0/sim_nSubZones)**NDIM

    sim_totForceInv = 1.d0/2.d0**sim_totForceSub
    sim_totForceSub = 2**(sim_totForceSub - 1)

    if (sim_periBeta .eq. 0.d0) then
        call RuntimeParameters_get("sim_periDist", sim_periDist)
    else
        sim_periDist = sim_objRadius/sim_periBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
    endif
    if (sim_periodFac .gt. 0.d0) then
        a = sim_periDist/(1.d0 - sim_orbEcc)
        period = dsqrt(4.d0*PI**2.d0/newton/(sim_ptMass + sim_objMass)*a**3.d0)
        sim_periTime = sim_periodFac*period + sim_tRelax
    endif
    if (sim_startBeta .gt. 0.d0) then
        a = sim_periDist/(1.d0 - sim_orbEcc)
        period = dsqrt(4.d0*PI**2.d0/newton/(sim_ptMass + sim_objMass)*a**3.d0)
        start_dist = sim_objRadius/sim_startBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
        if (start_dist .gt. 2.d0*a - sim_periDist) then
            print *, 'Simulation_init', start_dist, 2.d0*a - sim_periDist,&
            sim_objRadius, sim_startBeta, sim_ptMass, sim_objMass
            call Driver_abortFlash('start_dist too large!')
        endif
        ecc_anom = dacos((a - start_dist)/a/sim_orbEcc)
        sim_periTime = abs(ecc_anom - sim_orbEcc*dsin(ecc_anom))*period/2.d0/PI + sim_tRelax
    endif

    sim_fixedPartTag = 0

    stvec = 0.d0
    bhvec = 0.d0

    if (gr_globalMe .eq. MASTER_PE) then
        call calc_orbit(0.d0, sim_objMass, sim_ptMass, obvec, ptvec)
        ptvec = ptvec - obvec
        if (sim_fixedParticle .eq. 1) then
            stvec = -ptvec
        else
            bhvec = ptvec
        endif
        print *, 'stvec', stvec
        print *, 'bhvec', bhvec
    endif

    call MPI_BCAST(bhvec, 6, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    call MPI_BCAST(stvec, 6, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                

    if (dr_restart) then
        call NameValueLL_getInt(io_scalar, "fixedparttag", sim_fixedPartTag, .true., ierr)
        if (ierr /= NORMAL) then
            sim_fixedPartTag = 0
            if (gr_globalMe .eq. MASTER_PE) &
                print *, 'Could not find fixedparttag, searching by mass instead'
            do i=1, localnp
                ! A bit worrisome, should be saving tags to checkpoint.
                if (particles_local(ipm,i) .eq. sim_ptMass .and. sim_fixedParticle .eq. 1) then
                    sim_fixedPartTag = particles_local(iptag,i)
                endif
                if (particles_local(ipm,i) .eq. sim_objMass .and. sim_fixedParticle .eq. 2) then
                    sim_fixedPartTag = particles_local(iptag,i)
                endif
            enddo
            if (sim_fixedPartTag .ne. 0) then
                nsend = 0
                do i = 1, gr_globalNumProcs
                    if (i-1 .eq. gr_globalMe) cycle
                    nsend = nsend + 1
                    call MPI_ISEND(sim_fixedPartTag, 1, FLASH_INTEGER, i-1, 86, gr_globalComm, reqs(nsend), ierr)
                    print *, 'sent to', i-1
                enddo
                print *, gr_globalMe, 'sent all data'
                if (nsend .gt. 0) call MPI_WAITALL(nsend, reqs, stats, ierr)
                print *, gr_globalMe, 'finished waiting'
            else
                print *, 'i am', gr_globalMe
                call MPI_IRECV(sim_fixedPartTag, 1, FLASH_INTEGER, MPI_ANY_SOURCE, 86, gr_globalComm, reqr, ierr)
                call MPI_WAIT(reqr, statr, ierr)
                print *, gr_globalMe, 'received all data'
            endif
        endif
    else
        if (gr_globalMe .eq. MASTER_PE) then
            pno = pt_sinkCreateParticle(sim_xCenter + bhvec(1), sim_yCenter + bhvec(2), &
                sim_zCenter + bhvec(3), 0., 1, gr_globalMe)
            if (sim_fixedParticle .eq. 1) then
                sim_fixedPartTag = particles_local(iptag,pno)
            endif
            particles_local(ipm,pno) = sim_ptMass
            particles_local(ipvx,pno) = bhvec(4)
            particles_local(ipvy,pno) = bhvec(5)
            particles_local(ipvz,pno) = bhvec(6)

            if (sim_kind .ne. 'cylinder') then
                pno = pt_sinkCreateParticle(sim_xCenter + stvec(1), sim_yCenter + stvec(2), &
                    sim_zCenter + stvec(3), 0., 1, gr_globalMe)
                particles_local(ipm,pno) = sim_starPtMass
                particles_local(ipvx,pno) = stvec(4)
                particles_local(ipvy,pno) = stvec(5)
                particles_local(ipvz,pno) = stvec(6)
                if (sim_fixedParticle .eq. 2) then
                    sim_fixedPartTag = particles_local(iptag,pno)
                endif
                print *, "Fixed particle's tag: ", sim_fixedPartTag
            endif
        endif
        call MPI_BCAST(sim_fixedPartTag, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    endif

    if (gr_globalMe .eq. MASTER_PE) then
        write(logstr, fmt='(A30, 2ES15.8)') 'Start distance:', start_dist
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 6ES15.8)') 'Pt. mass start pos:', ptvec
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 2ES15.8)') 'Semimajor axis, eccentricity:', a, sim_orbEcc
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 2ES15.8)') 'Period, pericenter time:', period, sim_periTime
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 2ES15.8)') 'Obj. radius, pericenter dist:', sim_objRadius, sim_periDist
        call Logfile_stampMessage(logstr)

        if (dr_simTime .eq. sim_tInitial) then
            open(unit = 11, file = 'extras.dat', status = 'unknown')
            write(11, fmt='(ES22.15)') sim_periBeta
            write(11, fmt='(ES22.15)') period
            write(11, fmt='(ES22.15)') sim_periTime
            write(11, fmt='(ES22.15)') sim_objRadius
            write(11, fmt='(ES22.15)') sim_objMass
            write(11, fmt='(ES22.15)') sim_orbEcc
            write(11, fmt='(ES22.15)') sim_ptMass
            write(logstr, fmt='(I4)')  lrefine_max
            write(11, fmt='(A)') adjustl(logstr)
            write(11, fmt='(ES22.15)') sim_xMax - sim_xMin
            close(11) 
        endif
    endif
end subroutine Simulation_init
