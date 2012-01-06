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
    use Driver_data, ONLY: dr_simTime, dr_meshComm
    use Stir_interface, ONLY : Stir
    use Heat_interface, ONLY : Heat
    use Burn_interface, ONLY : Burn
    use Cool_interface, ONLY : Cool
    use Heatexchange_interface, ONLY : Heatexchange
    use EnergyDeposition_interface, ONLY : EnergyDeposition
    use Deleptonize_interface, ONLY : Deleptonize
    use Simulation_data, ONLY: sim_smallX, &
        sim_tRelax, sim_relaxRate, sim_fluffDampCoeff, sim_fluffDampCutoff, sim_accRadius, sim_accCoeff, &
        sim_fluidGamma
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize, Grid_findExtrema
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
    use Gravity_data, ONLY: grv_ptvec, grv_obvec, grv_ptmass, grv_exactvec, grv_optmass, grv_momacc, &
        grv_angmomacc, grv_eneracc, grv_massacc, grv_comPeakCut, grv_totmass, grv_totmass
    use Grid_data, ONLY: gr_smalle, gr_meshMe
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
    real    :: tinitial, vol, ldenscut, denscut, extrema
  
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
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

    call RuntimeParameters_get('tinitial',tinitial)
    if (dr_simTime .lt. tinitial + sim_tRelax) then
        relax_rate = (dr_simTime - tinitial)/sim_tRelax*(1.0 - sim_relaxRate) + sim_relaxRate
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
                        solnData(VELX_VAR:VELZ_VAR,i,j,k) = solnData(VELX_VAR:VELZ_VAR,i,j,k) - grv_exactvec(4:6)
                        !reduce by constant factor

                        solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k)*relax_rate
                        solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)*relax_rate
                        solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k)*relax_rate

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
    else
        pt_pos = grv_exactvec - grv_obvec + grv_ptvec !Used for calculating properties of mass accreted onto point mass
        relax_rate = sim_fluffDampCoeff

        ldenscut = -huge(0.d0)
        do lb = 1,blockCount
           call Grid_findExtrema (blockList(lb), DENS_VAR, 1, extrema)
           if (extrema .gt. ldenscut) ldenscut = extrema
        enddo

        call MPI_ALLREDUCE(ldenscut, denscut, 1, FLASH_REAL, & 
                           MPI_MAX, dr_meshComm, ierr)

        denscut = denscut*grv_comPeakCut

        ! Calculate the average velocity of the peak to subtract it
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

            vol = (xCoord(2) - xCoord(1))**3.d0
            do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                        if (solnData(DENS_VAR,i,j,k) .lt. sim_fluffDampCutoff) cycle
                        tot_mass = tot_mass + vol*solnData(DENS_VAR,i,j,k)
                        tot_mom = tot_mom + vol*solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR:VELZ_VAR,i,j,k)
                        if (solnData(DENS_VAR,i,j,k) .lt. denscut) cycle
                        peak_mass = peak_mass + vol*solnData(DENS_VAR,i,j,k)
                        peak_mom = peak_mom + vol*solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR:VELZ_VAR,i,j,k)
                    enddo
                enddo
            enddo

            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo

        call MPI_ALLREDUCE(tot_mass, gtot_mass, 1, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        call MPI_ALLREDUCE(tot_mom, gtot_mom, 3, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        call MPI_ALLREDUCE(peak_mass, gpeak_mass, 1, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        call MPI_ALLREDUCE(peak_mom, gpeak_mom, 3, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        
        tot_avg_vel = gtot_mom/gtot_mass
        peak_avg_vel = gpeak_mom/gpeak_mass

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
                        if (solnData(DENS_VAR,i,j,k) .lt. sim_fluffDampCutoff) then
                            !reduce by constant factor
                            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k)*relax_rate
                            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)*relax_rate
                            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k)*relax_rate
                        else
                            ! Ensure total center of mass has no motion
                            solnData(VELX_VAR:VELZ_VAR,i,j,k) = solnData(VELX_VAR:VELZ_VAR,i,j,k) - peak_avg_vel
                        endif

                        solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                            0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                    enddo
                enddo
            enddo

            !print *, 'before', grv_obvec(4:6)
            !print *, 'peak_avg_vel', peak_avg_vel
            !print *, 'after', grv_obvec(4:6)
            !call Driver_abortFlash('done')

            if (sim_accRadius .ne. 0.d0) then
                do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                    zz = zCoord(k)
                    z2 = (zz - (grv_exactvec(3) - grv_obvec(3) + grv_ptvec(3)))**2
                    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                        yy = yCoord(j)
                        y2 = z2 + (yy - (grv_exactvec(2) - grv_obvec(2) + grv_ptvec(2)))**2
                        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                            xx = xCoord(i)
                            dist = dsqrt(y2 + (xx - (grv_exactvec(1) - grv_obvec(1) + grv_ptvec(1)))**2)
                            if (dist .le. sim_accRadius) then
                                vol = (xCoord(2) - xCoord(1))**3.d0
                                mass_acc = vol*(solnData(DENS_VAR,i,j,k) - max(dexp(-((sim_accRadius-dist)/sim_accRadius)**2.d0)*solnData(DENS_VAR,i,j,k), sim_fluffDampCutoff))
                                tot_com_acc = tot_com_acc + mass_acc*((/ xx, yy, zz /) - pt_pos(1:3))
                                tot_mom_acc = tot_mom_acc + mass_acc*(solnData(VELX_VAR:VELZ_VAR,i,j,k) - pt_pos(4:6))
                                tot_angmom_acc = tot_angmom_acc + mass_acc*&
                                    (/ (yy-pt_pos(2))*(solnData(VELZ_VAR,i,j,k)-pt_pos(6)) - &
                                        (zz-pt_pos(3))*(solnData(VELY_VAR,i,j,k)-pt_pos(5)), &
                                       (zz-pt_pos(3))*(solnData(VELX_VAR,i,j,k)-pt_pos(4)) - &
                                        (xx-pt_pos(1))*(solnData(VELZ_VAR,i,j,k)-pt_pos(6)), &
                                       (xx-pt_pos(1))*(solnData(VELY_VAR,i,j,k)-pt_pos(5)) - &
                                        (yy-pt_pos(2))*(solnData(VELX_VAR,i,j,k)-pt_pos(4)) /)
                                tot_mass_acc = tot_mass_acc + mass_acc
                                tot_ener_acc = tot_ener_acc + mass_acc*(solnData(EINT_VAR,i,j,k) + &
                                    0.5d0*sum((solnData(VELX_VAR:VELZ_VAR,i,j,k) - pt_pos(4:6))**2.d0))
                                solnData(DENS_VAR,i,j,k) = max(dexp(-((sim_accRadius-dist)/sim_accRadius)**2.d0)*solnData(DENS_VAR,i,j,k), sim_fluffDampCutoff)
                                solnData(EINT_VAR,i,j,k) = max(dexp(-((sim_accRadius-dist)/sim_accRadius)**2.d0)**(sim_fluidGamma - 1.d0)*&
                                    solnData(EINT_VAR,i,j,k), gr_smalle)
                                !solnData(VELX_VAR:VELZ_VAR,i,j,k) = dexp(-((sim_accRadius-dist)/sim_accRadius)**2.d0)*solnData(VELX_VAR:VELZ_VAR,i,j,k)
                                solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                    0.5d0*(solnData(VELX_VAR,i,j,k)**2.d0 + &
                                           solnData(VELY_VAR,i,j,k)**2.d0 + &
                                           solnData(VELZ_VAR,i,j,k)**2.d0)
                            endif
                        enddo
                    enddo
                enddo
            endif
  
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo

        ! Add subtracted velocity to tracking point
        !grv_obvec(4:6) = grv_obvec(4:6) + peak_avg_vel

        call MPI_ALLREDUCE(tot_mass_acc, gtot_mass_acc, 1, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        call MPI_ALLREDUCE(tot_ener_acc, gtot_ener_acc, 1, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        call MPI_ALLREDUCE(tot_com_acc, gtot_com_acc, 3, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        call MPI_ALLREDUCE(tot_mom_acc, gtot_mom_acc, 3, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)
        call MPI_ALLREDUCE(tot_angmom_acc, gtot_angmom_acc, 3, FLASH_REAL, MPI_SUM, dr_meshComm, ierr)

        !if (gr_meshMe .eq. MASTER_PE) then
        !    print *, 'Mass accreted this step: ', gtot_mass_acc
        !    print *, 'COM accreted mass: ', gtot_com_acc/gtot_mass_acc
        !    print *, 'Mom accreted mass: ', gtot_mom_acc/gtot_mass_acc
        !    print *, 'Net shift: ', grv_ptvec(1:3) - (grv_ptmass*grv_ptvec(1:3) + &
        !        (gtot_com_acc/gtot_mass_acc - grv_exactvec(1:3) + grv_obvec(1:3) - grv_ptvec(1:3))*gtot_mass_acc) / (grv_ptmass + gtot_mass_acc)
        !    print *, 'Net vel: ', grv_ptvec(4:6) - (grv_ptmass*grv_ptvec(4:6) + &
        !        (gtot_mom_acc/gtot_mass_acc - grv_exactvec(4:6) + grv_obvec(4:6) - grv_ptvec(4:6))*gtot_mass_acc) / (grv_ptmass + gtot_mass_acc)
        !endif

        grv_ptvec(1:3) = (grv_ptmass*grv_ptvec(1:3) + gtot_com_acc) / (grv_ptmass + gtot_mass_acc)
        grv_ptvec(4:6) = (grv_ptmass*grv_ptvec(4:6) + gtot_mom_acc) / (grv_ptmass + gtot_mass_acc)
        grv_obvec(1:3) = (grv_totmass*grv_obvec(1:3) - gtot_com_acc + &
            gtot_mass_acc*(grv_exactvec(1:3) - pt_pos(1:3))) / (grv_totmass - gtot_mass_acc)
        grv_obvec(4:6) = (grv_totmass*grv_obvec(4:6) - gtot_mom_acc + &
            gtot_mass_acc*(grv_exactvec(4:6) - pt_pos(4:6))) / (grv_totmass - gtot_mass_acc)
        grv_optmass = grv_ptmass
        grv_ptmass = grv_ptmass + gtot_mass_acc
        grv_massacc = grv_massacc + gtot_mass_acc
        grv_momacc = grv_momacc + gtot_mom_acc
        grv_angmomacc = grv_angmomacc + gtot_angmom_acc
        grv_eneracc = grv_eneracc + gtot_ener_acc

        call Stir(blockCount, blockList, dt) 
        call Burn(blockCount, blockList, dt) 
        call Heatexchange(blockCount, blockList, dt) 
        call EnergyDeposition(blockCount, blockList, dt, dr_simTime)
        call Deleptonize(blockCount, blockList, dt, dr_simTime)
        call Heat(blockCount, blockList, dt, dr_simTime) 
        call Cool(blockCount, blockList, dt, dr_simTime) 
    endif
  
    return
end subroutine Driver_sourceTerms
