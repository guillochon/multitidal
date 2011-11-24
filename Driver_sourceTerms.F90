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
    use Driver_data, ONLY: dr_simTime
    use Stir_interface, ONLY : Stir
    use Heat_interface, ONLY : Heat
    use Burn_interface, ONLY : Burn
    use Cool_interface, ONLY : Cool
    use Heatexchange_interface, ONLY : Heatexchange
    use EnergyDeposition_interface, ONLY : EnergyDeposition
    use Deleptonize_interface, ONLY : Deleptonize
    use Simulation_data, ONLY: sim_smallX, &
        sim_tRelax, sim_relaxRate, sim_fluffDampCoeff, sim_fluffDampCutoff, sim_accRadius, sim_accCoeff, &
        sim_rhoAmbient, sim_fluidGamma
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize
    use Eos_interface, ONLY : Eos_wrapped, Eos
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
    use Gravity_data, ONLY: grv_factor, grv_ptvec, grv_obvec, grv_ptmass, grv_exactvec
    use Grid_data, ONLY: gr_smalle
    implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    real, intent(IN)    :: dt
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    integer, OPTIONAL, intent(IN) :: pass
    
    integer  ::  i, j, k, l, put, lb, ierr
    real     ::  xx, yy, zz, xcenter, ycenter, zcenter, period
    real     ::  c, avg_dens, avg_temp, avg_pres, dist, vel_sc, mass, y2, z2
    real     ::  tot_avg_dens, tot_avg_temp, tot_avg_pres, tot_mass
    integer  ::  istat, myPE
  
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ
    real, dimension(:,:,:,:),pointer :: solnData
    real, dimension(EOS_NUM) :: eosData
    real, dimension(MDIM) :: avg_vel, tot_avg_vel
    integer,dimension(MDIM) :: axis
  
    logical :: gcell = .true.
  
    integer :: cnt, tot_cnt
    real :: min_cell_size, relax_rate
    
    real :: tinitial, adj_cfactor

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

            do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                        ! Ensure total center of mass has no motion
                        !solnData(VELX_VAR:VELZ_VAR,i,j,k) = solnData(VELX_VAR:VELZ_VAR,i,j,k) - grv_exactvec(4:6)

                        if (solnData(DENS_VAR,i,j,k) .lt. sim_fluffDampCutoff) then
                            !reduce by constant factor
                            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k)*relax_rate
                            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)*relax_rate
                            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k)*relax_rate
                        endif

                        solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                            0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                    enddo
                enddo
            enddo

            if (sim_accRadius .ne. 0.d0) then
                do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                    z2 = (zCoord(k) - (grv_exactvec(3) - grv_obvec(3) + grv_ptvec(3)))**2
                    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                        y2 = z2 + (yCoord(j) - (grv_exactvec(2) - grv_obvec(2) + grv_ptvec(2)))**2
                        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                            dist = dsqrt(y2 + (xCoord(i) - (grv_exactvec(1) - grv_obvec(1) + grv_ptvec(1)))**2)
                            if (dist .le. sim_accRadius) then
                                solnData(DENS_VAR,i,j,k) = max((1.d0-dexp(-(dist/sim_accRadius)**2.d0))*solnData(DENS_VAR,i,j,k), sim_rhoAmbient)
                                solnData(EINT_VAR,i,j,k) = max((1.d0-dexp(-(dist/sim_accRadius)**2.d0))**(sim_fluidGamma - 1.d0)*&
                                    solnData(EINT_VAR,i,j,k), gr_smalle)
                                solnData(VELX_VAR:VELZ_VAR,i,j,k) = (1.d0-dexp(-(dist/sim_accRadius))**2.d0)*solnData(VELX_VAR:VELZ_VAR,i,j,k)
                                solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                    0.5*(solnData(VELX_VAR,i,j,k)**2. + &
                                         solnData(VELY_VAR,i,j,k)**2. + &
                                         solnData(VELZ_VAR,i,j,k)**2.)
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
