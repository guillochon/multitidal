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



subroutine Driver_sourceTerms(blockCount, blockList, dt) 

    use Driver_data, ONLY: dr_simTime
    use Stir_interface, ONLY : Stir
    use Heat_interface, ONLY : Heat
    use Burn_interface, ONLY : Burn
    use Cool_interface, ONLY : Cool
    use Simulation_data, ONLY: sim_smallX, sim_table, sim_tableRows, sim_tableCols, &
        sim_tRelax, sim_relaxRate, sim_BHRadius
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize, Grid_getMyPE
    use Eos_interface, ONLY : Eos_wrapped, Eos
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, Mtot
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
    use Gravity_data, ONLY: grv_factor, grv_ptxpos, grv_ptypos, grv_ptzpos
    implicit none
#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  
    real, intent(IN)    :: dt
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    
    integer  ::  i, j, k, l, put, lb, ierr
    real     ::  xx, yy, zz, xcenter, ycenter, zcenter
    real     ::  c, avg_dens, avg_temp, avg_pres, dist, vel_sc
    integer  ::  istat, myPE
  
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ
    real, dimension(:,:,:,:),pointer :: solnData
    real, dimension(:,:,:,:),allocatable :: velCopy
    real, dimension(EOS_NUM) :: eosData
    real, dimension(MDIM) :: avg_vel
    integer,dimension(MDIM) :: axis
  
    logical :: gcell = .true.
  
    integer cnt
    real :: min_cell_size, relax_rate
    
    real :: tinitial, adj_cfactor, newx, newy, newz
  
    call RuntimeParameters_get('tinitial',tinitial)
    call gr_mpoleCenterOfMass(DENS_VAR)
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
            allocate(velCopy(VELX_VAR:VELZ_VAR,&
                             blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS), &
                             blkLimits(LOW, JAXIS):blkLimits(HIGH, JAXIS), &
                             blkLimits(LOW, KAXIS):blkLimits(HIGH, KAXIS)),stat=istat)
            velCopy(VELX_VAR:VELZ_VAR,:,:,:) = solnData(VELX_VAR:VELZ_VAR,blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS),&
                                                                          blkLimits(LOW, JAXIS):blkLimits(HIGH, JAXIS),&
                                                                          blkLimits(LOW, KAXIS):blkLimits(HIGH, KAXIS))
            do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                        !Arnett 1994
                        !velCopy(VELX_VAR,i,j,k) = 0.5*solnData(VELX_VAR,i,j,k) + 0.25*solnData(VELX_VAR,i-1,j,k) &
                        !    + 0.25*solnData(VELX_VAR,i+1,j,k)
                        !velCopy(VELY_VAR,i,j,k) = 0.5*solnData(VELY_VAR,i,j,k) + 0.25*solnData(VELY_VAR,i,j-1,k) &
                        !    + 0.25*solnData(VELY_VAR,i,j+1,k)
                        !velCopy(VELZ_VAR,i,j,k) = 0.5*solnData(VELZ_VAR,i,j,k) + 0.25*solnData(VELZ_VAR,i,j,k-1) &
                        !    + 0.25*solnData(VELZ_VAR,i,j,k+1)

                        !reduce by constant factor
                        velCopy(VELX_VAR,i,j,k) = velCopy(VELX_VAR,i,j,k)*relax_rate
                        velCopy(VELY_VAR,i,j,k) = velCopy(VELY_VAR,i,j,k)*relax_rate
                        velCopy(VELZ_VAR,i,j,k) = velCopy(VELZ_VAR,i,j,k)*relax_rate
                        
                        eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
                        eosData(EOS_TEMP) = solnData(TEMP_VAR,i,j,k)
                        call Eos(MODE_DENS_TEMP,1,eosData,solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
                        solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
                        solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
                        solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
                        solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                            0.5*(velCopy(VELX_VAR,i,j,k)**2. + velCopy(VELY_VAR,i,j,k)**2. + velCopy(VELZ_VAR,i,j,k)**2.)
                        !solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                        !    0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                    enddo
                enddo
            enddo
  
            solnData(VELX_VAR:VELZ_VAR,blkLimits(LOW, IAXIS):blkLimits(HIGH, IAXIS),&
                                       blkLimits(LOW, JAXIS):blkLimits(HIGH, JAXIS),&
                                       blkLimits(LOW, KAXIS):blkLimits(HIGH, KAXIS)) = velCopy(VELX_VAR:VELZ_VAR,:,:,:)
            deallocate(velCopy)
  
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
            call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockList(lb))
        enddo
    else
        call RuntimeParameters_get('sim_xctr',xcenter)
        call RuntimeParameters_get('sim_yctr',ycenter)
        call RuntimeParameters_get('sim_zctr',zcenter)
        call Grid_getMinCellSize(min_cell_size)
        call gr_mpoleCenterOfMass(DENS_VAR)
        call PhysicalConstants_get('speed of light',c)

        call parabolic_orbit(dr_simTime, newx, newy, newz)
        newx = grv_ptxpos + newx
        newy = grv_ptypos + newy
        newz = grv_ptzpos + newz

        ! calculate the average properties near the black hole and reset the BH to these averages.
        avg_dens = 0.0
        avg_temp = 0.0
        avg_pres = 0.0
        cnt = 0
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
            do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                zz = zCoord(k)
                do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
                    yy = yCoord(j)
                    do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
                        xx = xCoord(i)
                        dist = sqrt((xx-newx)**2.+(yy-newy)**2.+(zz-newz)**2.)
                        if (dist .gt. abs(sim_BHRadius*2.*grv_factor/c**2.) .and. dist .lt. abs(sim_BHRadius*2.*grv_factor/c**2.) + min_cell_size) then
                            avg_dens = avg_dens + solnData(DENS_VAR,i,j,k)
                            avg_temp = avg_temp + solnData(TEMP_VAR,i,j,k)
                            avg_pres = avg_pres + solnData(PRES_VAR,i,j,k)
                            cnt = cnt + 1
                        endif
                    enddo
                enddo
            enddo
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo
  
        call MPI_ALLREDUCE(cnt, cnt, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(avg_dens, avg_dens, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(avg_temp, avg_temp, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(avg_pres, avg_pres, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
        avg_dens = avg_dens / cnt
        avg_temp = avg_temp / cnt
        avg_pres = avg_pres / cnt
        
        if (cnt .gt. 0) then
            call Grid_getMyPE(myPE)
            if (myPE .eq. MASTER_PE) then
                write(*,*) 'avg_dens, avg_pres, avg_temp'
                write(*,*) avg_dens, avg_pres, avg_temp
            endif

            xn = sim_smallX
            xn(H1_SPEC) = 1.0
            eosData(EOS_DENS) = avg_dens
            eosData(EOS_TEMP) = avg_temp
            eosData(EOS_PRES) = avg_pres
            call Eos(MODE_DENS_PRES,1,eosData,xn)
        endif

        !! calculate the velocity of the center of the star so it can be substracted from all grid cells.
        !avg_vel = 0
        !cnt = 0
        !do lb = 1, blockCount
        !    call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
        !    sizeX = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
        !    allocate(xCoord(sizeX),stat=istat)
        !    sizeY = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
        !    allocate(yCoord(sizeY),stat=istat)
        !    sizeZ = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
        !    allocate(zCoord(sizeZ),stat=istat)
  
        !    if (NDIM == 3) call Grid_getCellCoords&
        !                        (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
        !    if (NDIM >= 2) call Grid_getCellCoords&
        !                        (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
        !    call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
  
        !    call Grid_getBlkPtr(blockList(lb),solnData)
        !    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        !        zz = zCoord(k)
        !        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        !            yy = yCoord(j)
        !            do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
        !                xx = xCoord(i)
        !                if (sqrt((xx-Xcm)**2.+(yy-Ycm)**2.+(zz-Zcm)**2.) .le. 5*min_cell_size) then
        !                    avg_vel(1) = avg_vel(1) + solnData(VELX_VAR,i,j,k)
        !                    avg_vel(2) = avg_vel(2) + solnData(VELY_VAR,i,j,k)
        !                    avg_vel(3) = avg_vel(3) + solnData(VELZ_VAR,i,j,k)
        !                    cnt = cnt + 1
        !                endif
        !            enddo
        !        enddo
        !    enddo
        !    call Grid_releaseBlkPtr(blockList(lb), solnData)
        !    deallocate(xCoord)
        !    deallocate(yCoord)
        !    deallocate(zCoord)
        !enddo
  
        !call MPI_ALLREDUCE(cnt, cnt, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        !call MPI_ALLREDUCE(avg_vel, avg_vel, MDIM, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
        !avg_vel = avg_vel / cnt

        do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
  
            call Grid_getBlkPtr(blockList(lb),solnData)
            do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
               do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
                  do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
                     !Reset region near black hole (within a few Schwarzschild radii)
                     if (sqrt((xx-newx)**2.+(yy-newy)**2.+(zz-newz)**2.) .lt. abs(sim_BHRadius*2.*grv_factor/c**2.)) then
                         solnData(VELX_VAR,i,j,k) = 0.0d0
                         solnData(VELY_VAR,i,j,k) = 0.0d0
                         solnData(VELZ_VAR,i,j,k) = 0.0d0
                         solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = xn
                         solnData(DENS_VAR,i,j,k) = eosData(EOS_DENS)
                         solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
                         solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
                         solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
                         solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
                         solnData(GAME_VAR,i,j,k) = eosData(EOS_PRES)/&
                              (eosData(EOS_EINT)*eosData(EOS_DENS)) + 1.0
                     endif
                     !!Subtract the COM's velocity from all grid cells (keeps star in center of grid)
                     !solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) - avg_vel(1)
                     !solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) - avg_vel(2)
                     !solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) - avg_vel(3)
                     !!Add a fictitious force to move the center of mass to the center of the simulation
                     !vel_sc = min(sqrt(avg_vel(1)**2. + avg_vel(2)**2. + avg_vel(3)**2.) / &
                     !   sqrt((xcenter - Xcm)**2. + (ycenter - Ycm)**2. + (zcenter - Zcm)**2.) * dt, 1.0)
                     !solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + (xcenter - Xcm)/dt * vel_sc
                     !solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + (ycenter - Ycm)/dt * vel_sc
                     !solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) + (zcenter - Zcm)/dt * vel_sc
                     solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                         0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                  enddo
               enddo
            enddo
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockList(lb))
        enddo
  
        call Stir(blockCount, blockList, dt) 
        call Burn(blockCount, blockList, dt) 
        call Heat(blockCount, blockList, dt, dr_simTime) 
        call Cool(blockCount, blockList, dt, dr_simTime) 
    endif
  
    return
end subroutine Driver_sourceTerms
