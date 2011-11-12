!!****if* source/Driver/DriverMain/Driver_verifyInitDt
!!
!! NAME
!!  Driver_verifyInitDt
!!
!! SYNOPSIS
!!  Driver_verifyInitDt(integer (IN) :: myPE)
!!
!! DESCRIPTION
!!
!! The initial timestep "dt" is a runtime parameter for the simulations.
!! This routine makes sure that users haven't inadvertently provided
!! the initial value for dt that violates the Courant-Friedrichs-Levy
!! condition.
!!
!! ARGUMENTS
!!
!!  myPE - current processor
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_myPE or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_verifyInitDt(myPE)

  use Driver_data, ONLY : dr_restart, dr_dt, dr_dtInit, dr_dtOld
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getDeltas, &
    Grid_getBlkPtr, Grid_releaseBlkPtr
  use Hydro_interface, ONLY : Hydro_computeDt
  use Cosmology_interface, ONLY: Cosmology_computeDt
  implicit none       

#include "Flash.h"
#include "constants.h" 
  include "Flash_mpi.h"
  
  integer, intent(in) :: myPE
  real :: dtCheck  ,dtCFL
  integer :: localNumBlocks

  integer    :: dtMinLoc(5)
  integer :: i, ierr
  integer, dimension(MAXBLOCKS) :: blockList

  integer :: coordSize
  logical :: gcell = .true.
  real, dimension(MDIM) :: del


#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: xCoord, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) :: yCoord, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) :: zCoord, dz, uzgrid
#else
  real, allocatable,dimension(:)::&
       xCoord,dx,uxgrid,yCoord,dy,uygrid,zCoord,dz,uzgrid
#endif

  !arrays which hold the starting and ending indicies of a block
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC

  !!coordinate infomration to be passed into physics  
  real, pointer :: solnData(:,:,:,:)
  integer :: isize,jsize,ksize


  if (.not. dr_restart) then
     ! compute the CFL timestep for the simulation and compare it to the
     ! user specified initial timestep.  Scream loudly if there is a problem.

     !initialize values 
     dtCheck = huge(dtCheck)
     dtMinLoc(:) = 0
     
     call Grid_getListOfBlocks(LEAF,blockList,localNumBlocks)

!!     call Grid_fillGuardCells(myPE,CENTER,ALLDIR)

     do i = 1, localNumBlocks
        
        !There is some overhead in calling Hydro_computeDt.  Although it is a
        !pain to get the coordinates and solution data before calling the 
        !routine, this is just initialization.  Getting the coordinates inside
        !Hydro_computeDt would be much more costly during the run
        
        
        !!Get the coordinate information for all the
        call Grid_getBlkIndexLimits(blockList(i),blkLimits,blkLimitsGC)
        isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
        
#ifndef FIXEDBLOCKSIZE
        allocate(xCoord(isize))
        allocate(dx(isize))
        allocate(uxgrid(isize))
        allocate(yCoord(jsize))
        allocate(dy(jsize))
        allocate(uygrid(jsize))
        allocate(zCoord(ksize))
        allocate(dz(ksize))
        allocate(uzgrid(ksize))
#endif

!print*,'before coords in Driver_init'        
        coordSize = isize
        call Grid_getCellCoords(IAXIS,blockList(i),CENTER,gcell,xCoord,coordSize)

        coordSize = jsize
        call Grid_getCellCoords(JAXIS,blockList(i),CENTER,gcell,yCoord,coordSize)

        coordSize = ksize
        call Grid_getCellCoords(KAXIS,blockList(i),CENTER,gcell,zCoord,coordSize)

!print*,'after coords in Driver_init'
        
        call Grid_getDeltas(blockList(i), del)
        dx(:) = del(1)
        dy(:) = del(2)
        dz(:) = del(3)
        
        uxgrid(:) = 0
        uygrid(:) = 0
        uzgrid(:) = 0
        
        call Grid_getBlkPtr(blockList(i),solnData)

        call Hydro_computeDt ( blockList(i), myPE, &
             xCoord, dx, uxgrid, &
             yCoord, dy, uygrid, &
             zCoord, dz, uzgrid, &
             blkLimits,blkLimitsGC,  &
             solnData,      &
             dtCheck, dtMinLoc)

        call Grid_releaseBlkPtr(blockList(i),solnData)

#ifndef FIXEDBLOCKSIZE
        deallocate(xCoord)
        deallocate(dx)
        deallocate(uxgrid)
        deallocate(yCoord)
        deallocate(dy)
        deallocate(uygrid)
        deallocate(zCoord)
        deallocate(dz)
        deallocate(uzgrid)
#endif

     end do

     ! find the minimum across all processors, store it in dtCFL on MasterPE
     call MPI_AllReduce(dtCheck, dtCFL, 1, FLASH_REAL, MPI_MIN, &
          MPI_COMM_WORLD, ierr)

#define TIMESTEP_SLOW_START_FACTOR 0.1
     if (dr_dtInit > TIMESTEP_SLOW_START_FACTOR*dtCFL) then
        
        if (MyPE .EQ. MASTER_PE) then
           print *, '***********************************************************'
           print *, ' Warning: The initial timestep is too large.'
           print *, '   initial timestep = ', dr_dtInit
           print *, '   CFL timestep     = ', dtCFL
           print *, ' Resetting dtinit to TIMESTEP_SLOW_START_FACTOR*dtcfl.'
           print *, '***********************************************************'
           print *, ' '
        endif
        
        dr_dt = TIMESTEP_SLOW_START_FACTOR*dtCFL
        
     else
        
        dr_dt = dr_dtInit
        
     endif

     dr_dtOld = dr_dt
     !print *, dr_dt, "dr_dt initial final"
     call Cosmology_computeDt(dtCheck)

  endif
     
  return
end subroutine Driver_verifyInitDt








