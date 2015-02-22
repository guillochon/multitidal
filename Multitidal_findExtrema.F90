!!****if* source/Grid/GridSolvers/Multipole/Multitidal_findExtrema
!!
!! NAME
!!
!!  Multitidal_findExtrema()
!!
!! SYNOPSIS
!!
!!  Multitidal_findExtrema(integer,intent(in) : blockID,
!!                     integer,intent(in) : ivar,
!!                     integer,intent(in) : flag,
!!                     real,intent (out)  : extrema)
!!
!! DESCRIPTION
!!
!! Find maximum or minimum value for a given variable.
!!
!!***


subroutine Multitidal_findExtrema (ivar,flag,extrema)
  
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr,&
                             Grid_getBlkBoundBox,&
                             Grid_getBlkIndexLimits,&
                             Grid_getListOfBlocks
  use Grid_data, ONLY : gr_meshComm

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: ivar, flag
  real,intent(INOUT)  :: extrema
  integer :: blockCount
  integer, dimension(MAXBLOCKS) :: blockList
  real :: extrema_local
  
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  real, pointer, dimension(:,:,:,:)      :: solnData
  integer            :: i, j, k, imax, jmax, kmax, imin, jmin, kmin
  integer            :: blockID, lb, ierr
  
  !=====================================================================
  
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  if (flag .eq. 1) then
      extrema_local = -huge(0.d0)
  else
      extrema_local = huge(0.d0)
  endif

  do lb = 1, blockCount
     blockID = blockList(lb)
     call Grid_getBlkBoundBox(blockID, bndBox)
     call Grid_getBlkPtr(blockID, solnData)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     kmax = blkLimits(HIGH,KAXIS)
     kmin = blkLimits(LOW,KAXIS)  
     jmax = blkLimits(HIGH,JAXIS)
     jmin = blkLimits(LOW,JAXIS)
     imax = blkLimits(HIGH,IAXIS)
     imin = blkLimits(LOW,IAXIS)

     do k = kmin, kmax
        do j = jmin, jmax
           do i = imin, imax
               if (flag .eq. 1) then
                   if (solnData(ivar,i,j,k) .gt. extrema_local) extrema_local = solnData(ivar,i,j,k)
               else
                   if (solnData(ivar,i,j,k) .lt. extrema_local) extrema_local = solnData(ivar,i,j,k)
               endif
           enddo
        enddo
     enddo
  enddo

  call MPI_ALLREDUCE(extrema_local, extrema, 1, FLASH_REAL, MPI_MAX, gr_meshComm, ierr)
  
  !==========================================================================
  call Grid_releaseBlkPtr(blockID, solnData)
  
  return
end subroutine Multitidal_findExtrema
