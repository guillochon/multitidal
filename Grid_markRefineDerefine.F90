!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine(integer(IN) :: myPE)
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!!  myPE : my processor id.
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_myPE or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine(MyPE)

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,gr_refine_val_cutoff,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_blkList, gr_lrefineLimit
  use tree, ONLY : newchild, refine, derefine, stay, nodetype,&
       lrefine,lrefine_max, parent, nchild,child
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getListOfBlocks,&
                             Grid_releaseBlkPtr,Grid_fillGuardCells
  use Driver_interface, ONLY: Driver_getSimTime
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_data, ONLY: sim_starRadius

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  integer,intent(IN) :: MyPE
  
  real, dimension(:,:,:,:), pointer :: solnData
  real :: ref_cut,deref_cut,ref_filter,ref_val_cut,t
  real :: xcenter,ycenter,zcenter
  real :: maxptime,shrinkdur,box,zbox
  integer       :: l,iref,blkCount,lb,i,j
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask
  real :: maxvals(MAXBLOCKS),maxvals_parent(MAXBLOCKS)
  integer :: nsend,nrecv,ierr
  integer :: reqr(MAXBLOCKS),reqs(MAXBLOCKS)
  integer :: statr(MPI_STATUS_SIZE,MAXBLOCKS)
  integer :: stats(MPI_STATUS_SIZE,MAXBLOCKS)
  integer :: max_refine

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  call Grid_fillGuardCells(MyPE,CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  !do l = 1,gr_numRefineVars
  !   iref = gr_refine_var(l)
  !   ref_cut = gr_refine_cutoff(l)
  !   deref_cut = gr_derefine_cutoff(l)
  !   ref_filter = gr_refine_filter(l)
  !   ref_val_cut = gr_refine_val_cutoff(l)
  !   call gr_markRefineDerefine(MyPE,iref,ref_cut,deref_cut,ref_filter,&
  !       ref_val_cut)
  !end do

  call Driver_getSimTime(t)
  call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList,blkCount)

  if (gr_lrefineLimit .ne. 0) then
      max_refine = gr_lrefineLimit
  else
      max_refine = lrefine_max
  endif

  call RuntimeParameters_get('xmax',xcenter)      
  call RuntimeParameters_get('ymax',ycenter)      
  call RuntimeParameters_get('zmax',zcenter)      
  xcenter = xcenter / 2.
  ycenter = ycenter / 2.
  zcenter = zcenter / 2.

  if (t .eq. 0.0) then
      !write(*,*) 'entered sphere refine'
      call gr_markInRadius(xcenter,ycenter,zcenter,1.1*sim_starRadius,max_refine,0)
  else
      do l = 1,gr_numRefineVars
         ref_val_cut = gr_refine_val_cutoff(l)
         call gr_markVarThreshold(iref,ref_val_cut,0,max_refine)
      enddo

      !call gr_markInRectangle(xcenter-box,xcenter+box,ycenter-box,ycenter+box,&
      !   zcenter-zbox,zcenter+zbox,8,1,1)

      do i = 1, blkCount
         lb = gr_blkList(i)
         call Grid_getBlkPtr(lb,solnData,CENTER)
         maxvals(lb) = maxval(solnData(iref,:,:,:))
         call Grid_releaseBlkPtr(lb,solnData)
      end do

!     Communicate maxvals of parents to their leaf children.
!     Maximally refined children collect messages from parents.

      maxvals_parent(:) = 0.0
      nrecv = 0
      do i = 1, blkCount
         lb = gr_blkList(i)
         if (nodetype(lb) == LEAF .AND. lrefine(lb) == lrefine_max) then
            if(parent(1,lb).gt.-1) then
               if (parent(2,lb).ne.MyPE) then
                  nrecv = nrecv + 1
                  call MPI_IRecv(maxvals_parent(lb),1, &
                     FLASH_REAL, &
                     parent(2,lb), &
                     lb, &
                     MPI_COMM_WORLD, &
                     reqr(nrecv), &
                     ierr)
               else
                  maxvals_parent(lb) = maxvals(parent(1,lb))
               end if
            end if
         end if
      end do

      ! parents send maxvals to children

      nsend = 0
      do i = 1, blkCount
         lb = gr_blkList(i)
         if (nodetype(lb) == PARENT_BLK .AND. lrefine(lb) == lrefine_max-1) then
            do j = 1,nchild
               if(child(1,j,lb).gt.-1) then
                  if (child(2,j,lb).ne.MyPE) then
                     nsend = nsend + 1
                     call MPI_ISend(maxvals(lb), &
                        1, &
                        FLASH_REAL, &
                        child(2,j,lb), &  ! PE TO SEND TO
                        child(1,j,lb), &  ! THIS IS THE TAG
                        MPI_COMM_WORLD, &
                        reqs(nsend), &
                        ierr)
                  end if
               end if
            end do
         end if
      end do

      if (nsend.gt.0) then
         call MPI_Waitall (nsend, reqs, stats, ierr)
      end if
      if (nrecv.gt.0) then
         call MPI_Waitall (nrecv, reqr, statr, ierr)
      end if

!!      maxvals_parent(:) = 0.0  ! <-- uncomment line for previous behavior
      do i = 1, blkCount
         lb = gr_blkList(i)
         if (nodetype(lb) == LEAF) then
            if (lrefine(lb) .gt. max_refine) then
                refine(lb) = .false.
                derefine(lb) = .true.
            elseif (maxvals(lb) < gr_refine_val_cutoff(1)) then
                refine(lb)   = .false.
!                if (maxvals_parent(lb) < gr_refine_val_cutoff .AND. .NOT. stay(lb)) derefine(lb)   = .true.
                if (maxvals_parent(lb) < gr_refine_val_cutoff(1)) derefine(lb)   = .true.
            endif
         end if
      enddo

  endif

  return
end subroutine Grid_markRefineDerefine

