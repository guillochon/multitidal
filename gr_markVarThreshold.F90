!!****if* source/Grid/GridMain/paramesh/gr_markVarThreshold
!!
!! NAME
!!  gr_markVarThreshold
!!
!! SYNOPSIS
!!  gr_markVarThreshold(  integer(in) :: Var
!!                         real(in)   :: var_th, 
!!                         integer(in) :: icmp, 
!!                         integer(in) :: lref )
!!
!! PURPOSE
!!  Refine all blocks for which a given variable (Var) exceeds or falls
!!  below some threshold (var_th).  The direction of the threshold is
!!  controlled by the parameter icmp.  Either blocks are brought
!!  up to a specific level of refinement or each block is refined once.
!!
!! ARGUMENTS
!!  Var -    the variable of interest
!!
!!  var_th  -     the limit on the variable
!! 
!!  icmp  -  icmp < 0  (de)refine if the variable is less than var_th
!!         icmp >= 0 (de)refine if the variable is greater then var_th
!! 
!!   lref -       If > 0, bring all qualifying blocks to this level of refinement.
!!
!!               If <= 0, refine qualifying blocks once.
!!
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!
!!***

!!REORDER(5):unk

subroutine gr_markVarThreshold (Var, var_th, icmp,lref)

  use paramesh_dimensions
  use tree
  use physicaldata, ONLY : unk
  use Grid_data, ONLY : gr_meshMe, gr_meshComm

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"
! Arguments

  integer, intent(IN) :: Var
  real,    intent(IN) :: var_th
  integer, intent(IN) :: icmp, lref

! Local data

  integer :: b, llref, nsend, nrecv, ierr, j
  logical :: Grid_mark
  logical :: unmark
  real val(maxblocks),val_par(maxblocks)
  integer reqr(maxblocks),reqs(maxblocks*nchild)
  integer :: statr(MPI_STATUS_SIZE,maxblocks)
  integer :: stats(MPI_STATUS_SIZE,maxblocks*nchild)

!-------------------------------------------------------------------------------

  unmark = .false.
  if (lref < 0) unmark = .true.
  llref = abs(lref)

  do b = 1, lnblocks
    val(b) = 0.
    if (nodetype(b).eq.LEAF.or.nodetype(b).eq.PARENT_BLK) then
      if (unmark) then
        if (icmp < 0) then
          val(b) = maxval(unk(var,:,:,:,b))
        else
          val(b) = minval(unk(var,:,:,:,b))
        endif
      else
        if (icmp < 0) then
          val(b) = minval(unk(var,:,:,:,b))
        else
          val(b) = maxval(unk(var,:,:,:,b))
        endif
      endif
    endif
  enddo

  val_par(1:lnblocks) = 0.
  nrecv = 0
  do b = 1,lnblocks
     if(parent(1,b).gt.-1) then
        if (parent(2,b).ne.gr_meshMe) then
           nrecv = nrecv + 1
           call MPI_IRecv(val_par(b),1, &
                MPI_DOUBLE_PRECISION, &
                parent(2,b), &
                b, &
                gr_meshComm, &
                reqr(nrecv), &
                ierr)
        else
           val_par(b) = val(parent(1,b))
        end if
     end if
  end do
 
  ! parents send val to children

  nsend = 0
  do b = 1,lnblocks
     do j = 1,nchild
        if(child(1,j,b).gt.-1) then
           if (child(2,j,b).ne.gr_meshMe) then
              nsend = nsend + 1
              call MPI_ISend(val(b), &
                   1, &
                   MPI_DOUBLE_PRECISION, &
                   child(2,j,b), &  ! PE TO SEND TO
                   child(1,j,b), &  ! THIS IS THE TAG
                   gr_meshComm, &
                   reqs(nsend), &
                   ierr)
           end if
        end if
     end do
  end do

  if (nsend.gt.0) then
     call MPI_Waitall (nsend, reqs, stats, ierr)
  end if
  if (nrecv.gt.0) then
     call MPI_Waitall (nrecv, reqr, statr, ierr)
  end if

  do b = 1,lnblocks
    if (nodetype(b).eq.LEAF .or. nodetype(b).eq.PARENT_BLK) then
      if (unmark) then
        if (icmp < 0) then
          Grid_mark = (val(b) <= var_th .and. val_par(b) <= var_th)
        else
          Grid_mark = (val(b) >= var_th .and. val_par(b) >= var_th)
        endif
      else
        if (icmp < 0) then
          Grid_mark = (val(b) <= var_th)
        else
          Grid_mark = (val(b) >= var_th)
        endif
      endif

      if (Grid_mark) then
        if (unmark) then
          if (lrefine(b) > llref ) then
            refine(b)   = .false.
            derefine(b) = .true.
          else if (lrefine(b) == llref) then
            refine(b) = .false.
          else if (llref <= 0) then
            derefine(b) = .true.
          endif
        else
          if (lrefine(b) < llref ) then
            refine(b)   = .true.
            derefine(b) = .false.
          else if (lrefine(b) == llref) then
            derefine(b) = .false.
          else if (llref <= 0) then
            refine(b) = .true.
          endif
          if (nodetype(b).eq.LEAF) then
            stay(b) = .true.
          endif
        endif
      endif
    endif
  enddo

!-------------------------------------------------------------------------------

  return
end subroutine gr_markVarThreshold
