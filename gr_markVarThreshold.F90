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

  use tree, ONLY : refine, derefine, lrefine, lnblocks, nodetype
  use physicaldata, ONLY : unk
  implicit none
#include "constants.h"
! Arguments

  integer, intent(IN) :: Var
  real,    intent(IN) :: var_th
  integer, intent(IN) :: icmp, lref

! Local data

  integer :: b, llref
  logical :: Grid_mark, Grid_mark2
  logical :: unmark = .false.

!-------------------------------------------------------------------------------

! Loop over all leaf-node blocks.

  if (lref < 0) unmark = .true.
  llref = abs(lref)

  do b = 1, lnblocks
    if (nodetype(b) == LEAF) then

! Compare the variable against the threshold.

      if (unmark) then
        if (icmp < 0) then
          Grid_mark = (maxval(unk(var,:,:,:,b)) < var_th)
          Grid_mark2 = (maxval(unk(var,:,:,:,b)) < 8.d0*var_th)
        else
          Grid_mark = (minval(unk(var,:,:,:,b)) > var_th)
          Grid_mark2 = (minval(unk(var,:,:,:,b)) > var_th/8.d0)
        endif
      else
        if (icmp < 0) then
          Grid_mark = (minval(unk(var,:,:,:,b)) < var_th)
          Grid_mark2 = (minval(unk(var,:,:,:,b)) < 8.d0*var_th)
        else
          Grid_mark = (maxval(unk(var,:,:,:,b)) > var_th)
          Grid_mark2 = (maxval(unk(var,:,:,:,b)) > var_th/8.d0)
        endif
      endif

! If the test passed, Grid_mark the block for refinement.

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
        endif

      endif

      if (Grid_mark2) then

        if (unmark) then
          if (lrefine(b) >= llref ) then
            refine(b)   = .false.
          endif
        else
          if (lrefine(b) <= llref ) then
            derefine(b) = .false.
          endif
        endif

      endif

    endif
  enddo

!-------------------------------------------------------------------------------

  return
end subroutine gr_markVarThreshold
