!!****if* source/Grid/GridSolvers/Multipole_experimental/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!!
!! SYNOPSIS
!!
!!   Grid_solvePoisson(integer(IN) :: iSoln,
!!                     integer(IN) :: iSrc, 
!!                  integer(6)(IN) :: bcTypes,
!!                   real(2,6)(IN) :: bcValues,
!!                     real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the multipole
!!   summation method for isolated boundary problems.  Periodic
!!   problems are not supported; the multipole method should be
!!   used only for approximately spherical matter distributions.
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - index to variable containing potential
!!  iSrc     - index to variable containing density
!!  bcTypes  - boundary types along various faces,
!!             only used in verifying that they are isolated
!!  bcValues - the values to boundary conditions, currently not used
!!  poisfact - Poisson factor to be used in calculation
!!
!!***

subroutine Grid_solvePoisson (iSoln,             &
                              iSrc,              &
                              bcTypes,           &
                              bcValues,          &
                                         poisfact)

  use Grid_interface,    ONLY: GRID_PDE_BND_ISOLATED, &
                               Grid_getListOfBlocks
  use Timers_interface,  ONLY: Timers_start,          &
                               Timers_stop
  use gr_mpoleData,      ONLY: blockCount,            &
                               blockList,             &
                               dumpMoments,           &
                               printRadialInfo,       &
                               !Tot_Moment_R,          &
                               Moment_R,              &
                               max_Q,max_LM,          &
                               old_max_Q
                               !tot_max_Q,             &
                               !old_tot_max_Q
  use Gravity_data,      ONLY: grv_mode
  use gr_mpoleInterface, ONLY: & !gr_mpoleCopyTotMoments,&
                               gr_mpoleDeallocateRadialArrays,&
                               gr_mpoleAllocateRadialArrays,&
                               gr_mpoleCenterOfMass,&
                               gr_mpoleRadialSampling,&
                               !gr_mpoleDeallocateTotMoments,&
                               gr_mpolePrintRadialInfo,&
                               gr_mpolePotentials,&
                               gr_mpoleDumpMoments
  use Grid_data,         ONLY: gr_meshMe

  implicit none

#include "constants.h"

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(6)
  real,    intent(in)    :: bcValues(2,6)
  real,    intent(inout) :: poisfact
!  
!
!    ...Check the grid boundaries. Abort, if not isolated
!
!
  if (bcTypes(1) /= GRID_PDE_BND_ISOLATED) then
      call Driver_abortFlash ("FATAL: Multipole Poisson solver requires isolated boundaries")
  end if
!  
!
!    ...Start timer.
!
!
  call Timers_start ("Multipole Solver")
!  
!
!    ...Preliminary chores for setting up the run.
!
!
  call Grid_getListOfBlocks         (LEAF, blockList, blockCount)

  !grv_mode = 1
  if (allocated(Moment_R)) call gr_mpoleDeallocateRadialArrays ()
  !call gr_mpoleCenterOfMass           (iSrc)
  !call gr_mpoleRadialSampling         ()
  !call gr_mpoleAllocateRadialArrays   ()
  !call gr_mpoleSetDampingFactors      ()
  !call gr_mpoleMoments                (iSrc)
  !if (allocated(Tot_Moment_R)) call gr_mpoleDeallocateTotMoments()
  !call gr_mpoleCopyTotMoments         ()

  grv_mode = 2
  !call gr_mpoleDeallocateRadialArrays ()
  call gr_mpoleCenterOfMass           (iSrc)
  call gr_mpoleRadialSampling         ()
  call gr_mpoleAllocateRadialArrays   ()
  call gr_mpoleSetDampingFactors      ()

  !if (gr_meshMe .eq. MASTER_PE) then
  !    print *, 'cur. tot, peak, old tot, peak Q', max_Q, tot_max_Q, old_max_Q, old_tot_max_Q
  !endif

!
!
!     ...Print radial info if requested by user.
!
!
  if (printRadialInfo) then
      call gr_mpolePrintRadialInfo ()
  end if
!
!
!    ...Do the work.
!
!
  call gr_mpoleMoments                (iSrc)
  call gr_mpolePotentials     (iSrc, iSoln, poisfact)

!
!
!    ...Dump the moments if requested by the user.
!
!
  if (dumpMoments) then
      call gr_mpoleDumpMoments ()
  end if
!
!
!    ...Final chores.
!
!
! JFG: DISABLING THIS LINE SO MOMENTS ARE STILL AVAILABLE TO ORBIT_UPDATE
!  call gr_mpoleDeallocateRadialArrays ()
!  
!
!    ...End timer.
!
!
  call Timers_stop ("Multipole Solver")
!
!
!    ...Check the potentials obtained. This routine is only internal
!       and should only be activated by the author of the code for
!       debugging and accuracy check purposes.
!
!
!  call gr_mpolePotential_exact   (iSrc,iSoln, poisfact)
!  
!
!    ...Ready!
!
!
  return
end subroutine Grid_solvePoisson
