!!****f* source/physics/Gravity/Gravity_computeDt
!!
!! NAME
!!
!!  Gravity_computeDt
!!  
!! SYNOPSIS
!!
!!  Gravity_computeDt(integer(IN)        :: blockID,
!!                    
!!                    real (OUT)         :: dt_grav,
!!                    integer(:)(INOUT)  :: dt_minloc(5))
!!
!! DESCRIPTION
!!
!!  Compute the timestep limiter due to the gravitational solver.
!!
!! ARGUMENTS
!!
!!  MyPE:          local processor number
!!  dt_grav:       Will Return the limiting timestep. Should be
!!                 set to a large value (1.D99) on input.
!!  dt_minloc(5):  An array to receive information about which
!!                 processor, block, and zone was responsible
!!                 for setting the limiting timestep.  The order
!!                 is i, j, k, b, p, where (i,j,k) = zone
!!                 indices, b = local block ID, and p = PE #.
!!                 This routine should only modify these values
!!                 if it changes dt_grav.
!!  blockID:       The local ID of the block to compute the
!!                 limiter on.
!!
!!***

subroutine Gravity_computeDt (blockID, dt_grav, dt_minloc)

!==============================================================================
  use Grid_interface, ONLY : Grid_getMinCellSize
  use Gravity_data, ONLY : grv_obvec, grv_ptvec

  implicit none
  
  integer, intent(IN)    ::  blockID
  
  integer, intent(INOUT) ::  dt_minloc(5)
  real,intent(OUT)       ::  dt_grav

  real :: mcs
  
  call Grid_getMinCellSize(mcs)

  !dt_grav = min(huge(dt_grav), mcs/sqrt(sum((grv_obvec(4:6) - grv_ptvec(4:6))**2.d0)))
  dt_grav = huge(dt_grav)
  
  return

end subroutine Gravity_computeDt
