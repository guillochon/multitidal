!!****f* source/physics/sourceTerms/Cool/Cool_computeDt
!!
!!
!! NAME
!!  
!!  Cool_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Cool_computeDt ( integer(IN) : blockID, 
!!                   
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_check, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for heating source term solver.
!! 
!!
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***

subroutine Cool_computeDt (blockID, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_check, dt_minloc )
#include "Flash.h"
#include "constants.h"

  implicit none

  integer, intent(IN) :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dt_check
  integer,INTENT(INOUT)    :: dt_minloc(5)
  real, pointer :: solnData(:,:,:,:) 

  return
end subroutine Cool_computeDt


