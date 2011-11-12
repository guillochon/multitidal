!!****if* source/physics/Hydro/HydroMain/split/PPM/Gravity_computeDt
!!
!! NAME
!!
!!  Gravity_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Gravity_computeDt(integer(IN):: blockID, 
!!                  integer(IN) :: myPE, 
!!                  real(IN) :: x(:), 
!!                  real(IN) :: dx(:), 
!!                  real(IN) :: uxgrid(:),
!!                  real(IN) :: y(:), 
!!                  real(IN) :: dy(:), 
!!                  real(IN) :: uygrid(:), 
!!                  real(IN) :: z(:), 
!!                  real(IN) :: dz(:), 
!!                  real(IN) :: uzgrid(:), 
!!                  integer(IN) :: blkLimits(2,MDIM)
!!                  integer(IN) :: blkLimitsGC(2,MDIM)
!!                  real,pointer ::  solnData(:,:,:,:),   
!!                  real,(INOUT) ::   dtCheck, 
!!                  integer(INOUT) :: dtMinLoc(:) )
!!
!! DESCRIPTION
!!
!!  Computes the timestep limiter for the hydrodynamical solver.  For pure
!!  hydrodynamics, the Courant-Fredrichs-Lewy criterion is used.  The sound
!!  speed is computed and together with the velocities, is used to constrain
!!  the timestep such that no information can propagate more than one zone
!!  per timestep.
!!
!!
!! ARGUMENTS
!!
!!  blockID -       local block ID
!!  myPE -          local processor number
!!  x, y, z -       coordinates
!!  dx, dy, dz -    deltas in each {x, y z} directions
!!  uxgrid, uygrid, uzgrid - velocity of grid expansion in {x, y z} directions
!!  blkLimits -    the indices for the interior endpoints of the block
!!  blkLimitsGC - the indices for endpoints including the guardcells
!!  solnData -      the physical, solution data from grid
!!  dtCheck -      variable to hold timestep constraint
!!  dtMinLoc(5) -  array to hold location of cell responsible for minimum dt:
!!                 dtMinLoc(1) = i index
!!                 dtMinLoc(2) = j index
!!                 dtMinLoc(3) = k index
!!                 dtMinLoc(4) = blockID
!!                 dtMinLoc(5) = myPE
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Gravity_computeDt (blockID, myPE, &
                           blkLimits,blkLimitsGC,        &
                           solnData,   &
                           dtCheck, dtMinLoc )
     
    
#include "Flash.h"
#include "constants.h"

    use Driver_interface, ONLY : Driver_abortFlash
    use Driver_data, ONLY : dr_dtOld
    use Gravity_data, ONLY : grv_cfl

    implicit none


    integer, intent(IN) :: blockID, myPE
    integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
    real,INTENT(INOUT)    :: dtCheck
    integer,INTENT(INOUT)    :: dtMinLoc(5)
    real, pointer :: solnData(:,:,:,:) 
    
    integer :: i, j, k, temploc(5)
    real    :: dt_temp, dt_ltemp
    
!==============================================================================

    dt_temp    = 0.
    temploc(:) = 0
    
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)        
       do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)                         
          do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             dt_ltemp = abs(solnData(GPOT_VAR,i,j,k) - solnData(GPOL_VAR,i,j,k))/abs(solnData(GPOT_VAR,i,j,k))/dr_dtOld
             
             if (dt_ltemp > dt_temp) then
                dt_temp    = dt_ltemp
                temploc(1) = i
                temploc(2) = j
                temploc(3) = k
                temploc(4) = blockID
                temploc(5) = MyPE
             endif

          enddo           
       enddo
    enddo

    if (dt_temp .eq. 0.d0) then
        dt_temp = 1.d20
    else
        dt_temp = grv_cfl / dt_temp
    endif
    
    if (dt_temp < dtCheck) then
        dtCheck = dt_temp
        dtMinLoc = temploc
    endif
    
    return
end subroutine Gravity_computeDt


