!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro_computeDt
!!
!! NAME
!!
!!  Hydro_computeDt
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
!!  Computes the timestep limiter for point mass gravity.
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


subroutine Gravity_computeDt ( dtCheck )
     
  
#include "Flash.h"
#include "constants.h"

  use Hydro_data, ONLY : hy_cfl
  use Grid_interface, ONLY: Grid_getMinCellSize
  use Driver_data, ONLY: dr_simTime
  use Simulation_data, ONLY: sim_tRelax, sim_starRadius
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Gravity_data, ONLY: grv_ptmass
  use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, Mtot

  implicit none

  real,INTENT(INOUT)    :: dtCheck
  
  real    :: dt_temp, min_cell_size, vel, tinitial, r2, newx, newy, newz
             
  
!==============================================================================

  call RuntimeParameters_get('tinitial',tinitial)
  if (dr_simTime .lt. tinitial + sim_tRelax) return

  call parabolic_orbit(dr_simTime, newx, newy, newz)
  call Grid_getMinCellSize(min_cell_size)
  r2 = (newx - Xcm)**2. + (newy + Ycm)**2. + (newz + Zcm)**2.

  dt_temp = hy_cfl*min_cell_size/vel*max(1.0,r2/(sim_starRadius*(2.*grv_ptmass/Mtot)**(1./3.))**2.)
     
  if (dt_temp < dtCheck) then
     dtCheck = dt_temp
  endif
  
  return
end subroutine Gravity_computeDt


