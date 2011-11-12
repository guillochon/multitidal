!!****if* source/physics/Gravity/GravityMain/Poisson/Multipole/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!
!! 
!! SYNOPSIS
!!
!!  Gravity_init(integer(IN) :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!!  ARGUMENTS
!!
!!  myPE - local processor number
!!
!!***

subroutine Gravity_init(myPE)

  use Gravity_data
  use Grid_interface, ONLY: Grid_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none
  real,save :: newton

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: myPE
  character(len=MAX_STRING_LENGTH) :: strGeometry

  ! Everybody should know these
  grv_myPE = myPE
  call Grid_getNumProcs(grv_numProcs)


  call PhysicalConstants_get("Newton", newton)
  
  call RuntimeParameters_get("geometry", strGeometry)
  call RuntimeParameters_mapStrToInt(strGeometry, grav_geometry)
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  call RuntimeParameters_get("ptxpos", grv_ptxpos)
  call RuntimeParameters_get("ptypos", grv_ptypos)
  call RuntimeParameters_get("ptzpos", grv_ptzpos)
  call RuntimeParameters_get("ptmass", grv_ptmass)
  call RuntimeParameters_get("sim_periDist", peri_dist)
  call RuntimeParameters_get("sim_periTime", peri_time)
  call RuntimeParameters_get("useGravity", useGravity)
  call RuntimeParameters_get("updateGravity", updateGravity)

  grav_poisfact = 4. * PI * newton
  grv_factor = -newton * grv_ptmass
  grv_thresh = 1.0D-10

  call gr_mpoleCenterOfMass(DENS_VAR)

  return
end subroutine Gravity_init
