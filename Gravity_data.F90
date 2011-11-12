!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!  
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  This modules stores the data for the Gravity unit
!!
!!***

module Gravity_data

#include "constants.h"

  character(len=MAX_STRING_LENGTH), save :: grav_boundary_type !string boundary condition
  integer, save :: grav_boundary  !integer boundary condition

  integer, save :: grav_geometry  !mesh geometry
  integer, save :: grv_myPE, grv_numProcs


  logical, save :: useGravity, updateGravity
  logical, save :: grav_temporal_extrp !extrapolate or otherwise rescale

  real,    save :: grav_poisfact
  real,    save :: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_ptmass, grv_factor,&
      peri_dist, peri_time, grv_thresh

end module Gravity_data
