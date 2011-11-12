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
  integer, save :: grv_myPE, grv_numProcs, grv_dynRefineMax


  logical, save :: useGravity, updateGravity
  logical, save :: grav_temporal_extrp !extrapolate or otherwise rescale

  real,    save :: grav_poisfact
  integer, save :: grv_mode
  double precision, save :: grv_ptmass, grv_factor,&
      peri_dist, peri_time, grv_thresh, orb_ecc, peri_beta, period_fac, &
      start_beta, grv_comCutoff, grv_ener, grv_tot_ener, grv_bound, grv_denscut
  double precision, save :: orb_t, orb_dt, grv_cfl, grv_totmass
  double precision, dimension(6), save :: grv_ptvec, grv_obvec, grv_optvec, grv_oobvec, grv_boundvec, grv_exactvec

end module Gravity_data
