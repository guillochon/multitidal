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
  integer, save :: grv_dynRefineMax
  integer, save :: grv_meshMe, grv_meshNumProcs, grv_meshComm, grv_bndMassIters
  integer, save :: grv_commSize=1

  logical, save :: useGravity, updateGravity
  logical, save :: grav_temporal_extrp !extrapolate or otherwise rescale
  logical, save :: grav_unjunkPden

  real,    save :: grav_poisfact
  integer, save :: grv_mode
  double precision, save :: grv_ptmass, grv_factor, grv_periDist, grv_thresh, &
      grv_comCutoff, grv_comPeakCut, grv_ener, grv_tot_ener, grv_bound, grv_denscut
  double precision, save :: orb_t, orb_dt, grv_totmass, grv_totmass0, grv_periTime
  double precision, dimension(6), save :: grv_ptvec, grv_obvec, grv_optvec, grv_oobvec, &
      grv_boundvec, grv_exactvec, grv_mpolevec, grv_oexactvec, grv_ompolevec
  double precision, dimension(6), save :: grv_hptvec, grv_hobvec, grv_peakvec
  double precision, dimension(3), save :: grv_obaccel, grv_ptaccel, grv_hptaccel, grv_optaccel, grv_oobaccel, &
      grv_momacc, grv_angmomacc, grv_mpoleaccel, grv_ompoleaccel, grv_o2mpoleaccel, grv_o2obaccel
  double precision, save :: grv_orbTol, grv_orbMinForce, grv_ototmass, grv_optmass, &
      grv_eneracc, grv_massacc, grv_dtOld, grv_dt2Old, grv_cfl, grv_orbStepFrac
  double precision, dimension(3,3), save :: grv_rotMat, grv_invRotMat
  logical, save :: grv_orb3D

end module Gravity_data
