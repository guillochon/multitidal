!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"

  !! *** Runtime Parameters *** !!

  double precision, save :: sim_pAmbient, sim_rhoAmbient
  double precision, save :: sim_fluidGamma, sim_xCenter, sim_yCenter, sim_zCenter
  double precision, save :: sim_smallX, sim_smallRho, sim_smallP, sim_pi
  double precision, save :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  integer, save          :: sim_nSubZones, sim_sinkRef
  integer, save          :: sim_tableRows, sim_tableCols
  integer, save          :: sim_maxBlocks
  double precision, save :: sim_tInitial, sim_tRelax, sim_relaxRate, sim_starRadius
  double precision, save :: sim_sinkRadius, sim_objMass, sim_objPolyN, sim_objCentDen
  double precision, save :: sim_fluffDampCoeff, sim_fluffDampCutoff

  !! *** Variables pertaining to this Simulation *** !!

  integer, parameter        :: np = 1000
  double precision, save    :: sim_inSubZones, sim_inSubzm1
  double precision, save    :: sim_inszd
  double precision, save :: obj_mu
  double precision, dimension(SPECIES_BEGIN:SPECIES_END) :: obj_xn
  double precision, dimension(np), save :: obj_radius, obj_rhop, obj_prss
  integer, save :: obj_ipos

end module Simulation_data
