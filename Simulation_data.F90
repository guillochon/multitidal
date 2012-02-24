!!****if* source/Simulation/SimulationMain/MultiTidalPoly/Simulation_data
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
!!  Stores the local data for Simulation setup: MultiTidalPoly
!!  
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"

  !! *** Runtime Parameters *** !!

  double precision, save :: sim_pAmbient, sim_rhoAmbient, &
                            sim_fluidGammac, sim_fluidGammae, &
                            sim_smallX, sim_smallRho, sim_smallP, &
                            sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  integer, save          :: sim_nSubZones
  integer, save          :: sim_tableRows, sim_tableCols
  integer, save          :: sim_maxBlocks, sim_ptMassRefine, sim_totForceSub
  double precision, save :: sim_tInitial, sim_tRelax, sim_relaxRate, sim_starRadius, &
                            sim_softenRadius, sim_accRadius, sim_objMass, sim_objCoreMass, sim_objPolyN, sim_objCentDen, &
                            sim_objPolyN2, sim_objRadius, sim_objCore, sim_objEnve, &
                            sim_accCoeff, sim_fluffDampCoeff, sim_fluffDampCutoff, sim_ptMass, &
                            sim_periBeta, sim_startBeta, sim_periodFac, &
                            sim_orbEcc, sim_objAbarCore, sim_startDistance, sim_ptMassRefRad, &
                            sim_totForceInv


  !! *** Variables pertaining to this Simulation *** !!

  integer, parameter        :: np = 1000
  double precision, save    :: sim_inSubInv
  double precision, save    :: sim_inszd
  double precision, save :: obj_muc, obj_mue
  double precision, dimension(SPECIES_BEGIN:SPECIES_END) :: obj_xn
  integer, dimension(SPECIES_BEGIN:SPECIES_END) :: speciesMask
  double precision, dimension(np), save :: obj_radius, obj_rhop, obj_prss
  double precision, save :: gammac, gammae
  integer, save :: obj_ipos, obj_ipoi
  logical, save :: sim_useInitialPeakDensity

  double precision, parameter :: sim_msun = 1.9889225d33

end module Simulation_data
