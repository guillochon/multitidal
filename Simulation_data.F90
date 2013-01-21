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

  double precision, save :: sim_pAmbient, sim_tAmbient, sim_rhoAmbient, &
                            sim_fluidGamma, &
                            sim_smallX, sim_smallRho, sim_smallP, &
                            sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  integer, save          :: sim_nSubZones
  integer, save          :: sim_maxBlocks, sim_ptMassRefine, sim_totForceSub
  double precision, save :: sim_tInitial, sim_tRelax, sim_relaxRate, sim_objRadius, &
                            sim_softenRadius, sim_accRadius, sim_objCentDens, sim_objMass, &
                            sim_accCoeff, sim_fluffDampCoeff, sim_fluffDampCutoff, sim_ptMass, &
                            sim_periBeta, sim_startBeta, sim_periodFac, &
                            sim_orbEcc, sim_startDistance, sim_ptMassRefRad, &
                            sim_totForceInv, sim_rotFac, sim_rotAngle, sim_tSpinup
  double precision, save :: sim_xCenter, sim_yCenter, sim_zCenter

#ifdef LOADPROFILE
  integer, save                                       :: sim_tableRows, sim_tableCols
  double precision, dimension(:,:), allocatable, save :: sim_table
  character(len=MAX_STRING_LENGTH), save              :: sim_profFile
#else
  integer, parameter                                     :: np = 1000
  double precision, save                                 :: sim_objPolyN, sim_objCentDen, obj_mu
  double precision, dimension(SPECIES_BEGIN:SPECIES_END) :: obj_xn
  double precision, dimension(np), save                  :: obj_radius, obj_rhop, obj_prss
  integer, save                                          :: obj_ipos
#endif

  !! *** Variables pertaining to this Simulation *** !!

  double precision, save    :: sim_inSubInv
  double precision, save    :: sim_inszd
  logical, save :: sim_useInitialPeakDensity

  double precision, parameter :: sim_msun = 1.9889225d33

end module Simulation_data
