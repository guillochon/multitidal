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
  integer, save          :: sim_nSubZones, sim_fixedParticle, sim_fixedPartTag
  integer, save          :: sim_maxBlocks, sim_ptMassRefine, sim_totForceSub, sim_cylinderType
  double precision, save :: sim_tInitial, sim_tRelax, sim_relaxRate, sim_objRadius, &
                            sim_softenRadius, sim_accRadius, sim_objCentDens, sim_objMass, &
                            sim_accCoeff, sim_fluffDampCoeff, sim_fluffDampCutoff, sim_ptMass, &
                            sim_periBeta, sim_startBeta, sim_periodFac, sim_periTime, &
                            sim_orbEcc, sim_startDistance, sim_ptMassRefRad, &
                            sim_totForceInv, sim_rotFac, sim_rotAngle, sim_tSpinup, &
                            sim_powerLawScale, sim_powerLawExponent, sim_powerLawExtent, &
                            sim_powerLawMass, sim_powerLawTemperature, sim_periDist, &
                            sim_fluffRefineCutoff, sim_xRayFraction, sim_starPtMass, &
                            sim_smallT, sim_cylinderScale, sim_cylinderDensity, &
                            sim_cylinderTemperature, sim_cylinderRadius, sim_condCoeff, &
                            sim_cylinderNCells, sim_cylinderMDot, sim_windVelocity, &
                            sim_windMdot, sim_windTemperature, sim_windNCells, &
                            sim_windLaunchRadius, sim_windKernel, sim_tDelay
  double precision, save :: sim_xCenter, sim_yCenter, sim_zCenter

  double precision, dimension(:,:), allocatable, save    :: sim_table
  integer, parameter                                     :: np = 1000
  double precision, save                                 :: sim_objPolyN, obj_mu, obj_gamc
  double precision, dimension(SPECIES_BEGIN:SPECIES_END) :: obj_xn
  double precision, dimension(np), save                  :: obj_radius, obj_rhop, obj_prss
  integer, save                                          :: obj_ipos, sim_tableRows, sim_tableCols
  character(len=MAX_STRING_LENGTH), save                 :: sim_profFile, sim_kind

  !! *** Variables pertaining to this Simulation *** !!

  double precision, save    :: sim_inSubInv
  double precision, save    :: sim_inszd
  double precision, dimension(6), save :: bhvec, stvec
  logical, save :: sim_useInitialPeakDensity, sim_useRadialProfile

  double precision, parameter :: sim_msun = 1.9889225d33

end module Simulation_data
