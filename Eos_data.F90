!!****if* source/physics/Eos/EosMain/Helmholtz/Eos_data
!!
!! NAME
!!
!!  Eos_data
!!
!! 
!! SYNOPSIS
!!
!!  use Eos_data
!!
!! DESCRIPTION
!!
!!  General parameters (non-array) for EOS Helmholtz
!!
!! ARGUMENTS
!!
!!
!!*** 

module Eos_data

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"  
#include "Eos_map.h"

  ! maximum number of iterations for the Newton loop to find T from e
  integer, save :: eos_maxNewton
 
  ! general grid parameters
  integer, save :: eos_myPE, eos_numProcs 
  
  ! how accurately to do the Newton iteration to get T from e
  real, save :: eos_tol
  
  real, save :: eos_smallt, eos_larget
  
  real,save :: eos_gasConstant 
  real,save :: eos_smalle
  real,save :: eos_singleSpeciesA, eos_singleSpeciesZ  ! only used in Ye formulation
  real,save :: eos_gamma
  real,save :: eos_gammam1
  real, save :: eos_eintSwitch

  integer,save :: eos_hfetInit

  ! force the iterative solver to leave inputs alone (always true in MODE_DENS_TEMP)
  logical, save :: eos_forceConstantInput

  ! Coulomb multiplier 
  real, save :: eos_coulombMult
  ! abort if pressures become negative
  logical, save :: eos_coulombAbort
 
  integer,parameter :: EOSIMAX=271,EOSJMAX=101


  real, save :: eos_tlo, eos_tstpi
  real, save :: eos_dlo, eos_dstpi
  real,dimension(EOSJMAX),save :: eos_dt,eos_dtSqr,eos_dtInv,eos_dtSqrInv,eos_t
  real,dimension(EOSIMAX),save :: eos_dd,eos_ddSqr,eos_ddInv,eos_ddSqrInv,eos_d

!..for the helmholtz free energy tables
!..for the pressure derivative with density tables
!..for the chemical potential tables
!..for the number density tables
  real,save,dimension(EOSIMAX,EOSJMAX) :: eos_f,eos_fd, eos_ft,eos_fdd,&
                                          eos_ftt,eos_fdt,eos_fddt,&
                                          eos_fdtt, eos_fddtt, & 
                                          eos_dpdf,eos_dpdfd,eos_dpdft,&
                                          eos_dpdfdd,eos_dpdftt,eos_dpdfdt,&
                                          eos_ef,eos_efd,eos_eft,eos_efdd,&
                                          eos_eftt,eos_efdt, & 
                                          eos_xf,eos_xfd,eos_xft,eos_xfdd,&
                                          eos_xftt,eos_xfdt

  integer, save, dimension(1:EOSMAP_NUM_ROLES, 1:2, 1:5) :: eos_mapLookup
end module Eos_data
