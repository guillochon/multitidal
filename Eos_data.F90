!!****if* source/physics/Eos/EosMain/Multigamma/Eos_data
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
!!  This is the data module for the Gamma law Eos implementation with 
!!  multiple fluids/species. 
!!  It stores all the runtime parameters, and all the unit scope
!!  data. Some of the unit scope data is fecthed by the wrapper layer
!!  from elsewhere in the code and some is local unit data common to
!!  multiple functions in the unit 
!! 
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Multigamma Eos.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!   smalle[Real]  --- the smallest value for energy 
!!
!!***

module Eos_data

#include "Flash.h"
#include "Eos.h"
#include "Eos_map.h"

  real :: eos_gasConstant
  real :: eos_smalle
  real :: eos_gamma
  real :: eos_eintSwitch
  real :: eos_largee

#ifdef FIXEDBLOCKSIZE
  real,save,dimension(MAXCELLS*NSPECIES) :: eos_massFr
  real,save,dimension(EOS_NUM*MAXCELLS) :: eos_inOut
#else
  real,save, allocatable, dimension(:) :: eos_inOut,eos_massFr
#endif
  real, save, dimension(NSPECIES) ::  eos_gc, eos_gammam1j, eos_ggprodj, eos_ggprodinvj, eos_gam1invj
  integer, save, dimension(1:EOSMAP_NUM_ROLES, 1:2, 1:5) :: eos_mapLookup
end module Eos_data
