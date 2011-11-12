!!****if* source/physics/Eos/EosMain/Multigamma/Eos_init
!!
!! NAME
!!
!!  Eos_init
!!
!! 
!! SYNOPSIS
!!
!!  call Eos_init(integer(in) :: myPE)
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars needed
!!  by the EOS unit from the runtime parameters and physical
!!  constants facilities. This version is for use when multiple species
!!  are present. The Gamma's for different species are obtained from
!!  the Multispecies unit, initialized in Simulation_initSpecies.F90
!!
!! ARGUMENTS
!!
!!  myPE - current processor (ignored by this implementation)
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Gamma law Eos for multiple
!!   species with different abundances. 
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might over write these values with the 
!!   flash.par values for your specific run.  
!!
!!  NOTES
!!
!!  Gamma law Eos defines two mesh-based parameters GAMC_VAR and GAME_VAR in Flash.h
!!
!!
!!***

subroutine Eos_init(myPE)

  use Eos_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Driver_interface, ONLY: Driver_abortFlash

 
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"

  integer, intent(IN) :: myPE

  integer :: numCells, istat
  integer :: ispecies

  call PhysicalConstants_get("ideal gas constant", eos_gasConstant)

  call RuntimeParameters_get("smalle",eos_smalle)
  call RuntimeParameters_get("eintSwitch",eos_eintSwitch)
  call RuntimeParameters_get("eos_largee",eos_largee)
#ifndef EINT_VAR
  if (eos_eintSwitch > 0.0) then
     call Driver_abortFlash("[Eos_init] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif

  do ispecies = 1, NSPECIES
     call Multispecies_getProperty(SPECIES_BEGIN + ispecies - 1, GAMMA,  eos_gc(ispecies))
  end do

  ! Note that these are all ARRAYS of size NSPECIES
  eos_gammam1j   = 1. / (eos_gc - 1.)
  eos_ggprodj    = eos_gammam1j * eos_gasConstant
  eos_ggprodinvj = 1. / eos_ggprodj
  eos_gam1invj   = 1. / eos_gammam1j

  call eos_fillMapLookup()
  return
end subroutine Eos_init
