!!****if* source/Simulation/SimulationComposition/Burn/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed
!!  for setups that use nuclear networks.  The setups that want to use multispecies
!!  capabilities of the code for something other than nuclear burning and ionization
!!  should include their own custom implementation of this routine
!!
!!  This implementation of the routine relies on a textfile
!!  SpeciesList.txt to provide the species related informatio. The
!!  textfile contains the elements sorted by their atomic number in
!!  increasing order, and the isotopes of each element in turn sorted by
!!  of their atomic number, again in increasing order. The subroutine
!!  reads in the records in file, if the record corresponds to an
!!  isotope that is included in the setup, it sets the properties of
!!  the isotope in the multispecies database, and if the istope is not
!!  included, it goes on to read the next one. This process is
!!  repeated until all the species included in the setup have been
!!  found.
!!
!!  The format of SpeciesList.txt is as follow
!!  Column#         Variable                  Description
!!  ---------------------------------------------------------
!!  1               isotopeName               Sorted in increasing atomic number
!!  2               Z                         zbar, Atomic number; number of protons in nucleus
!!  3               A                         abar, total number of protons and neutrons in nucleus
!!  4               N                         Number of neutrons
!!  5               Eb                        binding energy
!!  6               ??
!!  7               ??
!!
!!
!!
!!  ARGUMENTS : There are no arguments in this subroutine
!!
!!  NOTE
!! 
!!***

subroutine Simulation_initSpecies()
  use Multispecies_interface, ONLY : Multispecies_setProperty
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_data, ONLY : sim_fluidGammac,sim_fluidGammae
  implicit none

#include "constants.h"
#include "Flash.h"
#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
  
!  integer, parameter :: SPEC_UNIT=2,SPEC_NUM=238
!  character(len=4)::isotopeName
!  real :: abar,zbar,bindEnergy,neutrons, electrons
!  integer :: i, isotope,count
  
  call RuntimeParameters_get('sim_fluidGammac',sim_fluidGammac)
  call RuntimeParameters_get('sim_fluidGammae',sim_fluidGammae)

!  open(unit=SPEC_UNIT,file="SpeciesList.txt")
!  count=0
!  i=0

!  do while((count<NSPECIES).and.(i<=SPEC_NUM))
!     i=i+1
!     read(SPEC_UNIT,*)isotopeName,zbar,abar,neutrons,bindEnergy
!     electrons = zbar

!     call Simulation_mapStrToInt(isotopeName,isotope,MAPBLOCK_UNK)
!     if(isotope /= NONEXISTENT) then
!        count=count+1
!        call Multispecies_setProperty(isotope, A, abar)
!        call Multispecies_setProperty(isotope, Z, zbar)
!        call Multispecies_setProperty(isotope, N, neutrons)
!        call Multispecies_setProperty(isotope, E, electrons)
!        call Multispecies_setProperty(isotope, EB, bindEnergy)
!        call Multispecies_setProperty(isotope, GAMMA, sim_fluidGamma)
!     end if
!  end do
!  close(SPEC_UNIT)
      call Multispecies_setProperty(H1_SPEC, A, 1.)
      call Multispecies_setProperty(HE4_SPEC, A, 4.)
      call Multispecies_setProperty(H1_SPEC, GAMMA, sim_fluidGammae)
      call Multispecies_setProperty(HE4_SPEC, GAMMA, sim_fluidGammae)
      call Multispecies_setProperty(CORE_SPEC, A, 5.16d0)
      call Multispecies_setProperty(CORE_SPEC, GAMMA, sim_fluidGammac)
#endif
end subroutine Simulation_initSpecies
