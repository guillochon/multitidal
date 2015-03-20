!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/Constant/MagneticResistivity
!!
!! NAME
!!  MagneticResistivity
!!
!! SYNOPSIS
!!  call MagneticResistivity(real,    intent(IN)  :: temp,
!!                      real,    intent(IN)  :: dens,
!!                      real,    intent(IN)  :: xn(NSPECIES),
!!                      real,    intent(OUT) :: magResist)
!!
!! DESCRIPTION
!!  Implementation for constant magnetic resistivity.
!!
!! ARGUMENTS
!!  temp      - Plasma temperature
!!  dens      - Plasma density
!!  xn        - Species
!!  magResist - Magnetic resistivity
!!
!! NOTE
!!  In previous versions, there used to be two runtime parameters:
!!  magnetic resistivity and magnetic viscosity.
!!  They respectively refer "eta" and "eta*c*c/(4*pi)" in CGS
!!  (or "eta/mu_0" in SI, where mu_0=4*pi*10^(-7) henry/meter).
!!  What it was done in the old way was to initialize magnetic
!!  viscosity (e.g., eta*c*c/(4*pi)) using the magnetic resistivity, eta.
!!  In this version, such distinctions between the magnetic resistivity
!!  and magnetic viscosity has been removed and we only use
!!  magnetic resistivity with proper scalings depending on unit system.
!!
!!
!!***

subroutine MagneticResistivity(temp,dens,xn,magResist)

  use MagneticResistivity_data, ONLY : mResistivity, mUnit, c
  use Simulation_data, ONLY : sim_magResistCutoff
  use Multitidal_interface, ONLY : Multitidal_findExtrema

#include "constants.h"
#include "Flash.h"

  implicit none

  !! Argument list -------------------------------
  real, intent(IN) :: temp, dens
  real, intent(IN), dimension(NSPECIES) :: xn
  real, intent(OUT):: magResist
  !! ----------------------------------------------

  real :: damp_dens

  call Multitidal_findExtrema(DENS_VAR, 1, damp_dens)

  !! Scale resistivity depending on units
  if (dens .lt. sim_magResistCutoff*damp_dens) then
      if (mUnit == "SI" .or. mUnit == "si" ) then
         magResist = mResistivity*1.e7/(4.*PI)
      elseif (mUnit == "CGS" .or. mUnit == "cgs" ) then
         magResist = mResistivity*c*c/(4.*PI)
      else !no unit
         magResist = mResistivity
      end if
  endif

end subroutine MagneticResistivity
