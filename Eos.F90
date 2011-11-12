!!****if* source/physics/Eos/EosMain/Multigamma/Eos
!!
!! NAME
!!
!!  Eos
!!
!! SYNOPSIS
!! 
!! 
!!  subroutine Eos(integer, intent(IN)           :: mode,
!!                 integer, intent(IN)           :: vecLen,
!!                 real, intent(INOUT)           :: eosData(:),
!!                 real, optional, intent(IN)    :: massFrac(:)
!!                 logical, optional, intent(IN) :: mask(:))
!!                 
!!
!! DESCRIPTION
!!
!!
!!  This routine applies the equation of state to thermodynamic quantities
!!  at one or more grid points.  The number of points is determined by the
!!  argument veclen.  Data should be packaged for this routine in the 1d
!!  array, eosData.  The data in eosData is organized as: 1:vecLen
!!  points contain the first variable, vecLen+1:2*vecLen points contain
!!  the second variable and so on. The number and 
!!  order of variables in the array is determined by the constants defined
!!  in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and internal energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  internal energy and temperature are generated as output. Note that
!!  internal energy is EINT_VAR, not ENER_VAR.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument) .true. means return the quantity, .false. means don't.
!!  
!!  This is a multigamma version, which means there are multiple species
!!  in the fluid, each in different abundances, and each with
!!  a different gamma.  This eos takes into account the contribution to
!!  the thermodynamic properties of the gas from each species
!!  appropriately.
!!  
!!  The argument, massFrac, holds the mass fractions in an order determined
!!  by the Multispecies unit.  
!!  
!!  
!! ARGUMENTS 
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on. 
!!
!!  massFrac : Contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!             An implementation that does not need derivative quantities should
!!             set the mask equal to .false.
!!
!!
!! EXAMPLE
!!
!! --- A single-point at a time example, does not calculate derivatives (based on Cellular Simulation)---
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES
!!  #include "Eos.h"         ! for EOS_VAR order
!!
!!  real  :: temp_zone, rho_zone, ptot, eint, gamma
!!  real, dimension(EOS_NUM)  :: eosData
!!  real, dimension(SPECIES_BEGIN:SPECIES_END) ::  massFraction  
!!  integer, dimension(2,MDIM)                 :: blockRange,blockExtent
!!
!!
!!  massFraction(:) = 1.0e-12        
!!  massFraction(C12_SPEC) = 1.0
!!
!!  .... initiale temp_zone, rho_zone
!!
!!  call Grid_getBlkIndexLimits(blockId,blockRange,blockExtent)
!!  do k = blockRange(LOW,KAXIS), blockRange(HIGH,KAXIS)
!!     do j = blockRange(LOW,JAXIS),blockRange(HIGH,JAXIS)
!!        do i = blockRange(LOW,IAXIS),blockRange(HIGH,IAXIS)
!!
!!           eosData(EOS_TEMP) = temp_zone
!!           eosData(EOS_DENS) = rho_zone
!!
!!           call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
!!
!!           ptot = eosData(EOS_PRES)
!!           eint = eosData(EOS_EINT)
!!           gamma = eosData(EOS_GAMC)
!!           
!!           call Grid_putPointData(blockId,CENTER,TEMP_VAR,EXTERIOR,iPosition,temp_zone)
!!           call Grid_putPointData(blockId,CENTER,DENS_VAR,EXTERIOR,iPosition,rho_zone)
!!           call Grid_putPointData(blockId,CENTER,PRES_VAR,EXTERIOR,iPosition,ptot)
!!           call Grid_putPointData(blockId,CENTER,EINT_VAR,EXTERIOR,iPosition,eint)
!!              if you want ENER_VAR, calculating it from EINT_VAR and kinetic energy
!!           call Grid_putPointData(blockId,CENTER,GAMC_VAR,EXTERIOR,iPosition,gamma)
!!           call Grid_putPointData(blockId,CENTER,GAME_VAR,EXTERIOR,iPosition,(ptot/(etot*sim_rhoAmbient) + 1.0))
!!
!!         enddo  ! end of k loop
!!     enddo     ! end of j loop
!!  enddo        ! end of i loop
!!
!! ------------------ Row at a time example, with derivates (based on Eos_unitTest) --------
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES, EOS_NUM
!!  #include "Eos.h"         ! for EOS_VAR order
!!  integer veclen, isize, jsize, ksize, i,j,k, e
!!  real, dimension(:), allocatable :: eosData
!!  real, dimension(:), allocatable :: massFrac
!!  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
!!  real, allocatable, dimension(:,:,:,:) :: derivedVariables
!!  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!!
!!   ! in the Eos_unitTest, this loops over all blocks.... here is a snippet from inside
!!     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!
!!    !  Allocate the necessary arrays for an entire block of data
!!    isize = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
!!    jsize = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
!!    ksize = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)
!!    vecLen=isize
!!    allocate(derivedVariables(isize,jsize,ksize,EOS_NUM))
!!    allocate(eosData(vecLen*EOS_NUM))
!!    allocate(massFrac(vecLen*NSPECIES))
!!    mask = .true.
!!
!!    ! indices into the first location for these variables
!!    pres = (EOS_PRES-1)*vecLen
!!    dens = (EOS_DENS-1)*vecLen
!!    temp = (EOS_TEMP-1)*vecLen
!!
!!
!!    call Grid_getBlkPtr(blockID,solnData)
!!    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
!!           do i = 1,vecLen
!!              massFrac((i-1)*NSPECIES+1:i*NSPECIES) = &
!!                   solnData(SPECIES_BEGIN:SPECIES_END,ib+i-1,j,k)
!!           end do
!!
!!           eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
!!           eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
!!           ! Eos Helmholtz needs a good initial estimate of temperature no matter what the mode
!!           eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
!!
!!           call Eos(MODE_DENS_PRES,vecLen,eosData,massFrac,mask)
!!
!!           do e=EOS_VARS+1,EOS_NUM
!!              m = (e-1)*vecLen
!!              derivedVariables(1:vecLen,j-NGUARD,k-NGUARD,e) =  eosData(m+1:m+vecLen)
!!           end do
!!        end do
!!     end do
!!
!! NOTES
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, and MODE_DENS_PRES are defined in constants.h.
!!
!!  All routines calling this routine should include a 
!!  use Eos_interface 
!!  statement, preferable with "ONLY" attribute e.g.
!!  use Eos_interface, ONLY:  Eos
!!
!!  For Gamma and Multigamma routines, the entropy and entropy derivatives 
!!  calculations have not been confirmed to be correct.  Use with caution.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***


subroutine Eos(mode, vecLen, eosData, massFrac, mask)

!==============================================================================
  use Eos_data, ONLY: eos_gammam1j, eos_gasConstant, eos_ggprodj, eos_gamma, eos_largee
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  real, optional, INTENT(in),dimension(vecLen*NSPECIES)    :: massFrac
  logical,  optional, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask
  
  real,dimension(NSPECIES) :: weight
  real :: rt,abarValue, abarInv, zbarValue, zbarFrac
  integer :: specieStart, specieEnd
  integer :: dens, temp, pres, eint, abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, c_v, c_p, gamc, pel, ne, eta
  integer :: i

  
!==============================================================================

  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  entr = (EOS_ENTR-1)*vecLen

  if(.not.present(massFrac)) then
     call Driver_abortFlash("[Eos] Multigamma needs mass fractions")
  end if


! mode:  temperature and density given
  if (mode == MODE_DENS_TEMP) then
     do i = 1,vecLen
        !  Add these temporary variables to prevent execution warnings under -debug
        ! Sadly, doesn't solve it -- still need to create a temporary array tempArray = massFrac(specieStart:specieEnd)        
        !  Not done due to copying slowdown

        specieStart = (i-1)*NSPECIES + 1
        specieEnd = i*NSPECIES

        call Multispecies_getSumInv(A, abarInv ,massFrac(specieStart:specieEnd))
        abarValue = 1.e0 / abarInv

        ! zbar for flame modules
        call Multispecies_getSumFrac(Z,zbarFrac,massFrac(specieStart:specieEnd))
        zbarValue = abarValue*zbarFrac

        weight = massFrac(specieStart:specieEnd)*eos_gammam1j
        call Multispecies_getSumInv(A,rt, weight)

        eosData(abar+i) = abarValue
        eosData(zbar+i) = zbarValue
        eosData(gamc+i) = 1.0e0 + 1.0e0/(rt*eosData(abar+i))
        eosData(eint+i) = eos_gasConstant*eosData(temp+i)/((eosData(gamc+i)-1.e0)*eosData(abar+i))
        if (eosData(eint+i) .gt. eos_largee) then
            eosData(eint+i) = eos_largee
            eosData(temp+i) = eosData(eint+i) *(eosData(gamc+i) - 1.0) * eosData(abar+i) / eos_gasConstant
        endif
        eosData(pres+i) = eos_gasConstant*eosData(dens+i)*eosData(temp+i)/eosData(abar+i)
        eosData(entr+i) = (eosData(pres+i)/eosData(dens+i) + eosData(eint+i))/eosData(temp+i)
     end do

! mode:  density and internal energy
  elseif (mode == MODE_DENS_EI) then

     do i = 1,vecLen 
        specieStart = (i-1)*NSPECIES + 1
        specieEnd = i*NSPECIES

        call Multispecies_getSumInv(A, abarInv, massFrac(specieStart:specieEnd))
        abarValue = 1.e0 / abarInv
       
        call Multispecies_getSumFrac(Z,zbarFrac,massFrac(specieStart:specieEnd))
        zbarValue = abarValue*zbarFrac

        weight = massFrac(specieStart:specieEnd)*eos_gammam1j
        call Multispecies_getSumInv(A, rt, weight)

        eosData(abar+i) = abarValue
        eosData(zbar+i) = zbarValue
        eosData(gamc+i) = 1.0e0 + 1.0e0/(rt*eosData(abar+i))
        if (eosData(eint+i) .gt. eos_largee) then
            eosData(eint+i) = eos_largee
        endif
        eosData(pres+i) = eosData(dens+i)*eosData(eint+i)*(eosData(gamc+i)-1.e0)
        eosData(temp+i) = eosData(eint+i)*(eosData(gamc+i)-1.e0) * eosData(abar+i)/eos_gasConstant
        eosData(entr+i) = (eosData(pres+i)/eosData(dens+i) + eosData(eint+i))/eosData(temp+i)
     end do

! mode:  density and pressure
  elseif (mode == MODE_DENS_PRES) then
     do i = 1,vecLen
        specieStart = (i-1)*NSPECIES + 1
        specieEnd = i*NSPECIES

        call Multispecies_getSumInv(A, abarInv, massFrac(specieStart:specieEnd))
        abarValue = 1.e0 / abarInv

        call Multispecies_getSumFrac(Z,zbarFrac,massFrac(specieStart:specieEnd))
        zbarValue = abarValue*zbarFrac

        weight = massFrac(specieStart:specieEnd)*eos_gammam1j
        call Multispecies_getSumInv(A, rt, weight)

        eosData(abar+i) = abarValue
        eosData(zbar+i) = zbarValue
        eosData(gamc+i) = 1.0e0 + 1.0e0 / (rt*eosData(abar+i))
        eosData(eint+i) = eosData(pres+i) / ( ( eosData(gamc+i) - 1.0 ) * eosData(dens+i) )
        if (eosData(eint+i) .gt. eos_largee) then
            eosData(eint+i) = eos_largee
            eosData(pres+i) = eosData(dens+i)*eosData(eint+i)*(eosData(gamc+i)-1.e0)
        endif
        eosData(temp+i) = eosData(eint+i) *(eosData(gamc+i) - 1.0) * eosData(abar+i) / eos_gasConstant
        eosData(entr+i) = (eosData(pres+i)/eosData(dens+i) + eosData(eint+i))/eosData(temp+i)
     end do

  else 
     call Driver_abortFlash("[Eos] Unrecognized input mode given to Eos")     
  endif
  
  if(present(mask)) then
     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        !! flash2 equation is dpt = gasconstant*density/abar
        eosData(dpt+1:dpt+vecLen) = eos_gasConstant*eosData(dens+1:dens+vecLen)/eosData(abar+1:abar+vecLen)
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
        !! flash2 equation is dpd =gasconstant*temperature/abar 
        eosData(dpd+1:dpd+vecLen) = eos_gasConstant*eosData(temp+1:temp+vecLen)/eosData(abar+1:abar+vecLen)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
         eosData(det+1:det+vecLen) = eos_gasConstant / eosData(abar+1:abar+vecLen)* & 
                (eosData(gamc+1:gamc+vecLen) - 1.0)
     end if
     if(mask(EOS_DED))then 
        ded = (EOS_DED-1)*vecLen
        eosData(ded+1:ded+vecLen) = 0.
     end if

     ! Entropy derivatives
     if (mask(EOS_DST)) then
        if (mask(EOS_DET) .AND. mask(EOS_DPT)) then
           det = (EOS_DET-1)*vecLen
           dpt = (EOS_DPT-1)*vecLen
           dst = (EOS_DST-1)*vecLen
           eosData(dst+1:dst+vecLen) = &
         &       ((eosData(dpt+1:dpt+vecLen) / eosData(dens+1:dens+vecLen) + eosData(det+1:det+vecLen)) - &
         &       (eosData(pres+1:pres+vecLen) / eosData(dens+1:dens+vecLen) + eosData(eint+1:eint+vecLen)) / &
         &       eosData(temp+1:temp+vecLen) ) / eosData(temp+1:temp+vecLen)
        else
           call Driver_abortFlash("[Eos] Cannot calculate EOS_DST without EOS_DET and EOS_DPT")
        end if
     end if
     if (mask(EOS_DSD)) then
        if (mask(EOS_DED) .AND. mask(EOS_DPD)) then
           dsd = (EOS_DSD-1)*vecLen
           ded = (EOS_DED-1)*vecLen
           dpd = (EOS_DPD-1)*vecLen
           eosData(dsd+1:dsd+vecLen) = &
        &       ( ((eosData(dpd+1:dpd+vecLen) - eosData(pres+1:pres+vecLen)/eosData(dens+1:dens+vecLen)) / &
        &          eosData(dens+1:dens+vecLen)) + eosData(ded+1:ded+vecLen)) / eosData(temp+1:temp+vecLen)
         else
           call Driver_abortFlash("[Eos] Cannot calculate EOS_DSD without EOS_DED and EOS_DPD")
        end if
     end if


     if(mask(EOS_PEL))then 
        pel = (EOS_PEL-1)*vecLen
        eosData(pel+1:pel+vecLen) = 0.
     end if
     if(mask(EOS_NE))then 
        ne = (EOS_NE-1)*vecLen
        eosData(ne+1:ne+vecLen) = 0.
     end if
     if(mask(EOS_ETA))then 
        eta = (EOS_ETA-1)*vecLen
        eosData(eta+1:eta+vecLen) = 0.
     end if
     
     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           eosData(c_v+1:c_v+vecLen) = eosData(det+1:det+vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if
     
     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           eosData(c_p+1:c_p+vecLen) = eos_gamma*eosData(c_v+1:c_v+vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
     end if
  end if
  return
end subroutine Eos

!!..no matter what the input mode compute the entropy
!!..ignore the -chemical_potential*number_density part for now
!!$  dens_inv = 1.0e0/eosData(dens+1:+vecLen)
!!$  temp_inv = 1.0e0/eosData(temp+1:+vecLen)
!!$  stot     = (pres*dens_inv + eosData(eint+1:+vecLen))*temp_inv 
!!$  dstotdd  = (eosData(EOS_DPD)*dens_inv - pres*dens_inv*dens_inv + eosData(EOS_DED))*temp_inv
!!$  dstotdt  = (eosData(EOS_DPT)*dens_inv + eosData(EOS_DET))*temp_inv  - (pres*dens_inv + eosData(eint+1:+vecLen)) * temp_inv*temp_inv 
!!$  



