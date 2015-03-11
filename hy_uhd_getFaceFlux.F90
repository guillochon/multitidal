!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_getFaceFlux
!!
!! NAME
!!
!!  hy_uhd_getFaceFlux
!!
!! SYNOPSIS
!!
!!  hy_uhd_getFaceFlux( integer(IN) :: blockID,
!!                      integer(IN) :: blkLimits(2,MDIM)
!!                      integer(IN) :: blkLimitsGC(2,MDIM), 
!!                      integer(IN) :: datasize(MDIM),
!!                      real(IN)    :: del(MDIM),
!!                      real(OUT)   :: xflux(:,:,:,:),
!!                      real(OUT)   :: yflux(:,:,:,:),
!!                      real(OUT)   :: zflux(:,:,:,:),
!!                      real, pointer, dimension(:,:,:,:) :: scrch_Ctr,
!!                      real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig,
!!                      logical,optional(IN) :: lastCall )
!!
!! ARGUMENTS
!!
!!  blockID           - a current block ID
!!  blkLimits         - an array that holds the lower and upper indices of the section
!!                      of block without the guard cells
!!  blkLimitsGC       - an array that holds the lower and upper indices of the section
!!                      of block with the guard cells
!!  datasize          - data size for boundary extrapolated data, boundary data
!!  del               - grid deltas
!!  xflux,yflux,zflux - face fluxes at each {x,y,z} direction
!!  scrch_Ctr         - Pointer to the scrch array
!!  hy_SpcR,hy_SpcL,hy_SpcSig - Pointers for Species and mass scalar recon.
!!  lastCall          - if true then store flux data in scratch array
!!
!! DESCRIPTION
!!
!!  This routine computes high-order Godunov fluxes at cell interface centers 
!!  for each spatial direction using a choice of Riemann solvers.
!!  Choices of Riemann solvers are Roe-type, HLL(E), HLLC, HLLD, 
!!  Hybrid (HLLC+HLL for hydro; HLLD+HLL for MHD), local Lax-Friedrichs, and Marquina solvers.
!!
!!*** 

!!REORDER(4):U, scrch_Ctr, scrch_Ptr, scrch_[XYZ], [xyz]flux

!#define COMPUTE_DT_FLUX

#include "Flash.h"
subroutine hy_uhd_getFaceFlux ( blockID,blkLimits,blkLimitsGC,datasize,del,&
                                xflux,yflux,zflux,scrch_Ctr,hy_SpcR,hy_SpcL,hy_SpcSig,lastCall)

  use Hydro_data,                    ONLY : hy_nref,hy_kref,          &
                                            hy_RiemannSolver,         &
                                            hy_useDiffuse,            &
                                            hy_useViscosity,          &
                                            hy_useConductivity,       &
                                            hy_use_avisc,             &
                                            hy_cvisc,                 &
                                            hy_updateHydroFluxes,     &
                                            hy_addThermalFlux,        &
                                            hy_geometry,              &
                                            hy_useAuxEintEqn,         &
                                            hy_hydroComputeDtOption

  use hy_uhd_interface,              ONLY : hy_uhd_addViscousFluxes,  &
                                            hy_uhd_addThermalFluxes,  &
                                            hy_uhd_Roe, &
                                            hy_uhd_LLF, &
                                            hy_uhd_HLL, &
                                            hy_uhd_HLLC,&
                                            hy_uhd_Marquina,&
                                            hy_uhd_MarquinaModified,&
                                            hy_uhd_setMinTimeStep

  use hy_uhd_slopeLimiters,          ONLY : signum
  use Grid_interface,                ONLY : Grid_getBlkPtr, &
                                            Grid_releaseBlkPtr, &
                                            Grid_getCellCoords
  use Eos_interface,                 ONLY : Eos
  use Conductivity_interface,        ONLY : Conductivity
  use Viscosity_interface,           ONLY : Viscosity
  use Timers_interface,              ONLY : Timers_start, Timers_stop

#if defined(FLASH_USM_MHD) 
  use Hydro_data,                    ONLY : hy_mref, &
                                            hy_forceHydroLimit,&
                                            hy_useResistivity, &
                                            hy_useMagneticResistivity,&
                                            hy_useBiermann
  use hy_uhd_interface,              ONLY : hy_uhd_addResistiveFluxes, &
                                            hy_uhd_HLLD,&
                                            hy_uhd_addBiermannBatteryTerms
  use MagneticResistivity_interface, ONLY : MagneticResistivity
#endif


  implicit none

#include "constants.h"
#include "Eos.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------------
  integer, intent(IN)  :: blockID
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM), intent(IN)         :: datasize
  real,    dimension(MDIM), intent(IN)         :: del

#ifdef FIXEDBLOCKSIZE
  real, intent(OUT) :: xflux (NFLUXES,&
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC)
  real, intent(OUT) :: yflux (NFLUXES,&
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC)
  real, intent(OUT) :: zflux (NFLUXES,&
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC)
#else
  real, intent(OUT) :: xflux (NFLUXES,&
      blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real, intent(OUT) :: yflux (NFLUXES,&
      blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real, intent(OUT) :: zflux (NFLUXES,&
      blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
#endif
  real, pointer, dimension(:,:,:,:) :: scrch_Ctr
  real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig  
  logical, optional, intent(IN) :: lastCall
  !! ---------------------------------------------------------

  integer :: i0,imax,j0,jmax,k0,kmax,jbeg,jend,kbeg,kend, i,j,k, ierr
  real, pointer, dimension(:,:,:,:) :: scrch_Ptr, U
  !U contains (DENS,VELX,VELY,VELZ,(MAGX,MAGY,MAGZ),PRES + GAMC,GAME,EINT,TEMP)
  real, dimension(HY_VARINUMMAX) :: VL,VR 
  real, dimension(NSPECIES)    :: speciesArr

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC) &
                  :: viscDynamic,viscKinematic,cond,dcff

#if defined(FLASH_USM_MHD) 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC) &
                  :: magResist
#endif

#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) &
                  :: viscDynamic,viscKinematic,cond,dcff

#if defined(FLASH_USM_MHD) 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) &
                  :: magResist
#endif

#endif
  real    :: cvisc
  integer :: k2,k3,kGrav,kUSM
  real    :: presPlus, presMinus

#ifdef FLASH_UHD_HYDRO
  !CD: Add these variables so we can maintain a single omp parallel
  !directive for unsplit hydro and MHD simulations.  They are not used.
  real, dimension(1,1,1) :: magResist
  real, save :: hy_mref = 0
  logical, save :: hy_useResistivity = .false.
  logical, save :: hy_useMagneticResistivity = .false.
  logical, save :: hy_useBiermann = .false.
#endif

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCenter, xLeft, xRight
  real, dimension(GRID_JHI_GC) :: yCenter
#else
  real, dimension(datasize(IAXIS)) :: xCenter, xLeft, xRight
  real, dimension(datasize(JAXIS)) :: yCenter
#endif
  real :: dPdr, rvol, alpha
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: speed, dy, dz
  integer :: ispu,isph
  integer, dimension(MDIM) :: normDir,tranDir


  xflux = 0.
  if (NDIM < 3) then
     zflux = 0.
     if (NDIM < 2) yflux = 0.
  endif

  kGrav=0
#ifdef GRAVITY
  kGrav=1
#endif

  kUSM = 0
#ifdef FLASH_USM_MHD
  if (.NOT. hy_forceHydroLimit) kUSM = 1
#endif


  call Grid_getBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR) 
  call Grid_getBlkPtr(blockID,U,CENTER)

  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = blkLimits(LOW, KAXIS)
  kmax = blkLimits(HIGH,KAXIS)

  if (NDIM == 1) then
     jbeg= 3
     jend=-1
     kbeg= 3
     kend=-1
     k2 = 0
     k3 = 0
  elseif (NDIM == 2) then
     jbeg= j0
     jend= jmax
     kbeg= 3
     kend=-1
     k2 = 1
     k3 = 0
  else
     jbeg= j0
     jend= jmax
     kbeg= k0
     kend= kmax
     k2 = 1
     k3 = 1
  endif


  if (hy_useDiffuse) then
     ! call Timers_start("get diffusion")
     !! Initialize

     do k=kbeg-2,kend+2
        do j=jbeg-2,jend+2
           do i=i0-2,imax+2

              !! copy species to a temporary array
              speciesArr(:) = U(SPECIES_BEGIN:SPECIES_END,i,j,k)

              if (hy_useViscosity) then
                 !! Get viscosity
                 call Viscosity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                      speciesArr,viscDynamic(i,j,k),viscKinematic(i,j,k))
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 !! Get heat conductivity
                 call Conductivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                      speciesArr,cond(i,j,k),dcff(i,j,k),2)
              endif
#if defined(FLASH_USM_MHD) 
              if (hy_useMagneticResistivity) then
                 !! Get magnetic viscosity
                 !!call MagneticResistivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                 !!     speciesArr,magResist(i,j,k))
                 call MagneticResistivity(U(:,i,j,k),magResist(i,j,k))

              endif
#endif
           enddo
        enddo
     enddo

     !! For non-ideal fluxes
#if defined(FLASH_USM_MHD) 
     magResist  = magResist/hy_mref
#endif
     viscDynamic= viscDynamic/hy_nref
  endif

  if (hy_geometry /= CARTESIAN) then
     call Grid_getCellCoords(IAXIS,blockID, CENTER,    .true.,xCenter, dataSize(IAXIS))
     call Grid_getCellCoords(JAXIS,blockID, CENTER,    .true.,yCenter, dataSize(JAXIS))
     call Grid_getCellCoords(IAXIS,blockID, LEFT_EDGE, .true.,xLeft,   dataSize(IAXIS))
     call Grid_getCellCoords(IAXIS,blockID, RIGHT_EDGE,.true.,xRight,  dataSize(IAXIS))
  endif

  !! Compute intercell fluxes using the updated left & right states
  !! Calculate x-flux first
  normDir=0; normDir(DIR_X)=1
  tranDir=2; tranDir(DIR_X)=0
  tranDir=tranDir*kUSM
  tranDir(2)=tranDir(2)*k2
  tranDir(3)=tranDir(3)*k3

  do k=blkLimits(LOW,KAXIS)-tranDir(DIR_Z),blkLimits(HIGH,KAXIS)+tranDir(DIR_Z)+normDir(DIR_Z)
     do j=blkLimits(LOW,JAXIS)-tranDir(DIR_Y),blkLimits(HIGH,JAXIS)+tranDir(DIR_Y)+normDir(DIR_Y)
        do i=blkLimits(LOW,IAXIS)-tranDir(DIR_X),blkLimits(HIGH,IAXIS)+tranDir(DIR_X)+normDir(DIR_X)

           if (hy_updateHydroFluxes) then

              VL(HY_DENS:HY_END_VARS-kGrav)=scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i-1,j,k)
              VR(HY_DENS:HY_END_VARS-kGrav)=scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,  j,k)

#ifdef BDRY_VAR
              ! solid internal boundary
              ! Cell i and i-1:
              if (U(BDRY_VAR,i,j,k) > 0.0 .and. U(BDRY_VAR,i-1,j,k) < 0.0) then
                 VR(HY_DENS:HY_END_VARS-kGrav) = VL(HY_DENS:HY_END_VARS-kGrav)
                 VR(HY_VELX) = -VL(HY_VELX)
              end if

              if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i-1,j,k) > 0.0) then
                 VL(HY_DENS:HY_END_VARS-kGrav) = VR(HY_DENS:HY_END_VARS-kGrav)
                 VL(HY_VELX) = -VR(HY_VELX)
              end if
#endif

              if (hy_RiemannSolver == ROE) then
                 call hy_uhd_Roe(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_Roe", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == HLL) then
                 call hy_uhd_HLL(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == HLLC) then
                 call hy_uhd_HLLC(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_X)

#if defined(FLASH_USM_MHD) 
              elseif (hy_RiemannSolver == HLLD) then
                 call hy_uhd_HLLD(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_X)
#endif
              elseif (hy_RiemannSolver == MARQ) then
                 call hy_uhd_Marquina(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == MARM) then
                 call hy_uhd_MarquinaModified(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == LLF) then
                 call hy_uhd_LLF(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_LLF", VL, VR, i,j,k, DIR_X)

              elseif (hy_RiemannSolver == HYBR) then
                 if (U(SHOK_VAR,i-1,j,k) + U(SHOK_VAR,i,j,k) > 0.) then
                    !! use diffusive HLL solver for a local strong shock/rarefaction region
                    call hy_uhd_HLL(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                    if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_X)
                 else
#ifdef FLASH_UHD_HYDRO
                    ! Hydro
                    call hy_uhd_HLLC(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                    if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_X)
#else
                    ! MHD (with change suggested by DW, JFG)
                    call hy_uhd_HLLC(DIR_X,VL,VR,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                    if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_X)
#endif
                 endif
              endif

              if (hy_hydroComputeDtOption == 1) then
                 !! Call for dt calculation
                 call hy_uhd_setMinTimeStep(blockID,i,j,k,del(DIR_X),speed)
              endif

              !! Artificial viscosity as in PPM, Colella and Woodward, 1984.
#ifdef BDRY_VAR
              if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i-1,j,k) < 0.0) then
#endif
              if (hy_use_avisc) then
                 if (NDIM == 1) then
                    cvisc = hy_cvisc*max(    -(U(VELX_VAR,i,  j,  k  )-U(VELX_VAR,i-1,j,  k  )),0.)
                 elseif (NDIM == 2) then
                    cvisc = hy_cvisc*max(-(    U(VELX_VAR,i,  j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                         0.25*(U(VELY_VAR,i,  j+1,k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                               U(VELY_VAR,i-1,j+1,k  )-U(VELY_VAR,i-1,j-1,k  ))*del(DIR_X)/del(DIR_Y)),&
                                         0.)
                 else
                    cvisc = hy_cvisc*max(-(     U(VELX_VAR,i,  j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                         0.25*( U(VELY_VAR,i,  j+1,k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                                U(VELY_VAR,i-1,j+1,k  )-U(VELY_VAR,i-1,j-1,k  ))*del(DIR_X)/del(DIR_Y) + &
                                         0.25*( U(VELZ_VAR,i,  j,  k+1)-U(VELZ_VAR,i,  j,  k-1) + &
                                                U(VELZ_VAR,i-1,j,  k+1)-U(VELZ_VAR,i-1,j,  k-1))*del(DIR_X)/del(DIR_Z)), &
                                         0.)
                 endif
              
                 xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k) = &
                      xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k) &
                      +cvisc*(/U(DENS_VAR,i-1,j,k)                    -U(DENS_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(VELX_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(VELX_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(VELY_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(VELY_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(VELZ_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(VELZ_VAR,i,j,k)&
                              ,U(DENS_VAR,i-1,j,k)*U(ENER_VAR,i-1,j,k)-U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)&
#if defined(FLASH_USM_MHD) 
                              ,0.&
                              ,U(MAGY_VAR,i-1,j,k)                    -U(MAGY_VAR,i,j,k)&
                              ,U(MAGZ_VAR,i-1,j,k)                    -U(MAGZ_VAR,i,j,k)&
#endif
                             /)
              endif
#ifdef BDRY_VAR
              endif
#endif
              !! Flux for internal energy density
              !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
              !! Note that there is a typo in Li's paper related on the left and right states.
              if (hy_useAuxEintEqn) then
                 if (xflux(HY_DENS_FLUX,i,j,k) > 0.) then
                    xflux(HY_PRES_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)/VL(HY_DENS)
                    xflux(HY_EINT_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VL(HY_EINT)
                 else
                    xflux(HY_PRES_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)/VR(HY_DENS)
                    xflux(HY_EINT_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VR(HY_EINT)
                 endif

              endif

#ifdef FLASH_UHD_3T
              if (xflux(HY_DENS_FLUX,i,j,k) > 0.) then
                 xflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VL(HY_EELE:HY_ERAD)
              else
                 xflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*VR(HY_EELE:HY_ERAD)
              endif

#endif

#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
              do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
                 isph= ispu-NPROP_VARS
                 if (xflux(HY_DENS_FLUX,i,j,k) < 0.) then
                    xflux(HY_END_FLUX+isph,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*hy_SpcL(isph,i,  j,k,DIR_X)
                 else
                    xflux(HY_END_FLUX+isph,i,j,k) = xflux(HY_DENS_FLUX,i,j,k)*hy_SpcR(isph,i-1,j,k,DIR_X)
                 endif
              enddo
#endif
#endif


           endif ! end of if (hy_updateHydroFluxes) then

           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes&
                      (blockID,blkLimitsGC,i,j,k,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),viscDynamic,DIR_X)
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 call hy_uhd_addThermalFluxes&
                      (blockID,blkLimitsGC,i,j,k,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),cond,DIR_X)
              endif

#if defined(FLASH_USM_MHD) 
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes&
                      (blockID,blkLimitsGC,i,j,k,xflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),magResist,DIR_X)
              endif
#endif
           endif

        enddo
     enddo
  enddo


!#if defined(FLASH_USM_MHD) 
#ifdef FLASH_USM_MHD
  if (present(lastCall)) then
     if (lastCall) then
#endif
        if (hy_useAuxEintEqn .OR. hy_geometry /= CARTESIAN) then
           !! Obtain an averaged pressure for internal energy update in hy_uhd_unsplitUpdate
           do k=blkLimits(LOW,KAXIS)-tranDir(DIR_Z),blkLimits(HIGH,KAXIS)+tranDir(DIR_Z)+normDir(DIR_Z)
              do j=blkLimits(LOW,JAXIS)-tranDir(DIR_Y),blkLimits(HIGH,JAXIS)+tranDir(DIR_Y)+normDir(DIR_Y)
                 do i=blkLimits(LOW,IAXIS)-tranDir(DIR_X),blkLimits(HIGH,IAXIS)+tranDir(DIR_X)+normDir(DIR_X)

                    !! For non-cartesian geometries
                    if (hy_geometry /= CARTESIAN) then
                       select case(hy_geometry) ! First, select whether y or z is phi-direction
                       case(CYLINDRICAL,POLAR)
                          alpha = 1.
                       case(SPHERICAL)
                          alpha = 2.
                       end select
                    endif !end of non-Cartesian support

                    if (hy_geometry /= CARTESIAN) then

                       presPlus = scrch_Ctr(XP05_SCRATCH_CENTER_VAR,i,j,k)
                       presMinus =scrch_Ctr(XN05_SCRATCH_CENTER_VAR,i,j,k)

#if defined(FLASH_USM_MHD) 
                       presPlus = presPlus + &
                            0.5*(scrch_Ctr(XP06_SCRATCH_CENTER_VAR,i,j,k)**2+&
                            scrch_Ctr(XP07_SCRATCH_CENTER_VAR,i,j,k)**2+&
                            scrch_Ctr(XP08_SCRATCH_CENTER_VAR,i,j,k)**2)

                       presMinus= presMinus + &
                            0.5*(scrch_Ctr(XN06_SCRATCH_CENTER_VAR,i,j,k)**2+&
                            scrch_Ctr(XN07_SCRATCH_CENTER_VAR,i,j,k)**2+&
                            scrch_Ctr(XN08_SCRATCH_CENTER_VAR,i,j,k)**2)
#endif
                       scrch_Ptr(VAR2_SCRATCH_CENTER_VAR,i,j,k) = &
                            (xLeft(i)*presMinus + xRight(i)*presPlus)*0.5/xCenter(i)

                    else
                       scrch_Ptr(VAR2_SCRATCH_CENTER_VAR,i,j,k) = &
                            0.5*( scrch_Ctr(XN05_SCRATCH_CENTER_VAR,i,j,k) &
                                 +scrch_Ctr(XP05_SCRATCH_CENTER_VAR,i,j,k))
                    endif

                 enddo
              enddo
           enddo
        endif ! end of if (hy_useAuxEintEqn)

        !#ifdef FLASH_USM_MHD
#if defined(FLASH_USM_MHD) || defined(FLASH_UHD_NEED_SCRATCHVARS)
        !keeping the above line as an option to use scrch arrays for flux conservation
#ifndef FLASH_GRID_UG
        scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,&
                                          i0:imax+1,jbeg-2:jend+2,kbeg-2:kend+2) &
         = xflux(HY_DENS_FLUX:HY_END_FLUX,i0:imax+1,jbeg-2:jend+2,kbeg-2:kend+2)
#endif
        !#endif
#endif


!#if defined(FLASH_USM_MHD) 
#ifdef FLASH_USM_MHD
     endif !end of if (lastCall)
  endif !end of if-present(lastCall)
#endif



#if NDIM >= 2

  !! Calculate y-flux
  normDir=0; normDir(DIR_Y)=1
  tranDir=2; tranDir(DIR_Y)=0
  tranDir=tranDir*kUSM
  tranDir(2)=tranDir(2)*k2
  tranDir(3)=tranDir(3)*k3

  do k=blkLimits(LOW,KAXIS)-tranDir(DIR_Z),blkLimits(HIGH,KAXIS)+tranDir(DIR_Z)+normDir(DIR_Z)
     do j=blkLimits(LOW,JAXIS)-tranDir(DIR_Y),blkLimits(HIGH,JAXIS)+tranDir(DIR_Y)+normDir(DIR_Y)
        do i=blkLimits(LOW,IAXIS)-tranDir(DIR_X),blkLimits(HIGH,IAXIS)+tranDir(DIR_X)+normDir(DIR_X)

        if (hy_updateHydroFluxes) then

           VL(HY_DENS:HY_END_VARS-kGrav)=scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j-1,k)
           VR(HY_DENS:HY_END_VARS-kGrav)=scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,  k)

#ifdef BDRY_VAR
           ! solid internal boundary
           if (U(BDRY_VAR,i,j,k) > 0.0 .and. U(BDRY_VAR,i,j-1,k) < 0.0) then
              VR(HY_DENS:HY_END_VARS-kGrav) = VL(HY_DENS:HY_END_VARS-kGrav)
              VR(HY_VELY) = -VL(HY_VELY)
           end if

           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j-1,k) > 0.0) then
              VL(HY_DENS:HY_END_VARS-kGrav) = VR(HY_DENS:HY_END_VARS-kGrav)
              VL(HY_VELY) = -VR(HY_VELY)
           end if
#endif

           if (hy_RiemannSolver == ROE) then
              call hy_uhd_Roe(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Roe", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == HLL) then
              call hy_uhd_HLL(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == HLLC) then
              call hy_uhd_HLLC(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Y)

#if defined(FLASH_USM_MHD) 
           elseif (hy_RiemannSolver == HLLD) then
              call hy_uhd_HLLD(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Y)
#endif
           elseif (hy_RiemannSolver == MARQ) then
              call hy_uhd_Marquina(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == MARM) then
              call hy_uhd_MarquinaModified(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == LLF) then
              call hy_uhd_LLF(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_LLF", VL, VR, i,j,k, DIR_Y)

           elseif (hy_RiemannSolver == HYBR) then
              if (U(SHOK_VAR,i,j-1,k) + U(SHOK_VAR,i,j,k) > 0.) then
                 !! use diffusive HLL solver for a local strong shock/rarefaction region
                 call hy_uhd_HLL(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Y)
              else
#ifdef FLASH_UHD_HYDRO
                 ! for hydro
                 call hy_uhd_HLLC(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Y)
#else
                 ! MHD (with change suggested by DW, JFG)
                 call hy_uhd_HLLC(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Y)
#endif
              endif
           endif


           if (hy_hydroComputeDtOption == 1) then
              !! Call for dt calculation
              if (hy_geometry == CARTESIAN .OR. hy_geometry == CYLINDRICAL) then 
                 ! the 'y' coordinate is not angular in Cartesian and cylindrical coords
                 dy = del(DIR_Y)
              else ! Angular coordinates in 2D: Spherical or Polar
                 dy = xCenter(i)*del(DIR_Y)
              endif
              call hy_uhd_setMinTimeStep(blockID,i,j,k,dy,speed)
           endif


           !! Artificial viscosity as in PPM, Colella and Woodward, 1984.
#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j-1,k) < 0.0) then
#endif
           if (hy_use_avisc) then
              if (NDIM == 2) then
                 cvisc = hy_cvisc*max(-(0.25*( U(VELX_VAR,i+1,j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                               U(VELX_VAR,i+1,j-1,k  )-U(VELX_VAR,i-1,j-1,k  ))*del(DIR_Y)/del(DIR_X) + &
                                               U(VELY_VAR,i,  j,  k  )-U(VELY_VAR,i,  j-1,k  )), &
                                      0.)
              else
                 cvisc = hy_cvisc*max(-(0.25*( U(VELX_VAR,i+1,j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                               U(VELX_VAR,i+1,j-1,k  )-U(VELX_VAR,i-1,j-1,k  ))*del(DIR_Y)/del(DIR_X) + &
                                               U(VELY_VAR,i,  j,  k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                        0.25*( U(VELZ_VAR,i,  j,  k+1)-U(VELZ_VAR,i,  j,  k-1) + &
                                               U(VELZ_VAR,i,  j-1,k+1)-U(VELZ_VAR,i,  j-1,k-1))*del(DIR_Y)/del(DIR_Z)), &
                                      0.)

              endif

              yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k) = &
                   yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k) &
                   +cvisc*(/U(DENS_VAR,i,j-1,k)                    -U(DENS_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(VELX_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(VELX_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(VELY_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(VELY_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(VELZ_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(VELZ_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j-1,k)*U(ENER_VAR,i,j-1,k)-U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)&
#if defined(FLASH_USM_MHD) 
                           ,U(MAGX_VAR,i,j-1,k)                    -U(MAGX_VAR,i,j,k)&
                           ,0.&
                           ,U(MAGZ_VAR,i,j-1,k)                    -U(MAGZ_VAR,i,j,k)&
#endif
                           /)

           endif
#ifdef BDRY_VAR
           endif
#endif

           if (hy_useAuxEintEqn) then
              !! Flux for internal energy density
              !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
              if (yflux(HY_DENS_FLUX,i,j,k) > 0.) then
                 yflux(HY_PRES_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)/VL(HY_DENS)
                 yflux(HY_EINT_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VL(HY_EINT)
              else
                 yflux(HY_PRES_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)/VR(HY_DENS)
                 yflux(HY_EINT_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VR(HY_EINT)
              endif
           endif

#ifdef FLASH_UHD_3T
           if (yflux(HY_DENS_FLUX,i,j,k) > 0.) then
              yflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VL(HY_EELE:HY_ERAD)
           else
              yflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*VR(HY_EELE:HY_ERAD)
           endif
#endif

#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
           do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
              isph= ispu-NPROP_VARS
              if (yflux(HY_DENS_FLUX,i,j,k) < 0.) then
                 yflux(HY_END_FLUX+isph,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*hy_SpcL(isph,i,j,  k,DIR_Y)
              else
                 yflux(HY_END_FLUX+isph,i,j,k) = yflux(HY_DENS_FLUX,i,j,k)*hy_SpcR(isph,i,j-1,k,DIR_Y)
              endif
           enddo
#endif
#endif

        endif !end of if (hy_updateHydroFluxes) then

           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes&
                      (blockID,blkLimitsGC,i,j,k,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),viscDynamic,DIR_Y)
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 call hy_uhd_addThermalFluxes&
                      (blockID,blkLimitsGC,i,j,k,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),cond,DIR_Y)
              endif

#if defined(FLASH_USM_MHD) 
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes&
                      (blockID,blkLimitsGC,i,j,k,yflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),magResist,DIR_Y)
              endif
#endif
           endif

        enddo
     enddo
  enddo



!#if defined(FLASH_USM_MHD) 
#ifdef FLASH_USM_MHD
  if (present(lastCall)) then
     if (lastCall) then
#endif
        if (hy_useAuxEintEqn) then
           !! Average out pressures at j+1/2 (YP05) and j-1/2 (YN05) and store them to
           !! scrch_Ctr(YN05_SCRATCH_CENTER_VAR,i,j,k) before the YP* array is wiped out and
           !! reused for storing fluxes.
           scrch_Ctr(YN05_SCRATCH_CENTER_VAR,:,:,:) = &
                (scrch_Ctr(YN05_SCRATCH_CENTER_VAR,:,:,:) +  scrch_Ctr(YP05_SCRATCH_CENTER_VAR,:,:,:))*0.5
        endif !end of if (hy_useAuxEintEqn) then
!#ifdef FLASH_USM_MHD
#if defined(FLASH_USM_MHD) || defined(FLASH_UHD_NEED_SCRATCHVARS)
!keeping the above line as an option to use scrch arrays for flux conservation
#ifndef FLASH_GRID_UG
        !scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:) = 0.
        !! re-use YP array for storing fluxes
        scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,&
                                          i0-2:imax+2,jbeg:jend+1,kbeg-2:kend+2)   &
         = yflux(HY_DENS_FLUX:HY_END_FLUX,i0-2:imax+2,jbeg:jend+1,kbeg-2:kend+2)
#endif
!#endif
#endif
                
!#if defined(FLASH_USM_MHD) 
#ifdef FLASH_USM_MHD
     endif !end of if (lastCall)
  endif !end of if-present(lastCall)
#endif


#if NDIM == 3
  !! Calculate z-flux
  normDir=0; normDir(DIR_Z)=1
  tranDir=2; tranDir(DIR_Z)=0
  tranDir=tranDir*kUSM
  tranDir(2)=tranDir(2)*k2
  tranDir(3)=tranDir(3)*k3

  do k=blkLimits(LOW,KAXIS)-tranDir(DIR_Z),blkLimits(HIGH,KAXIS)+tranDir(DIR_Z)+normDir(DIR_Z)
     do j=blkLimits(LOW,JAXIS)-tranDir(DIR_Y),blkLimits(HIGH,JAXIS)+tranDir(DIR_Y)+normDir(DIR_Y)
        do i=blkLimits(LOW,IAXIS)-tranDir(DIR_X),blkLimits(HIGH,IAXIS)+tranDir(DIR_X)+normDir(DIR_X)

        if (hy_updateHydroFluxes) then

           VL(HY_DENS:HY_END_VARS-kGrav)=scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k-1)
           VR(HY_DENS:HY_END_VARS-kGrav)=scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)

#ifdef BDRY_VAR
           ! solid internal boundary
           if (U(BDRY_VAR,i,j,k) > 0.0 .and. U(BDRY_VAR,i,j,k-1) < 0.0) then
              VR(HY_DENS:HY_END_VARS-kGrav) = VL(HY_DENS:HY_END_VARS-kGrav)
              VR(HY_VELZ) = -VL(HY_VELZ)
           end if

           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j,k-1) > 0.0) then
              VL(HY_DENS:HY_END_VARS-kGrav) = VR(HY_DENS:HY_END_VARS-kGrav)
              VL(HY_VELZ) = -VR(HY_VELZ)
           end if
#endif

           if (hy_RiemannSolver == ROE) then
              call hy_uhd_Roe(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Roe", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == HLL) then
              call hy_uhd_HLL(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == HLLC) then
              call hy_uhd_HLLC(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Z)

#if defined(FLASH_USM_MHD) 
           elseif (hy_RiemannSolver == HLLD) then
              call hy_uhd_HLLD(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Z)
#endif
           elseif (hy_RiemannSolver == MARQ) then
              call hy_uhd_Marquina(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == MARM) then
              call hy_uhd_MarquinaModified(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_Marquina", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == LLF) then
              call hy_uhd_LLF(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
              if(ierr /= 0) call do_error("hy_uhd_LLF", VL, VR, i,j,k, DIR_Z)

           elseif (hy_RiemannSolver == HYBR) then
              if (U(SHOK_VAR,i,j,k-1) + U(SHOK_VAR,i,j,k) > 0.) then
                 !! use diffusive HLL solver for a local strong shock/rarefaction region
                 call hy_uhd_HLL(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLL", VL, VR, i,j,k, DIR_Z)
              else
#ifdef FLASH_UHD_HYDRO
                 ! for hydro
                 call hy_uhd_HLLC(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLC", VL, VR, i,j,k, DIR_Z)
#else
                 ! MHD (with change suggested by DW, JFG)
                 call hy_uhd_HLLC(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),speed,ierr)
                 if(ierr /= 0) call do_error("hy_uhd_HLLD", VL, VR, i,j,k, DIR_Z)
#endif
              endif
           endif
 
           if (hy_hydroComputeDtOption == 1) then
              !! Call for dt calculation
              if (hy_geometry == CARTESIAN) then 
                 ! the 'y' coordinate is not angular in Cartesian and cylindrical coords
                 dz = del(DIR_Z)
              elseif (hy_geometry == CYLINDRICAL) then
                 dz = xCenter(i)*del(DIR_Z)
              else ! Angular coordinates in 3D: Spherical or Polar
                 dz = xCenter(i)*sin(yCenter(j))*del(DIR_Z) ! z is phi
              endif
              call hy_uhd_setMinTimeStep(blockID,i,j,k,dz,speed)
           endif

           !! Artificial viscosity as in PPM, Colella and Woodward, 1984.
#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) < 0.0 .and. U(BDRY_VAR,i,j,k-1) < 0.0) then
#endif
           if (hy_use_avisc) then
              cvisc = hy_cvisc*max(-(0.25*( U(VELX_VAR,i+1,j,  k  )-U(VELX_VAR,i-1,j,  k  ) + &
                                            U(VELX_VAR,i+1,j,  k-1)-U(VELX_VAR,i-1,j,  k-1))*del(DIR_Z)/del(DIR_X) + &
                                     0.25*( U(VELY_VAR,i,  j+1,k  )-U(VELY_VAR,i,  j-1,k  ) + &
                                            U(VELY_VAR,i,  j+1,k-1)-U(VELY_VAR,i,  j-1,k-1))*del(DIR_Z)/del(DIR_Y) + &
                                            U(VELZ_VAR,i,  j,  k  )-U(VELZ_VAR,i,  j,  k-1)),&
                                   0.)

              zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k) = &
                   zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k) &
                   +cvisc*(/U(DENS_VAR,i,j,k-1)                    -U(DENS_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(VELX_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(VELX_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(VELY_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(VELY_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(VELZ_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(VELZ_VAR,i,j,k)&
                           ,U(DENS_VAR,i,j,k-1)*U(ENER_VAR,i,j,k-1)-U(DENS_VAR,i,j,k)*U(ENER_VAR,i,j,k)&
#if defined(FLASH_USM_MHD) 
                           ,U(MAGX_VAR,i,j,k-1)                    -U(MAGX_VAR,i,j,k)&
                           ,U(MAGY_VAR,i,j,k-1)                    -U(MAGY_VAR,i,j,k)&
                           ,0. &
#endif
                           /)
           endif
#ifdef BDRY_VAR
           endif
#endif

           if (hy_useAuxEintEqn) then
              !! Flux for internal energy density
              !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
              if (zflux(HY_DENS_FLUX,i,j,k) > 0.) then
                 zflux(HY_PRES_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)/VL(HY_DENS)
                 zflux(HY_EINT_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VL(HY_EINT)
              else
                 zflux(HY_PRES_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)/VR(HY_DENS)
                 zflux(HY_EINT_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VR(HY_EINT)
              endif
           endif

#ifdef FLASH_UHD_3T
           if (zflux(HY_DENS_FLUX,i,j,k) > 0.) then
              zflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VL(HY_EELE:HY_ERAD)
           else
              zflux(HY_EELE_FLUX:HY_ERAD_FLUX,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*VR(HY_EELE:HY_ERAD)
           endif
#endif


#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
           do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
              isph= ispu-NPROP_VARS
              if (zflux(HY_DENS_FLUX,i,j,k) < 0.) then
                 zflux(HY_END_FLUX+isph,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*hy_SpcL(isph,i,j,k,  DIR_Z)
              else
                 zflux(HY_END_FLUX+isph,i,j,k) = zflux(HY_DENS_FLUX,i,j,k)*hy_SpcR(isph,i,j,k-1,DIR_Z)
              endif
           enddo
#endif
#endif

        endif !end of if (hy_updateHydroFluxes) then

           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes&
                      (blockID,blkLimitsGC,i,j,k,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),viscDynamic,DIR_Z)
              endif

              if (hy_useConductivity .and. hy_addThermalFlux) then
                 call hy_uhd_addThermalFluxes&
                      (blockID,blkLimitsGC,i,j,k,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),cond,DIR_Z)
              endif

#if defined(FLASH_USM_MHD) 
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes&
                      (blockID,blkLimitsGC,i,j,k,zflux(F01DENS_FLUX:F01DENS_FLUX+HY_VARINUM-1,i,j,k),magResist,DIR_Z)
              endif
#endif
           endif

        enddo
     enddo
  enddo


!#if defined(FLASH_USM_MHD) 
#ifdef FLASH_USM_MHD
  if (present(lastCall)) then
     if (lastCall) then
#endif
        if (hy_useAuxEintEqn) then
           !! Average out pressures at k+1/2 (ZP05) and k-1/2 (ZN05) and store them to
           !! scrch_Ctr(ZN05_SCRATCH_CENTER_VAR,i,j,k) before the ZP* array is wiped out and
           !! reused for storing fluxes.
           scrch_Ctr(ZN05_SCRATCH_CENTER_VAR,:,:,:) = &
                (scrch_Ctr(ZN05_SCRATCH_CENTER_VAR,:,:,:) +  scrch_Ctr(ZP05_SCRATCH_CENTER_VAR,:,:,:))*0.5
        endif !end of if (hy_useAuxEintEqn)

!#ifdef FLASH_USM_MHD
#if defined(FLASH_USM_MHD) || defined(FLASH_UHD_NEED_SCRATCHVARS)
!keeping the above line as an option to use scrch arrays for flux conservation
#ifndef FLASH_GRID_UG
        !! initialize with zero
        !scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:) = 0.
        !! re-use ZP array for storing fluxes
        scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,&
                                          i0-2:imax+2,jbeg-2:jend+2,kbeg:kend+1)   &
         = zflux(HY_DENS_FLUX:HY_END_FLUX,i0-2:imax+2,jbeg-2:jend+2,kbeg:kend+1)
#endif
!#endif
#endif

      
!#if defined(FLASH_USM_MHD) 
#ifdef FLASH_USM_MHD
     endif !end of if (lastCall)
  endif !end of if-present(lastCall)
#endif


#endif /* end of #if NDIM==3 */
#endif /* end of #if NDIM >= 2 */


  if (hy_useAuxEintEqn) then
  !! Average out the n+1/2 pressures at the cell interfaces to the cell center on each cell
  !! and store them to scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,:,:,:). 
  !! These averaged pressures at n+1/2 are used later for the internal energy update.
  !! Note that scrch_Ptr(VAR2_SCRATCH_CENTER_VAR,:,:,:) holds the averaged pressures 
  !! in radial R-direction only for non-Cartesian cases for geometric source term 
  !! and should not be modified here.

  !! (1) store the averaged pressure in 'R-direction' to VAR2, and
  !! (2) store the averaged pressure in 'ALL directions' to VAR1
  scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,:,:,:) = scrch_Ptr(VAR2_SCRATCH_CENTER_VAR,:,:,:)
#if NDIM > 1
  scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,:,:,:) = &
       (scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,:,:,:)+scrch_Ctr(YN05_SCRATCH_CENTER_VAR,:,:,:))*0.5
#if NDIM > 2
  scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,:,:,:) = &
    (2.*scrch_Ptr(VAR1_SCRATCH_CENTER_VAR,:,:,:)+scrch_Ctr(ZN05_SCRATCH_CENTER_VAR,:,:,:))/3.
#endif
#endif
  endif ! end of if (hy_useAuxEintEqn) then
  
  !! Release pointer
  call Grid_releaseBlkPtr(blockID,scrch_Ptr,SCRATCH_CTR)
  call Grid_releaseBlkPtr(blockID,U,CENTER)


contains

  subroutine do_error(solver, VL, VR, i,j,k, dir)
    use Grid_interface, ONLY: Grid_getCellCoords
    use Driver_interface, ONLY: Driver_abortFlash
    implicit none
    
    character(len=*), intent(in) :: solver
    real, intent(IN) :: VL(:)
    real, intent(IN) :: VR(:)
    integer, intent(IN) :: i, j, k, dir

    print *, "LEFT STATE:"
    print *, "DENS: ", VL(HY_DENS)
    print *, "PRES: ", VL(HY_PRES)
    print *, "GAMC: ", VL(HY_GAMC)
    print *, "VELX: ", VL(HY_VELX)
    print *, "VELY: ", VL(HY_VELY)
    print *, "VELZ: ", VL(HY_VELZ)

    print *, ""
    print *, "RIGHT STATE:"
    print *, "DENS: ", VR(HY_DENS)
    print *, "PRES: ", VR(HY_PRES)
    print *, "GAMC: ", VR(HY_GAMC)
    print *, "VELX: ", VR(HY_VELX)
    print *, "VELY: ", VR(HY_VELY)
    print *, "VELZ: ", VR(HY_VELZ)

    print *, ""
    print *, "INTERFACE: ", i, j, k

    ! Note: We only allocate these arrays here, which is fine because
    !       the code is to be aborted immediately once this routine is called.
    allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
    allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
    allocate(zcent(blkLimitsGC(HIGH, KAXIS)))

    call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcent, blkLimitsGC(HIGH, IAXIS)) 
    call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycent, blkLimitsGC(HIGH, JAXIS))
    call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcent, blkLimitsGC(HIGH, KAXIS))

    print *, "NEIGBORING CELLS:"

    if(dir == DIR_X) then
       print *, "X DIRECTION"
       call write_cell_info(i-3,j,k)
       call write_cell_info(i-2,j,k)
       call write_cell_info(i-1,j,k)
       call write_cell_info(i+0,j,k)
       call write_cell_info(i+1,j,k)
       call write_cell_info(i+2,j,k)
    elseif(dir == DIR_Y) then
       print *, "Y DIRECTION"
       call write_cell_info(i,j-3,k)
       call write_cell_info(i,j-2,k)
       call write_cell_info(i,j-1,k)
       call write_cell_info(i,j+0,k)
       call write_cell_info(i,j+1,k)
       call write_cell_info(i,j+2,k)
    else
       print *, "Z DIRECTION"
       call write_cell_info(i,j,k-3)
       call write_cell_info(i,j,k-2)
       call write_cell_info(i,j,k-1)
       call write_cell_info(i,j,k+0)
       call write_cell_info(i,j,k+1)
       call write_cell_info(i,j,k+2)
    end if

    deallocate(xcent)
    deallocate(ycent)
    deallocate(zcent)


    call Driver_abortFlash( &
         "[hy_uhd_getFaceFlux]: Imaginary sound speed has obtained in " // &
         trim(solver) // " solver. " // &
         "Please try other (more diffusive) slope limiter, flux, order, cfl, etc. "//&
         "in order to increase numerical stability. LOOK AT THE LOG FILE")

  end subroutine do_error


  subroutine write_cell_info(i, j, k)
    implicit none

    integer, intent(in) :: i, j, k

    print *, ""
    print *, "CELL: ", i, j, k
    print *, "POSITION: ", xcent(i), ycent(j), zcent(k)
    print *, "DENS: ", U(DENS_VAR,i, j, k)
    print *, "PRES: ", U(PRES_VAR,i, j, k)
    print *, "GAMC: ", U(GAMC_VAR,i, j, k)

  end subroutine write_cell_info

End Subroutine hy_uhd_getFaceFlux
