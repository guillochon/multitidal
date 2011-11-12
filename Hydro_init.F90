!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro_init
!!
!! NAME
!!
!!  Hydro_init
!!
!!
!! SYNOPSIS
!!
!!  Hydro_init(integer(IN) :: myPE)
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  myPE -- local processor number
!!
!! PARAMETERS
!!
!!   These are the runtime parameters used in the split PPM Hydro 
!!   implementation.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have overwritten these values with the flash.par values
!!   for your specific run.  
!!
!!    cfl [REAL]
!!    geometry [STRING]
!!        Grid geometry
!!    eosMode [STRING]
!!        the default Eos mode, usually MODE_DENS_EI, 
!!        where density and energy are provided to 
!!        calculate pressure and temperature
!!    irenorm
!!    flux_correct [BOOLEAN]
!!    hybrid_riemann [BOOLEAN]
!!    smlrho [REAL]
!!        Cutoff value for density
!!    smallp [REAL]
!!        Cutoff value for pressure
!!    eintSwitch [REAL]  Defined in Eos Unit
!!    cvisc [REAL]
!!        Artificial viscosity constant
!!    dp_sh_md
!!    epsiln [REAL]
!!    useShockBurn [BOOLEAN]
!!    nriem [INTEGER]
!!    omg1
!!        PPM dissipation parameter omega1
!!    omg2
!!        PPM dissipation parameter omega2
!!    rieman_tol [REAL]
!!    small [REAL]
!!       Generic small value that can be used as floor where needed
!!    smallu [REAL]
!!       Cutoff value for velocity
!!    smallx [REAL]
!!       Cutoff value for abundances
!!    vgrid
!!    ppm_modifystates
!!    leveque
!!    igodu [INTEGER]
!!       Enable Guodunov method instead of PPM if set to 1
!!    iplm [INTEGER]
!!       Enable piecewise linear method instead of PPM if set to 1
!!    use_steepening [BOOLEAN]
!!    use_cma_flattening [BOOLEAN]
!!    use_cma_steepening [BOOLEAN]
!!    use_cma_advection [BOOLEAN]
!!
!!***

subroutine Hydro_init(myPE)

  !!These are all the runtime parameters.  First the logicals, then the
  !! integers, then the reals    

  use Hydro_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage, &
    Logfile_stampVarMask
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Grid_interface, ONLY: Grid_getNumProcs, Grid_setFluxHandling

  implicit none

  integer, intent(IN) :: myPE

  logical :: useSpecialFluxVar
  integer :: istat

#include "constants.h"
#include "Flash.h"  

  character(len=MAX_STRING_LENGTH) :: str_geometry,eosModeString
  
  ! Everybody should know these
  hy_myPE = myPE
  call Grid_getNumProcs(hy_numProcs)
  
  !!hydro_timestep
  call RuntimeParameters_get ("cfl", hy_cfl)

  !!**Hydro_sweep RuntimeParameters

  call RuntimeParameters_get ("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, hy_geometry)

  call RuntimeParameters_get ("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, hy_eosMode)
  if(hy_eosMode/=MODE_DENS_EI)&
       call Driver_abortFlash("Hydro : Wrong Eos mode for PPM")
  hy_useGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity", hy_useGravity)
#endif
  
  call RuntimeParameters_get("irenorm", hy_irenorm)
  call RuntimeParameters_get("flux_correct",   hy_fluxCorrect)
  call RuntimeParameters_get("hybrid_riemann", hy_hybridRiemann)

  !!**Hydro_updateSolution
  call RuntimeParameters_get("smlrho",hy_smlrho)
  call RuntimeParameters_get("smallp",hy_smallp)  
  call RuntimeParameters_get("eintSwitch",hy_eintSwitch)  
 
  !!**Hydro_1d needs
  call RuntimeParameters_get("cvisc", hy_cvisc)
  call RuntimeParameters_get("dp_sh_md", hy_dp_sh_md )
  call RuntimeParameters_get("epsiln", hy_epsiln)

! Define this ifdef variable in your unit's Config file with PPDEFINE
#ifdef FLASH_SOURCEBURN
  call RuntimeParameters_get("useShockBurn", hy_useShockBurn)
#else
  hy_useShockBurn = .false.
#endif
  call RuntimeParameters_get("nriem", hy_nriem)
  call RuntimeParameters_get("omg1", hy_omg1)
  call RuntimeParameters_get("omg2", hy_omg2)
  call RuntimeParameters_get("rieman_tol", hy_riemanTol)
  call RuntimeParameters_get("small", hy_small)
  !!smallp from update_soln
  call RuntimeParameters_get("smallu",hy_smallu)
  call RuntimeParameters_get("smallx",hy_smallx)
  !!smlrho from update_soln
  call RuntimeParameters_get("vgrid", hy_vgrid)

  !!**PPM inputs
  call RuntimeParameters_get("ppm_modifystates", hy_ppmModifystates)
  call RuntimeParameters_get("leveque",hy_leveque)
  call RuntimeParameters_get("igodu", hy_igodu)
  call RuntimeParameters_get("iplm", hy_iplm)
  call RuntimeParameters_get("use_steepening", hy_useSteepening)
  call RuntimeParameters_get("use_cma_flattening", hy_useCmaFlattening)
  call RuntimeParameters_get("use_cma_steepening", hy_useCmaSteepening)
  call RuntimeParameters_get("use_cma_advection", hy_useCmaAdvection)

  call RuntimeParameters_get("hy_fluxRepresentation", hy_fluxRepresentation)
  if (trim(hy_fluxRepresentation) == "auto") then
#ifdef FLASH_GRID_PARAMESH2
     hy_fluxRepresentation = "hybrid" !Paramesh2 always assumes this - KW
#else
     if (hy_geometry == CARTESIAN) then
        hy_fluxRepresentation = "hybrid"
     else
        hy_fluxRepresentation = "fluxes"
     end if
#endif
  end if

  if (hy_fluxRepresentation == "hybrid") then
     hy_useCellAreasForFluxes = .FALSE.
     call Grid_setFluxHandling('consv_flux_densities',status=istat)
  else if (hy_fluxRepresentation == "fluxes") then
     hy_useCellAreasForFluxes = .TRUE.
     call Grid_setFluxHandling('consv_fluxes',status=istat)
  else
     call Driver_abortFlash('Hydro_init: Runtime Parameter hy_fluxRepresentation must be either '//&
          '"hybrid" or "fluxes".')
  end if
#if NDIM > 1
  if (istat .NE. 0) then
     if (hy_fluxCorrect) then
        if (myPE .EQ. MASTER_PE) print*,'WARNING from Hydro_init: hy_fluxRepresentation '//&
             'was requested as "',hy_fluxRepresentation,'",',' but the Grid unit does not support this,'//&
             ' using the handling supported by Grid.'
        call Logfile_stampMessage(myPE,'WARNING from Hydro_init: hy_fluxRepresentation '//&
             'was requested as "'//hy_fluxRepresentation//'", but the Grid unit does not support this,'//&
             ' using the handling supported by Grid.')
        hy_useCellAreasForFluxes = .NOT. hy_useCellAreasForFluxes
     end if
  end if
#endif

  useSpecialFluxVar = hy_useCellAreasForFluxes

  if (useSpecialFluxVar) then
     hy_specialFluxVar = P_FLUX
  end if
  if (myPE .EQ. MASTER_PE) then
     print*,'Info: Hydro_init has set hy_specialFluxVar to ',hy_specialFluxVar
  end if
  

!! Determine some unit-wide variables that are directly derived from
!! Runtime parameters
  if (hy_cvisc == 0.0) then
     hy_transverseStencilWidth = 0
  else
     hy_transverseStencilWidth = 1
  end if

!! Determine the geometries of the individual dimensions

  if (hy_geometry == CARTESIAN)then
     hy_dirGeom(IAXIS) = XYZ
     hy_dirGeom(JAXIS) = XYZ
     hy_dirGeom(KAXIS) = XYZ
  elseif(hy_geometry == POLAR)then
     hy_dirGeom(IAXIS) = RAD_CYL
     hy_dirGeom(JAXIS) = PHI_CYL
     hy_dirGeom(KAXIS) = XYZ
  elseif(hy_geometry == CYLINDRICAL) then
     hy_dirGeom(IAXIS) = RAD_CYL
     hy_dirGeom(JAXIS) = XYZ
     hy_dirGeom(KAXIS) = PHI_CYL
  elseif(hy_geometry == SPHERICAL) then
     hy_dirGeom(IAXIS) = RAD_SPH
     hy_dirGeom(JAXIS) = THETA
     hy_dirGeom(KAXIS) = PHI_SPH
  else
     call Driver_abortFlash("unsupported geometry ")
  end if
  
  if (hy_vgrid /= 0.e0) then
     hy_movingGrid = .true.
  else
     hy_movingGrid = .false.
  endif

!! Now initialize the GC Mask

  hy_gcMask = .FALSE.

#ifdef FL_NON_PERMANENT_GUARDCELLS
  hy_gcMask(PRES_VAR) = .TRUE.
#else
  if (hy_hybridRiemann .AND. (hy_cvisc .ne. 0.0)) hy_gcMask(PRES_VAR) = .TRUE.
#endif

  hy_gcMask(DENS_VAR) = .TRUE.
  hy_gcMask(ENER_VAR) = .TRUE.
#ifdef EINT_VAR
  hy_gcMask(EINT_VAR) = .TRUE.
#endif
  hy_gcMask(TEMP_VAR) = .TRUE.     !for now - only used fo initial guess by Helmholtz Eos - KW

#ifdef GPOT_VAR
  hy_gcMask(GPOT_VAR) = .TRUE.
#endif
#ifdef GPOL_VAR
  hy_gcMask(GPOL_VAR) = .TRUE.
#endif
#ifdef GPO2_VAR
  hy_gcMask(GPO2_VAR) = .TRUE.
#endif

  hy_gcMask(VELX_VAR) = .TRUE.
#if NDIM >= 2
#ifdef VELY_VAR
  hy_gcMask(VELY_VAR) = .TRUE.
#endif
#endif
#if NDIM == 3
#ifdef VELZ_VAR
  hy_gcMask(VELZ_VAR) = .TRUE.
#endif
#endif
#if SPECIES_BEGIN <= UNK_VARS_END
  hy_gcMask(SPECIES_BEGIN:UNK_VARS_END) = .TRUE.
#endif

  hy_gcMaskSize=NUNK_VARS

  call Logfile_stampVarMask(hy_gcMask, .FALSE., '[Hydro_init]', 'gcNeed')

end subroutine Hydro_init
