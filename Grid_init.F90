!!****if* source/Grid/GridMain/paramesh/Grid_init
!!
!! NAME
!!  Grid_init
!!
!! SYNOPSIS
!!
!!  Grid_init(integer(IN) :: myPE,
!!            integer(IN) :: numProcs)
!!           
!!
!! DESCRIPTION
!!  Initialize Grid_data
!!
!! ARGUMENTS
!!
!!  myPE     -- current processor
!!  numProcs -- total number of processors
!!
!! PARAMETERS 
!!
!!  nblockx [INTEGER] 
!!     num initial blocks in x dir
!!  nblocky [INTEGER] 
!!     num initial blocks in y dir   
!!  nblockz [INTEGER] 
!!     num initial blocks in z dir   
!!  lrefine_max [INTEGER] 
!!      maximum AMR refinement level
!!  lrefine_min [INTEGER] 
!!      minimum AMR refinement level
!!  nrefs [INTEGER] 
!!      refine/derefine AMR grid every nrefs timesteps
!!
!!  refine_var_1 [INTEGER] 
!!     indicates first variable on which to refine
!!  refine_cutoff_1 [REAL] 
!!      threshold value to trigger refinement for refine_var_1
!!  derefine_cutoff_1 [REAL]
!!      threshold value to trigger derefinement for refine_var_1
!!  refine_filter_1 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_1
!!
!!  refine_var_2 [INTEGER] 
!!     indicates second variable on which to refine
!!  refine_cutoff_2 [REAL] 
!!      threshold value to trigger refinement for refine_var_2
!!  derefine_cutoff_2 [REAL]
!!      threshold value to trigger derefinement for refine_var_2
!!  refine_filter_2 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_2
!!
!!  refine_var_3 [INTEGER] 
!!     indicates third variable on which to refine (if needed)
!!  refine_cutoff_3 [REAL] 
!!      threshold value to trigger refinement for refine_var_3
!!  derefine_cutoff_3 [REAL]
!!      threshold value to trigger derefinement for refine_var_3
!!  refine_filter_3 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_3
!!
!!  refine_var_4 [INTEGER] 
!!     indicates fourth variable on which to refine (if needed)
!!  refine_cutoff_4 [REAL] 
!!      threshold value to trigger refinement for refine_var_4
!!  derefine_cutoff_4 [REAL]
!!      threshold value to trigger derefinement for refine_var_4
!!  refine_filter_4 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_4
!!
!!  flux_correct [BOOLEAN]
!!     turns on or off flux correction
!! small  [REAL]
!!   Generic small value that can be used as floor where needed
!! smlrho [REAL]  
!!   Cutoff value for density    
!! smallp [REAL]  
!!   Cutoff value for pressure
!! smalle [REAL]  
!!   Cutoff value for energy
!! smallt [REAL]  
!!   Cutoff value for temperature
!! smallu [REAL]  
!!   Cutoff value for velocity
!! smallx [REAL]  
!!   Cutoff value for abundances
!! eosMode[STRING]
!!   determines which variables to calculate from the ones
!!   defined. Possible values are "dens_ie", "dens_pres" and "dens_temp"
!! interpol_order [INTEGER]
!!   Order of interpolation, used in Paramesh2 "monotonic" interpolation
!!   for mesh prolongation
!! grid_monotone_hack [BOOLEAN]
!!   If .true., apply radical monotonicity constraints to interpolants,
!!   i.e., completely flatten them if they violate monotonicity.  Used
!!   in Paramesh2 "quadratic_cartesian" interpolation for mesh prolongation.
!! earlyBlockDistAdjustment [BOOLEAN]
!!   If .true., let Paramesh redistribute blocks
!!   across processors early, so that the block distribution chosen by
!!   Paramesh will be in effect when time evolution begins after restart.
!!   If earlyBlockDistAdjustment is .false., the block distribution enacted
!!   by the IO unit when it read a checkpoint file will normally still be
!!   in effect when time evolution begins after a restart.
!!   This flag is ignored if not restarting from a checkpoint.
!!
!!***

subroutine Grid_init(myPE, numProcs)

  use Grid_data
  use tree, ONLY : lrefine_min, lrefine_max, nfaces, nchild
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Simulation_interface, ONLY : Simulation_mapStrToInt, &
    Simulation_getVarnameType
  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: myPE,numProcs

  integer :: i, j, k, localNumBlocks, ii, numLeafBlks

  character(len=MAX_STRING_LENGTH) :: refine_var_string_1, &
       refine_var_string_2, refine_var_string_3, &
       refine_var_string_4, basename
  character(len=MAX_STRING_LENGTH) :: xl_bcString,xr_bcString
  character(len=MAX_STRING_LENGTH) :: yl_bcString,yr_bcString
  character(len=MAX_STRING_LENGTH) :: zl_bcString,zr_bcString
  character(len=MAX_STRING_LENGTH) :: eosModeString, grav_boundary_type
  real :: dx, dy, dz
  real, dimension(NDIM) :: rnb


!----------------------------------------------------------------------------------
! mesh geometry - moved here so Paramesh_init can use gr_geometry for some checking
!----------------------------------------------------------------------------------
  call RuntimeParameters_get("geometry",gr_str_geometry)
  call RuntimeParameters_mapStrToInt(gr_str_geometry, gr_geometry)

    !Store Grid copy of myPE and numProcs so accessor functions can use them.
  gr_myPE = myPE
  gr_numProcs = numProcs

  call Paramesh_init()

! The following renaming was done: "conserved_var" -> "convertToConsvdForMeshCalls". - KW
  call RuntimeParameters_get("convertToConsvdForMeshCalls", gr_convertToConsvdForMeshCalls)
  call RuntimeParameters_get("convertToConsvdInMeshInterp", gr_convertToConsvdInMeshInterp)
  if (gr_convertToConsvdInMeshInterp) then
#ifdef FLASH_GRID_PARAMESH2
        ! For PARAMESH3, if the new way of conversion to conserved
        ! form is requested, use the old way instead and generate
        ! appropriate warnings. - KW
     if (gr_convertToConsvdForMeshCalls) then
        if(myPE == MASTER_PE) print*,'WARNING : convertToConsvdInMeshInterp is ignored in PARAMESH2'
        call Logfile_stampMessage(MyPE,"WARNING : convertToConsvdInMeshInterp is ignored in PARAMESH2")
     else
        if(myPE == MASTER_PE) &
             print*,'WARNING : convertToConsvdInMeshInterp is not implemented in PARAMESH2'// &
               ', using convertToConsvdForMeshCalls logic instead.'
        call Logfile_stampMessage(MyPE, &
             'WARNING : convertToConsvdInMeshInterp is not implemented in PARAMESH2'// &
             ', using convertToConsvdForMeshCalls logic instead.')
        gr_convertToConsvdForMeshCalls = .TRUE.
     end if
#else
     if (gr_convertToConsvdForMeshCalls) then
        ! For PARAMESH3, if both ways of conversion to conserved form are requested,
        ! Let the new mechanism win and try to make sure the old one is not used. - KW
        if(myPE == MASTER_PE) &
             print*,'WARNING: convertToConsvdForMeshCalls ignored since convertToConsvdInMeshInterp is requested'
        call Logfile_stampMessage(MyPE, &
             "WARNING: convertToConsvdForMeshCalls ignored since convertToConsvdInMeshInterp is requested")
        gr_convertToConsvdForMeshCalls = .FALSE.
     end if
#endif
  end if

#ifndef FLASH_GRID_PARAMESH2
  call RuntimeParameters_get("enableMaskedGCFill", gr_enableMaskedGCFill)
#endif

  call RuntimeParameters_get("nrefs", gr_nrefs)
  call RuntimeParameters_get('lrefine_min', lrefine_min)
  call RuntimeParameters_get('lrefine_max', lrefine_max)
!  call RuntimeParameters_get('lrefine_limit', gr_lrefineLimit)

  call RuntimeParameters_get("smalle",gr_smalle)
  call RuntimeParameters_get("smlrho",gr_smallrho)
  call RuntimeParameters_get("smallx",gr_smallx) !
!  call RuntimeParameters_get("grid_monotone_hack", gr_monotone) ! for "quadratic_cartesian" interpolation
  call RuntimeParameters_get("interpol_order",gr_intpol) ! for "monotonic" interpolation
#ifdef GRID_WITH_MONOTONIC
  gr_intpolStencilWidth = 2     !Could possibly be less if gr_intpol < 2  - KW
#endif


  !get the boundary conditions stored as strings in the flash.par file
  call RuntimeParameters_get("xl_boundary_type", xl_bcString)
  call RuntimeParameters_get("xr_boundary_type", xr_bcString)
  call RuntimeParameters_get("yl_boundary_type", yl_bcString)
  call RuntimeParameters_get("yr_boundary_type", yr_bcString)
  call RuntimeParameters_get("zl_boundary_type", zl_bcString)
  call RuntimeParameters_get("zr_boundary_type", zr_bcString)

  !map the string boundary conditions to integer constants defined in constants.h
  call RuntimeParameters_mapStrToInt(xl_bcString,gr_domainBC(LOW,IAXIS))
  call RuntimeParameters_mapStrToInt(xr_bcString,gr_domainBC(HIGH,IAXIS))
  call RuntimeParameters_mapStrToInt(yl_bcString,gr_domainBC(LOW,JAXIS))
  call RuntimeParameters_mapStrToInt(yr_bcString,gr_domainBC(HIGH,JAXIS))
  call RuntimeParameters_mapStrToInt(zl_bcString,gr_domainBC(LOW,KAXIS))
  call RuntimeParameters_mapStrToInt(zr_bcString,gr_domainBC(HIGH,KAXIS))

  call RuntimeParameters_get("bndPriorityOne",gr_bndOrder(1))
  call RuntimeParameters_get("bndPriorityTwo",gr_bndOrder(2))
  call RuntimeParameters_get("bndPriorityThree",gr_bndOrder(3))

  !get the initial grid layout
  call RuntimeParameters_get("nblockx", gr_nBlockX) !number of initial blks in x dir
  call RuntimeParameters_get("nblocky", gr_nBlockY) !number of initial blks in y dir
  call RuntimeParameters_get("nblockz", gr_nblockZ) !number of initial blks in z dir  
  call RuntimeParameters_get("refine_on_particle_count",gr_refineOnParticleCount)

  call RuntimeParameters_get("min_particles_per_blk",gr_minParticlesPerBlk)
  call RuntimeParameters_get("max_particles_per_blk",gr_maxParticlesPerBlk)

#ifdef FLASH_GRID_PARAMESH2
  call RuntimeParameters_get( "msgbuf",i)
  gr_msgbuffer = (i==1)
#endif

!------------------------------------------------------------------------------
! mesh geometry       (gr_geometry was already set above)
!------------------------------------------------------------------------------
  !get the physical domain limits
  call RuntimeParameters_get('xmin', gr_imin)
  call RuntimeParameters_get('xmax', gr_imax)
  call RuntimeParameters_get('ymin', gr_jmin)
  call RuntimeParameters_get('ymax', gr_jmax)
  call RuntimeParameters_get('zmin', gr_kmin)
  call RuntimeParameters_get('zmax', gr_kmax)

  !Store computational domain limits in a convenient array.  Used later in Grid_getBlkBC.
  gr_globalDomain(LOW,IAXIS) = gr_imin
  gr_globalDomain(LOW,JAXIS) = gr_jmin
  gr_globalDomain(LOW,KAXIS) = gr_kmin
  gr_globalDomain(HIGH,IAXIS) = gr_imax
  gr_globalDomain(HIGH,JAXIS) = gr_jmax
  gr_globalDomain(HIGH,KAXIS) = gr_kmax


! Determine the geometries of the individual dimensions, and scale
! angle value parameters that are expressed in degrees to radians.
! This call must be made after gr_geometry, gr_domainBC, and gr_{j,k}{min,max}
! have been set based on the corresponding runtime parameters.
  call gr_initGeometry()


  call RuntimeParameters_get("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosMode)

  call RuntimeParameters_get("eosModeInit", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosModeInit)

  gr_eosModeNow = gr_eosModeInit ! may change after initialization is done

  call RuntimeParameters_get("earlyBlockDistAdjustment", gr_earlyBlockDistAdjustment)
  gr_justExchangedGC = .false.


 
  !get the refinement variables (ie. dens, pres, etc) as a string
  call RuntimeParameters_get( "refine_var_1", refine_var_string_1)
  call RuntimeParameters_get( "refine_var_2", refine_var_string_2)
  call RuntimeParameters_get( "refine_var_3", refine_var_string_3)
  call RuntimeParameters_get( "refine_var_4", refine_var_string_4)

  gr_refine_var = NONEXISTENT
  gr_numRefineVars=0

  !map the refinement vars to integer constants defined in Flash.h
  if(refine_var_string_1 /= "none") then
     call Simulation_mapStrToInt(refine_var_string_1, gr_refine_var(gr_numRefineVars+1),MAPBLOCK_UNK)
     if (gr_refine_var(gr_numRefineVars+1) .LE. 0) then
        if (MyPE==MASTER_PE) print*, 'WARNING: Unrecognized variable name in refine_var_1, treating it as "none"'
        call Logfile_stampMessage(MyPE,'WARNING: Unrecognized variable name in refine_var_1, treating it as "none"')
     else
        gr_numRefineVars=gr_numRefineVars+1
        call RuntimeParameters_get( "refine_cutoff_1", gr_refine_cutoff(gr_numRefineVars)  )
        call RuntimeParameters_get( "derefine_cutoff_1", gr_derefine_cutoff(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_filter_1",  gr_refine_filter(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_val_cutoff_1",  gr_refine_val_cutoff(gr_numRefineVars) )
     end if
  end if

  if(refine_var_string_2 /= "none") then
     call Simulation_mapStrToInt(refine_var_string_2, gr_refine_var(gr_numRefineVars+1),MAPBLOCK_UNK)
     if (gr_refine_var(gr_numRefineVars+1) .LE. 0) then
        if (MyPE==MASTER_PE) print*, 'WARNING: Unrecognized variable name in refine_var_2, treating it as "none"'
        call Logfile_stampMessage(MyPE,'WARNING: Unrecognized variable name in refine_var_2, treating it as "none"')
     else
        gr_numRefineVars=gr_numRefineVars+1
        call RuntimeParameters_get( "refine_cutoff_2", gr_refine_cutoff(gr_numRefineVars)  )
        call RuntimeParameters_get( "derefine_cutoff_2", gr_derefine_cutoff(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_filter_2",  gr_refine_filter(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_val_cutoff_2",  gr_refine_val_cutoff(gr_numRefineVars) )
     end if
  end if

  if(refine_var_string_3 /= "none") then
     call Simulation_mapStrToInt(refine_var_string_3, gr_refine_var(gr_numRefineVars+1),MAPBLOCK_UNK)
     if (gr_refine_var(gr_numRefineVars+1) .LE. 0) then
        if (MyPE==MASTER_PE) print*, 'WARNING: Unrecognized variable name in refine_var_3, treating it as "none"'
        call Logfile_stampMessage(MyPE,'WARNING: Unrecognized variable name in refine_var_3, treating it as "none"')
     else
        gr_numRefineVars=gr_numRefineVars+1
        call RuntimeParameters_get( "refine_cutoff_3", gr_refine_cutoff(gr_numRefineVars)  )
        call RuntimeParameters_get( "derefine_cutoff_3", gr_derefine_cutoff(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_filter_3",  gr_refine_filter(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_val_cutoff_3",  gr_refine_val_cutoff(gr_numRefineVars) )
     end if
  end if

  if(refine_var_string_4 /= "none") then
     call Simulation_mapStrToInt(refine_var_string_4, gr_refine_var(gr_numRefineVars+1),MAPBLOCK_UNK)
     if (gr_refine_var(gr_numRefineVars+1) .LE. 0) then
        if (MyPE==MASTER_PE) print*, 'WARNING: Unrecognized variable name in refine_var_4, treating it as "none"'
        call Logfile_stampMessage(MyPE,'WARNING: Unrecognized variable name in refine_var_4, treating it as "none"')
     else
        gr_numRefineVars=gr_numRefineVars+1
        call RuntimeParameters_get( "refine_cutoff_4", gr_refine_cutoff(gr_numRefineVars)  )
        call RuntimeParameters_get( "derefine_cutoff_4", gr_derefine_cutoff(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_filter_4",  gr_refine_filter(gr_numRefineVars) )
        call RuntimeParameters_get( "refine_val_cutoff_4",  gr_refine_val_cutoff(gr_numRefineVars) )
     end if
  end if
  call RuntimeParameters_get("lrefine_del", gr_lrefineDel)
  gr_maxRefine=lrefine_max
  if(gr_numRefineVars==0)then
     if(myPE == MASTER_PE) print*,'WARNING : Adaptive Grid did not find any refinement variables'
     call Logfile_stampMessage(MyPE,"WARNING : Adaptive Grid did not find any variable to refine")
  end if

#ifdef FLASH_PARTICLES
  call RuntimeParameters_get('useParticles',gr_useParticles)
  call RuntimeParameters_get('pt_maxPerProc',gr_maxParticlesPerProc)
#else
  gr_useParticles=.false.
#endif

  gr_allPeriodic = .true.
  do i = 1,NDIM
     if(gr_domainBC(LOW,i)/=PERIODIC)gr_allPeriodic=.false.
     if(gr_domainBC(HIGH,i)/=PERIODIC)gr_allPeriodic=.false.
  end do

  !Check if there are gravitational isolated boundary conditions
  !in order to determine which solvers to intialize.
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  gr_isolatedBoundaries = (grav_boundary_type=="isolated")

  gr_anyVarToConvert = .FALSE.
  do i = UNK_VARS_BEGIN,UNK_VARS_END
     gr_vars(i)=i
     call Simulation_getVarnameType(i, gr_vartypes(i))
     if (gr_vartypes(i) .eq. VARTYPE_PER_MASS) gr_anyVarToConvert = .TRUE.
  end do



  !! calculating deltas for each level of 
  !! refinement and putting them in the
  !! delta variable
  dx = gr_imax - gr_imin
  dy = gr_jmax - gr_jmin
  dz = gr_kmax - gr_kmin
  rnb = 0.0
  rnb(1) = dx/(1.0*NXB*gr_nBlockX)
#if NDIM > 1
  rnb(2) = dy/(1.0*NYB*gr_nBlockY)
#endif
#if NDIM > 2
  rnb(3) = dz/(1.0*NZB*gr_nBlockZ)
#endif  
  do i = 1,lrefine_max
     gr_delta(1:NDIM,i) = rnb
     rnb = rnb/2.0
  end do

      
  gr_minCellSize = ( gr_imax - gr_imin ) / ( gr_nblockX*NXB*2**(lrefine_max-1) ) 
  
  if (NDIM >= 2 .and. .not. gr_dirIsAngular(JAXIS) ) then
     gr_minCellSize = min( gr_minCellSize, ( gr_jmax - gr_jmin ) / &
          ( gr_nblockY*NYB*2**(lrefine_max-1) ) )
  end if
  if (NDIM == 3 .and. .not. gr_dirIsAngular(KAXIS) ) then
     gr_minCellSize = min( gr_minCellSize, (gr_kmax - gr_kmin ) / &
          ( gr_nblockZ*NZB*2**(lrefine_max-1) ) )
     
  end if

#ifdef FL_NON_PERMANENT_GUARDCELLS
  gr_blkPtrRefCount = 0 
  gr_blkPtrRefCount_fc = 0
  gr_lastBlkPtrGotten = 0 
  gr_lastBlkPtrGotten_fc = 0
#endif
  call gr_setDataStructInfo()
  call gr_ptInit()
  call gr_bcInit()

  !Initialize grid arrays used by IO
  allocate(gr_nToLeft(0:numProcs-1))
  allocate(gr_gid(nfaces+nchild+1, MAXBLOCKS))

end subroutine Grid_init
