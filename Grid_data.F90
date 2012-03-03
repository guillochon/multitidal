!!****if* source/Grid/GridMain/paramesh/Grid_data
!!
!! NAME
!!  Grid_data
!!
!! SYNOPSIS
!!
!!  use Grid_data
!!
!! DESCRIPTION 
!!  
!!  This includes the global integer identifier for a block, the grid geometry information
!!  
!!  
!! 
!! CREATE AD:04/12/04
!!
!! 
!!   Defining data structures for storing paramesh related infomation.  
!!   including function for updating the grid information
!!
!! MODIFIED AD:05/19/04
!!   
!!***

!!REORDER(5):scratch, scratch_ctr, scratch_facevar[xyz], gr_[xyz]flx

Module Grid_data

  implicit none

#include "constants.h"
#include "Flash.h"

!! Define the block information

!  Type define

  integer, parameter :: somethingBig = 20


  real,save,dimension(3,somethingBig) :: gr_delta
  integer, save :: gr_iloGc = GRID_ILO_GC
  integer, save :: gr_ihiGc = GRID_IHI_GC
  integer, save :: gr_jloGc = GRID_JLO_GC
  integer, save :: gr_jhiGc = GRID_JHI_GC
  integer, save :: gr_kloGc = GRID_KLO_GC
  integer, save :: gr_khiGc = GRID_KHI_GC
  integer, save :: gr_iguard = NGUARD
  integer, save :: gr_jguard = NGUARD 
  integer, save :: gr_kguard = NGUARD

  integer, save :: gr_ilo = GRID_ILO
  integer, save :: gr_ihi = GRID_IHI
  integer, save :: gr_jlo = GRID_JLO
  integer, save :: gr_jhi = GRID_JHI
  integer, save :: gr_klo = GRID_KLO
  integer, save :: gr_khi = GRID_KHI
  integer,save,dimension(MDIM)::gr_bndOrder
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vars
  integer,save,dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_vartypes
  logical, save :: gr_anyVarToConvert
  logical, save :: gr_justExchangedGC, gr_allPeriodic, gr_isolatedBoundaries
  logical, save :: gr_useParticles,gr_refineOnParticleCount,gr_refineOnPdens
  integer, save :: gr_minParticlesPerBlk, gr_maxParticlesPerBlk
  integer, save :: gr_globalComm, gr_globalMe, gr_globalNumProcs
  integer, save :: gr_meshComm, gr_meshMe, gr_meshNumProcs
  integer, save :: gr_meshAcrossComm, gr_meshAcrossMe, gr_meshAcrossNumProcs

  logical, save :: gr_useEnergyDeposition

  type gridBlock
     !!cornerID is integer coordinates of the lower left cornor
     !! (ie the smallest point) of a block
     integer,dimension(MDIM) :: cornerID
     !! atmost 2 neighbors, 2faces along
     !! each dimension, hence.
     real,dimension(3,GRID_IHI_GC) :: firstAxisCoords
     real,dimension(3,GRID_JHI_GC) :: secondAxisCoords
     real,dimension(3,GRID_KHI_GC) :: thirdAxisCoords
     integer :: blockType
  end type gridBlock

  type(gridBlock),save,dimension(MAXBLOCKS),target :: gr_oneBlock
  
#if NSCRATCH_GRID_VARS > 0
  real,target,dimension(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
                        GRID_IHI_GC+1,&
                        GRID_JHI_GC+1,&
                        GRID_KHI_GC+1,&
                        MAXBLOCKS) :: scratch
#else
  real,target,dimension(1,1,1,1,1) :: scratch
#endif

#if NSCRATCH_CENTER_VARS > 0
  real,target,dimension(SCRATCH_CENTER_VARS_BEGIN:SCRATCH_CENTER_VARS_END,&
                        GRID_IHI_GC,&
                        GRID_JHI_GC,&
                        GRID_KHI_GC,&
                        MAXBLOCKS) :: scratch_ctr
#else
  real,target,dimension(1,1,1,1,1):: scratch_ctr
#endif

#if NSCRATCH_FACEX_VARS > 0
  real,target,dimension(SCRATCH_FACEX_VARS_BEGIN:SCRATCH_FACEX_VARS_END,&
                        GRID_IHI_GC+1,&
                        GRID_JHI_GC,  &
                        GRID_KHI_GC,  &
                        MAXBLOCKS) :: scratch_facevarx
#else
  real, target,dimension(1,1,1,1,1):: scratch_facevarx
#endif

#if NSCRATCH_FACEY_VARS > 0
  real,target,dimension(SCRATCH_FACEY_VARS_BEGIN:SCRATCH_FACEY_VARS_END,&
                        GRID_IHI_GC,  &
                        GRID_JHI_GC+1,&
                        GRID_KHI_GC,  &
                        MAXBLOCKS) :: scratch_facevary
#else
  real, target,dimension(1,1,1,1,1):: scratch_facevary
#endif

#if NSCRATCH_FACEZ_VARS > 0
  real,target,dimension(SCRATCH_FACEZ_VARS_BEGIN:SCRATCH_FACEZ_VARS_END,&
                        GRID_IHI_GC,  &
                        GRID_JHI_GC,  &
                        GRID_KHI_GC+1,&
                        MAXBLOCKS) :: scratch_facevarz
#else
  real, target,dimension(1,1,1,1,1):: scratch_facevarz
#endif


  integer, save :: gr_eosMode
  integer, save :: gr_eosModeInit, gr_eosModeNow
  integer, save :: gr_oneRefLev=1 !! To be used with the multigrid
  integer ,save :: gr_nrefs
  logical ,save :: gr_convertToConsvdForMeshCalls
  logical ,save :: gr_convertToConsvdInMeshInterp
  logical ,save :: gr_earlyBlockDistAdjustment
  logical, save :: gr_monotone
  integer, save :: gr_intpol
  real, save :: gr_smallx

  integer, save :: gr_numRefineVars, gr_numRefineVarsMax
  integer,allocatable,dimension(:) ,save :: gr_refineVars
  real ,save :: gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
  real, save, dimension(LOW:HIGH,MDIM) :: gr_globalDomain
  integer,save,dimension(2,MDIM) :: gr_domainBC, gr_blkBC
  integer ,save :: gr_geometry
  integer, save, dimension(MDIM) :: gr_dirGeom
  logical, save, dimension(MDIM) :: gr_dirIsAngular
  logical, save :: gr_geometryOverride
  character(len=MAX_STRING_LENGTH) :: gr_str_geometry
  integer ,save :: gr_nblockX, gr_nblockY, gr_nblockZ
  integer ,save, dimension(MAXREFVARS) :: gr_refine_var    
  real,save,dimension(MAXCELLS,UNK_VARS_BEGIN:UNK_VARS_END) :: gr_oneRow
  real,dimension(MAXREFVARS), save::gr_refine_cutoff,&
       gr_derefine_cutoff,gr_refine_filter,gr_refine_val_cutoff
  integer,dimension(MAXREFVARS), save :: gr_refine_level
  real, save :: gr_smalle,gr_smallrho
  integer, save :: gr_blkCount, gr_matchCount
  integer ,dimension(MAXBLOCKS), save :: gr_blkList,gr_matchList
  integer, allocatable:: gr_mortons(:)
  real, save :: gr_minCellSize
  real, save, dimension(MDIM) :: gr_minCellSizes

  !below values needed to make data structures for IO output
  integer,save :: gr_globalOffset !stores beginning blk offset for a proc
  integer,save,allocatable,target,dimension(:,:) :: gr_gid  !holds neigh, child, parent info for checkpoint files
  integer, save :: gr_globalNumBlocks !
  integer, save, allocatable :: gr_nToLeft(:) !holds 
  real,save, dimension(NFLUXES,2,NYB,NZB,MAXBLOCKS) :: gr_xflx
  real,save, dimension(NFLUXES,NXB,2,NZB,MAXBLOCKS) :: gr_yflx
  real,save, dimension(NFLUXES,NXB,NYB,2,MAXBLOCKS) :: gr_zflx

#ifdef FLASH_GRID_PARAMESH2
  logical, save :: gr_msgbuffer 
  logical, parameter :: gr_enableMaskedGCFill = .FALSE.
#else
  logical, save :: gr_enableMaskedGCFill
#endif

#ifdef GRID_WITH_MONOTONIC
  integer, save :: gr_intpolStencilWidth
#else
  ! The following is for Paramesh3f with native interpolation
  integer, parameter :: gr_intpolStencilWidth = 1
#endif

#ifdef FLASH_PARTICLES


  !gr_blkParticleInfo is a multidimensional array holding various info
  !needed for updating the particle refinement
  !gr_blkParticleInfo(1,:) holds cornerID for 1st dim
  !gr_blkParticleInfo(2,:) holds cornerID for 2nd dim
  !gr_blkParticleInfo(3,:) holds cornerID for 3rd dim
  !gr_blkParticleInfo(4,:) keeps track of if an old block has been processed or not


  integer, save :: gr_blkParticleInfo(4,MAXBLOCKS)
  integer, save :: gr_maxParticlesPerProc
#endif

  integer, save :: gr_numDataStruct
  integer, save, dimension(NDATATYPES) :: gr_gridDataStruct,gr_gridDataStructSize

#ifdef FL_NON_PERMANENT_GUARDCELLS
  integer,save :: gr_blkPtrRefCount, gr_blkPtrRefCount_fc
  integer,save :: gr_lastBlkPtrGotten, gr_lastBlkPtrGotten_fc
  logical,save,dimension(NUNK_VARS) :: gr_ccMask
#if(NFACE_VARS>0)
  logical,save,dimension(3,NFACE_VARS) :: gr_fcMask
#else
  logical :: gr_fcMask
#endif

#endif

  integer,save :: gr_lrefineDel, gr_maxRefine
  logical,save :: gr_enforceMaxRefinement
  logical,save :: gr_lrefineMaxRedDoByLogR, gr_lrefineMaxRedDoByTime
  real,save :: gr_lrefineMaxRedRadiusSq
  real,save :: gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK
  real,save :: gr_lrefineMaxRedTimeScale, gr_lrefineMaxRedTRef, gr_lrefineMaxRedLogBase
  integer,save :: gr_restrictAllMethod

! Stuff for customRegion 
  logical, save :: gr_useQuietStart, gr_usePiston, gr_useHole
  real, save, dimension(LOW:HIGH,MDIM) :: gr_quietStartBnd
  real, save, dimension(LOW:HIGH,MDIM) :: gr_pistonBnd
  real, save :: gr_holeRad, gr_holeTime, gr_holeVel
  integer, save :: gr_holeBnd

  real, save :: gr_quietStartDens, gr_quietStartTemp
  real, save :: gr_pistonDens, gr_pistonVelx, gr_pistonVely, gr_pistonVelz
  real, save, allocatable, dimension(:,:,:,:) :: gr_customRegionData

end Module Grid_data





