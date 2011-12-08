!!****if* source/physics/Hydro/HydroMain/split/PPM/hy_ppm_sweep
!!
!! NAME
!!
!!  hy_ppm_sweep
!!
!! SYNOPSIS
!!
!!  call hy_ppm_sweep(integer (IN):: blockCount, 
!!                    integer (IN):: blockList(blockCount) ,
!!                    real    (IN):: timeEndAdv, 
!!                    real    (IN):: dt, 
!!                    real    (IN):: dtOld,            
!!                    integer (IN):: sweepDir )
!!
!! DESCRIPTION
!!
!!  Performs a hydro update of the state variables, species
!!  abundances, and mass scalars and then an update of dependant
!!  thermodynamic variables by doing a one dimensional sweep and
!!  update over a set of blocks, followed by a call to the eos on all
!!  cells of the blocks. sweepDir determines the direction of the one
!!  dimensional sweep.  blockList is an integer array of size
!!  blockCount that contains the blocks over which to sweep.  dt
!!  gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.  dtOld is the previous
!!  timestep taken.
!! 
!!  This routine performs a guardcell fill and for each block: 
!!   - applies an eos to the guard cells; 
!!   - computes fluxes (with the hydro solver with a call to
!!   hy_ppm_block)
!!   - if we're not doing flux correction (as controlled by the flux_correct
!!   runtime parameter), then update all the cell values from the fluxes (with a
!!   call to hy_ppm_updateSoln), otherwise, update just cells not on the
!!   boundaries, and save fluxes for cells on the boundary;
!!   - and finally, apply an eos to the block. 
!! 
!!  After the main block loop, if doing flux correction, have
!!  the Grid correct boundary fluxes for all blocks where approriate,
!!  and do another loop over blocks, updating the cell values for
!!  cells on the block boundaries using the corrected fluxes, and
!!  apply an eos on the block.
!! 
!! ARGUMENTS
!!
!!  blockCount -  number of blocks to advance
!!  blockList -   local block numbers of blocks to advance
!!  timeEndAdv -  end time
!!  dt -          timestep
!!  dtOld -       old timestep
!!  sweepDir -    direction of hydro sweep
!!
!!***


! solnData depends on the ordering on unk
!!REORDER(4): solnData, tempFlx


#ifdef DEBUG_ALL 
#define DEBUG_HYDRO
#endif
#define DEBUG_GRID_GCMASK

subroutine hy_ppm_sweep (  blockCount, blockList, &
                         timeEndAdv, dt, dtOld,  &
                         sweepDir )

  use Hydro_data, ONLY : hy_updateHydroFluxes, hy_useDiffuse,    &
                         hy_fluxCorrect,  hy_gravMass,    &
                         hy_irenorm, hy_hybridRiemann, &
                         hy_dirGeom, hy_eosMode, hy_transverseStencilWidth, &
                         hy_cvisc,hy_gcMask,hy_gcMaskSize, &
                         hy_specialFluxVars, &
                         hy_alwaysCallDetectShock, &
                         hy_meshMe, hy_meshNumProcs
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getDeltas, &
    Grid_getBlkBC, Grid_getBlkIndexLimits, Grid_getCellCoords, &
    Grid_getBlkPtr, Grid_putBlkData, Grid_releaseBlkPtr, Grid_putFluxData, &
    Grid_conserveFluxes, Grid_getFluxData, Grid_renormAbundance, &
    Grid_limitAbundance, &
    Grid_markCustomRegion, Grid_updateCustomRegion, Grid_fillCustomRegion
  use Eos_interface, ONLY : Eos_wrapped
  use Diffuse_interface, ONLY : Diffuse_therm, Diffuse_visc
  use Hydro_interface, ONLY : Hydro_detectShock
  use hy_ppm_interface, ONLY : hy_ppm_block, hy_ppm_updateSoln,&
                           hy_ppm_getTemporaryData, hy_ppm_putTemporaryData
  use IO_interface, ONLY: IO_writeCheckpoint
  implicit none

#include "constants.h"
#include "Flash.h"
#include "PPM.h"
#include "Eos.h"


  !! ----------------------
  !! ---- INPUT ARGUMENTS
  integer, INTENT(IN) :: blockCount
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld
  integer, INTENT(IN) :: sweepDir
  integer, INTENT(IN), dimension(blockCount) :: blockList


  real, pointer :: solnData(:,:,:,:)   
  integer isize,jsize,ksize,isizeGC,jsizeGC,ksizeGC
  integer,dimension(MDIM) :: size,startingPos
  !! ---- SOLUTION UPDATE/FLUX data structures
  !! These are set in hy_ppm_block from the 1d slices of
  !! calculated fluxes.   They are used to update the solution, and
  !! given to Grid package to for flux matching routines. 

#ifdef FIXEDBLOCKSIZE
  real, DIMENSION(NFLUXES, GRID_ILO_GC:GRID_IHI_GC, &
                             GRID_JLO_GC:GRID_JHI_GC, &
                             GRID_KLO_GC:GRID_KHI_GC) :: tempFlx

  !! These are used for a proper solution update.  We do not adjust
  !! the actual solnData in the kernel, but, return these values to update.
  real, DIMENSION(GRID_ILO_GC:GRID_IHI_GC,                   &
                  GRID_JLO_GC:GRID_JHI_GC,                   &
                  GRID_KLO_GC:GRID_KHI_GC) ::  tempArea,     &
                                               tempGrav1d_o, &
                                               tempGrav1d,   &
                                               tempDtDx,     &
                                               tempFict,     &
                                               tempAreaLeft, &
                                               tempDens,     & !Added by JFG
                                               tempDens2       !Added by JFG
  real, DIMENSION(GRID_ILO_GC:GRID_IHI_GC,                   &
                  GRID_JLO_GC:GRID_JHI_GC,                   &
                  GRID_KLO_GC:GRID_KHI_GC) ::  shock

  real, DIMENSION(GRID_ILO_GC:GRID_IHI_GC,                   &
                  GRID_JLO_GC:GRID_JHI_GC,                   &
                  GRID_KLO_GC:GRID_KHI_GC) ::  regionType

  logical, DIMENSION(GRID_ILO_GC:GRID_IHI_GC,                   &
                  GRID_JLO_GC:GRID_JHI_GC,                   &
                  GRID_KLO_GC:GRID_KHI_GC) ::  isGc
  real,DIMENSION(MAXCELLS) :: primaryCoord ,  &
                              primaryLeftCoord , &
                              primaryRghtCoord , &
                              primaryDx        , &
                              secondCoord      , &
                              thirdCoord       , &
                              radialCoord     , &
                              ugrid
 
  integer, parameter :: numCells = MAXCELLS
                                
#else
  real, allocatable, DIMENSION(:,:,:,:) :: tempFlx

  !! These are used for a proper solution update.  We do not adjust
  !! the actual solnData in the kernel, but, return these values to update.
  real, allocatable, DIMENSION(:,:,:) ::  tempArea,     &
                                          tempGrav1d_o, &
                                          tempGrav1d,   &
                                          tempDtDx,     &
                                          tempFict,     &
                                          tempAreaLeft
  real, allocatable, DIMENSION(:,:,:) ::  shock, regionType

  real,allocatable, DIMENSION(:) :: primaryCoord ,  &
                                    primaryLeftCoord , &
                                    primaryRghtCoord , &
                                    primaryDx        , &
                                    secondCoord      , &
                                    thirdCoord       , &
                                    radialCoord     , &
                                    ugrid
  logical,allocatable,dimension(:,:,:) :: isGc
  integer :: numCells
  

#endif

  integer :: i,j,k,blk
  integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC,eosRange,bcs
  real, dimension(MDIM) :: del

  integer :: numGuard
  integer :: updateMode = UPDATE_ALL 
  logical :: gcell = .true.
  logical :: doFluxCorrect
  integer :: igeom
  integer :: idir, level = 0
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged(MDIM) =.FALSE.
#endif


!!----------------------------
!!------End Data Declarations
!!----------------------------

!  call IO_writeCheckpoint()

  doFluxCorrect = (hy_fluxCorrect .AND. (hy_updateHydroFluxes .OR. hy_useDiffuse))
  if (doFluxCorrect) then
     updateMode = UPDATE_INTERIOR
  else
     updateMode = UPDATE_ALL
  end if

  if(sweepDir==SWEEP_X)idir=IAXIS
  if(sweepDir==SWEEP_Y)idir=JAXIS
  if(sweepDir==SWEEP_Z)idir=KAXIS

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged(idir)) then
     call Logfile_stampVarMask(hy_gcMask, .FALSE., '[hy_sweep]', 'gcNeed')
     gcMaskLogged(idir) = .TRUE.
  end if
#endif
  call Grid_fillGuardCells( CENTER, ALLDIR, eosMode=hy_eosMode,doEos=.FALSE.,&
       maskSize=hy_gcMaskSize, mask=hy_gcMask, makeMaskConsistent=.true.)


!DEV: We may want to restore this at some point.  DO NOT DELETE!
!!$  call Grid_fillGuardCells( CENTER, idir,minLayers=hy_transverseStencilWidth)

!!------------------------------------------------------------------
!!LOOP OVER THE BLOCKS AND RUN HYDRO
!!------------------------------------------------------------------

  call Timers_start("hy_ppm_sweep")

  do blk=1,blockCount

     call Grid_getDeltas(blockList(blk),del)
     call Grid_getBlkBC(blockList(blk),bcs)
     if (bcs(LOW,idir) == HYDROSTATIC_NVREFL .OR. bcs(LOW,idir) == HYDROSTATIC_F2_NVREFL) &
          bcs(LOW,idir) = REFLECTING
     if (bcs(HIGH,idir) == HYDROSTATIC_NVREFL .OR. bcs(HIGH,idir) == HYDROSTATIC_F2_NVREFL) &
          bcs(HIGH,idir) = REFLECTING
     call Grid_getBlkIndexLimits(blockList(blk),blkLimits,blkLimitsGC)

     


     isizeGC=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsizeGC=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksizeGC=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     

     isize=blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     jsize=blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     ksize=blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE

     numCells = max(isizeGC,jsizeGC)
     numcells = max(numCells,ksizeGC)

     allocate(tempFlx(NFLUXES,isizeGC,jsizeGC,ksizeGC))

     allocate(tempArea(isizeGC,jsizeGC,ksizeGC))
     allocate(tempGrav1d_o(isizeGC,jsizeGC,ksizeGC))
     allocate(tempGrav1d(isizeGC,jsizeGC,ksizeGC))
     allocate(tempDtDx(isizeGC,jsizeGC,ksizeGC))
     allocate(tempFict(isizeGC,jsizeGC,ksizeGC))
     allocate(tempAreaLeft(isizeGC,jsizeGC,ksizeGC))

     allocate(shock(isizeGC,jsizeGC,ksizeGC))
     allocate(regionType(isizeGC,jsizeGC,ksizeGC))

     allocate(primaryCoord(numCells))
     allocate(primaryLeftCoord(numCells))
     allocate(primaryRghtCoord(numCells))
     allocate(primaryDx(numCells))
     allocate(secondCoord(numCells))
     allocate(thirdCoord(numCells))
     allocate(radialCoord(numCells))
     allocate(ugrid(numCells))
     allocate(isGc(isizeGC,jsizeGC,ksizeGC))

#endif
     !! This logical array is the size of a single block
     !! Its values are true on guardcells, and false on the interior
     isGc=.true.
     isGc(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))=.false.
     !!Note that this assumes that SWEEP_X == 1, SWEEP_Y == 2, SWEEP_Z == 3.
     primaryDx = del(sweepDir)  
     primaryCoord = 0.0
     primaryLeftCoord = 0.0
     primaryRghtCoord = 0.0
     secondCoord = 0.0
     thirdCoord  = 0.0
     ugrid(:) = 0.0 !!FIX

     call Grid_getCellCoords&
          (IAXIS,blockList(blk),CENTER,gcell,radialCoord,isizeGC)

     if (sweepDir == SWEEP_X) then
        igeom = hy_dirGeom(IAXIS)
        primaryCoord = radialCoord
        call Grid_getCellCoords(IAXIS,blockList(blk),&
             LEFT_EDGE,gcell,primaryLeftCoord,isizeGC)
        call Grid_getCellCoords(IAXIS,blockList(blk),&
             RIGHT_EDGE,gcell,primaryRghtCoord,isizeGC)
        call Grid_getCellCoords(JAXIS,blockList(blk),&
             CENTER,gcell,secondCoord,jsizeGC)
        call Grid_getCellCoords(KAXIS,blockList(blk),&
             CENTER,gcell,thirdCoord,ksizeGC)
        numGuard = blkLimits(LOW,IAXIS)-1
     else if (sweepDir == SWEEP_Y) then
        igeom = hy_dirGeom(JAXIS)
        secondCoord = radialCoord
        call Grid_getCellCoords(JAXIS,blockList(blk),&
             CENTER,gcell,primaryCoord,jsizeGC)
        call Grid_getCellCoords(JAXIS,blockList(blk),&
             LEFT_EDGE,gcell,primaryLeftCoord,jsizeGC)
        call Grid_getCellCoords(JAXIS,blockList(blk),&
             RIGHT_EDGE,gcell,primaryRghtCoord,jsizeGC)
        call Grid_getCellCoords(KAXIS,blockList(blk),&
             CENTER,gcell,thirdCoord,ksizeGC)
        numGuard = blkLimits(LOW,JAXIS)-1
     else 
        igeom = hy_dirGeom(KAXIS)
        secondCoord = radialCoord
        call Grid_getCellCoords(KAXIS,blockList(blk),&
             CENTER,gcell,primaryCoord,ksizeGC)
        call Grid_getCellCoords(KAXIS,blockList(blk),&
             LEFT_EDGE,gcell,primaryLeftCoord,ksizeGC)
        call Grid_getCellCoords(KAXIS,blockList(blk),&
             RIGHT_EDGE,gcell,primaryRghtCoord,ksizeGC)
        call Grid_getCellCoords(JAXIS,blockList(blk),&
             CENTER,gcell,thirdCoord,jsizeGC)
        numGuard = blkLimits(LOW,KAXIS)-1
     end if
        
     !! Restore later when implemented, this is meant for moving grid.
     !! call Grid_getCellCoords(iugrid, iXCOORD, blockList(blk), ugrid)


     call Grid_getBlkPtr(blockList(blk),solnData)

     if (hy_cvisc .ne. 0.0) then !! zero all velocity in transverse directions
        if (NDIM == 1) then
           solnData(VELY_VAR:VELZ_VAR,:,:,:)=0.0
        else if (NDIM == 2) then
           solnData(VELZ_VAR,:,:,:) = 0.0
        end if
     else   !! zero velocity on guardcells of transverse directions
        if (NDIM == 1) then
           do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
              do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
                 do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                    if(isGc(i,j,k))solnData(VELY_VAR:VELZ_VAR,i,j,k)=0.0
                 end do
              end do
           end do
        else if (NDIM == 2) then
           do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
              do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
                 do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                    if(isGc(i,j,k))solnData(VELZ_VAR,i,j,k)=0.0
                 end do
              end do
           end do
        end if
     end if

!! ---------Setting up blkLimits to call EOS
     call Timers_start("eos gc")
     eosRange = blkLimits

     if(sweepDir==SWEEP_X) then
        eosRange(LOW,IAXIS) = blkLimitsGC(LOW,IAXIS)
        eosRange(HIGH,IAXIS) = blkLimits(LOW,IAXIS)-1
     else if(sweepDir==SWEEP_Y) then
        eosRange(LOW,JAXIS) = blkLimitsGC(LOW,JAXIS)
        eosRange(HIGH,JAXIS) = blkLimits(LOW,JAXIS)-1
     else if(sweepDir==SWEEP_Z) then
        eosRange(LOW,KAXIS) = blkLimitsGC(LOW,KAXIS)
        eosRange(HIGH,KAXIS) = blkLimits(LOW,KAXIS)-1
     end if

     call Eos_wrapped(hy_eosMode,eosRange,blockList(blk))

     if(sweepDir==SWEEP_X) then
        eosRange(LOW,IAXIS) = blkLimits(HIGH,IAXIS)+1
        eosRange(HIGH,IAXIS) = blkLimitsGC(HIGH,IAXIS)
     else if(sweepDir==SWEEP_Y) then
        eosRange(LOW,JAXIS) = blkLimits(HIGH,JAXIS)+1
        eosRange(HIGH,JAXIS) = blkLimitsGC(HIGH,JAXIS)
     else if(sweepDir==SWEEP_Z) then
        eosRange(LOW,KAXIS) = blkLimits(HIGH,KAXIS)+1
        eosRange(HIGH,KAXIS) = blkLimitsGC(HIGH,KAXIS)
     end if

     call Eos_wrapped(hy_eosMode,eosRange,blockList(blk))
     call Timers_stop("eos gc")

!! ---- SAVE OLD TEMPERATURE
     !! used in driver to limit timestep
#ifdef OTMP_SCRATCH_GRID_VAR
     size(IAXIS)=isize
     size(JAXIS)=jsize
     size(KAXIS)=ksize
     startingPos=1
     call Grid_putBlkData(blockList(blk),SCRATCH,OTMP_SCRATCH_GRID_VAR,&
                          INTERIOR, startingPos, &
          solnData(TEMP_VAR, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                        blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                        blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)),size)
#endif
!!--------- DO HYDRO ON A BLOCK

!! Insert call to Grid_markCustomRegion here - SMC
     call Grid_markCustomRegion(blocklist(blk),blkLimitsGC, regionType)
     call Grid_fillCustomRegion(blocklist(blk),blkLimits,blkLimitsGC, regionType,sweepDir)

!! inserting check on gr_doHydro as well - SMC
     if ( (hy_updateHydroFluxes .OR. hy_useDiffuse) .AND. ANY(regionType == GR_NORMAL) ) then

        shock(:,:,:) = 0.0
        if ((hy_hybridRiemann .AND. hy_updateHydroFluxes).OR.hy_alwaysCallDetectShock) then
           call Hydro_detectShock(solnData, shock,blkLimits,blkLimitsGC,&
                primaryCoord,secondCoord,thirdCoord)
        end if

        ! ADDED BY JFG TO ACCOMODATE TIME-DEPENDENT DENSITY CALCULATION
        tempDens = solnData(DENS_VAR,:,:,:)
#ifdef GPO2_VAR
        tempDens2 = solnData(ODEN_VAR,:,:,:)
#endif

        call hy_ppm_block(hy_meshMe, blockList(blk),sweepDir, dt, dtOld, &
                          blkLimits,blkLimitsGC,bcs,          &
                          numCells,numGuard,      &
                          primaryCoord ,     &
                          primaryLeftCoord , &
                          primaryRghtCoord , &
                          primaryDx        , &
                          secondCoord      , &
                          thirdCoord       , &
                          radialCoord     , &
                          ugrid            , &
                          tempArea,                   &
                          tempGrav1d_o,               &
                          tempGrav1d,                 &
                          tempDtDx,                   &
                          tempFict,                   &
                          tempAreaLeft,               &
                          tempFlx,       & 
                          shock, solnData)
        
        ! ADDED BY JFG TO ACCOMODATE TIME-DEPENDENT DENSITY CALCULATION
        solnData(ODEN_VAR,:,:,:) = tempDens
#ifdef GPO2_VAR
        solnData(ODE2_VAR,:,:,:) = tempDens2
#endif

     else
        tempFlx = 0.0

        !The following arrays are zeroed to avoid computations with
        !uninitialized values in hy_pmm_updateSoln.
        tempDtDx = 0.0
        tempGrav1d_o = 0.0
        tempGrav1d = 0.0
        tempFict = 0.0
        tempArea = 0.0

        !The following array is set to one so that Grid_getFluxData does not
        !divide by zero.
        tempAreaLeft = 1.0

     end if



!!--------- END OF HYDRO ON A BLOCK 

!!--------- Do Diffuse if desired: Give FLUX-BASED implementations the opportunity to
!!          add to the the flux arrays.  
!!          (These may actually be the only contributions to the flux arrays if
!!          hy_updateHydroFluxes is false.)

!! FORFUTURE: Diffuse Unit not fully implemented

     call Diffuse_therm(sweepDir, igeom, blockList(blk), numCells,&
                        blkLimits, blkLimitsGC, primaryLeftCoord, &
                        primaryRghtCoord, tempFlx, tempAreaLeft)

     call Diffuse_visc (sweepDir, igeom, blockList(blk), numCells,&
                        blkLimits, blkLimitsGC, primaryLeftCoord,primaryRghtCoord,&
                        tempFlx, tempAreaLeft,secondCoord,thirdCoord)
     
!!$     call Diffuse_species(sweepDir, igeom, blockList(blk), numCells,&
!!$                        blkLimits, blkLimitsGC, primaryLeftCoord,primaryRghtCoord,&
!!$                        tempFlx, tempFly, tempFlz)

!!--------- Updating the solution for this one block 
!!--------- update all values if there is no flux correction
!!--------- otherwise update only the interior

     if (hy_updateHydroFluxes .OR. hy_useDiffuse) then
        call Timers_start("hy_ppm_updateSoln")
        call hy_ppm_updateSoln(updateMode, &
                        sweepDir, dt,&
                        blkLimits,blkLimitsGC,numCells, &
                        tempArea, tempGrav1d_o, tempGrav1d, &
                        tempDtDx, tempFict,   &
                        tempFlx,       &
                        solnData )
        call Timers_stop("hy_ppm_updateSoln")
     end if


     if(.not. (doFluxCorrect)) then
        if(hy_irenorm==1) then
           call Grid_renormAbundance(blockList(blk),blkLimits,solnData)
        else
           call Grid_limitAbundance(blkLimits,solnData)
        end if

!!!        call Hydro_recalibrateEints(blkLimits, blockList(blk))

        call Timers_start("eos")
        call Eos_wrapped(hy_eosMode, blkLimits, blockList(blk))
        call Timers_stop("eos")
     end if

     call Grid_releaseBlkPtr(blockList(blk),solnData)


!!--------- FLUX CONSERVATION
     if (doFluxCorrect) then

        size(1) = isizeGC
        size(2) = jsizeGC
        size(3) = ksizeGC

        call Grid_putFluxData(blockList(blk),sweepDir,tempFlx,size, hy_specialFluxVars, tempAreaLeft)
        call hy_ppm_putTemporaryData(sweepDir,blockList(blk),size,&
             tempArea, tempDtDx, tempGrav1d, tempGrav1d_o, tempFict, tempAreaLeft)

     end if

     ! Insert call to updateCustomRegion here.  This may need to be moved (?) - SMC
     call Grid_updateCustomRegion(blocklist(blk),blkLimitsGC, regionType)

#ifndef FIXEDBLOCKSIZE
     deallocate(tempFlx)
     
     deallocate(tempArea)
     deallocate(tempGrav1d_o)
     deallocate(tempGrav1d)
     deallocate(tempDtDx)
     deallocate(tempFict)
     deallocate(tempAreaLeft)
     
     deallocate(shock)
     deallocate(regionType)     
     deallocate(primaryCoord)
     deallocate(primaryLeftCoord)
     deallocate(primaryRghtCoord)
     deallocate(primaryDx)
     deallocate(secondCoord)
     deallocate(thirdCoord)
     deallocate(radialCoord)
     deallocate(ugrid)
     deallocate(isGc)
#endif

  end do




! Do this part only if refining and flux correcting

  if(doFluxCorrect) then

     call Timers_start("Grid_conserveFluxes")
     call Grid_conserveFluxes(sweepDir, level)
     call Timers_stop("Grid_conserveFluxes")

     updateMode = UPDATE_BOUND


     do blk = 1,blockCount


        call Grid_getBlkIndexLimits(blockList(blk),blkLimits,blkLimitsGC)
        
        isizeGC=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jsizeGC=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        ksizeGC=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
        
#ifndef FIXEDBLOCKSIZE
        
        numCells = max(isizeGC,jsizeGC)
        numcells = max(numCells,ksizeGC)
        
        allocate(tempFlx(NFLUXES,isizeGC,jsizeGC,ksizeGC))
        
        allocate(tempArea(isizeGC,jsizeGC,ksizeGC))
        allocate(tempGrav1d_o(isizeGC,jsizeGC,ksizeGC))
        allocate(tempGrav1d(isizeGC,jsizeGC,ksizeGC))
        allocate(tempDtDx(isizeGC,jsizeGC,ksizeGC))
        allocate(tempFict(isizeGC,jsizeGC,ksizeGC))
        allocate(regionType(isizeGC,jsizeGC,ksizeGC))
        allocate(tempAreaLeft(isizeGC,jsizeGC,ksizeGC))
#endif
        call Grid_markCustomRegion(blockList(blk),blkLimitsGC,regionType)
        call Grid_fillCustomRegion(blockList(blk),blkLimits,blkLimitsGC,regionType,sweepDir)
        
        call hy_ppm_getTemporaryData(sweepDir,blockList(blk),size,&
             tempArea,tempDtDx,tempGrav1d,tempGrav1d_o,tempFict, tempAreaLeft)
        call Grid_getFluxData(blockList(blk),sweepDir, tempFlx,size, hy_specialFluxVars, tempAreaLeft)
        
        
        call Grid_getBlkPtr(blockList(blk),solnData)
        
        call Timers_start("hy_ppm_updateSoln")
        
        call hy_ppm_updateSoln(updateMode, &
             sweepDir, dt,&
             blkLimits,blkLimitsGC,numCells, &
             tempArea, tempGrav1d_o, tempGrav1d, &
             tempDtDx, tempFict,   &
             tempFlx,   &
             solnData )
        call Timers_stop("hy_ppm_updateSoln")
        
        if(hy_irenorm==1) then
           call Grid_renormAbundance(blockList(blk),blkLimits,solnData)
        else
           call Grid_limitAbundance(blkLimits,solnData)
        end if
        
!        call IO_writeCheckpoint()

!!!        call Hydro_recalibrateEints(blkLimits, blockList(blk))

        call Timers_start("eos")
        call Eos_wrapped(hy_eosMode, blkLimits, blockList(blk))
        call Timers_stop("eos")
        
        call Grid_releaseBlkPtr(blockList(blk),solnData)
        call Grid_updateCustomRegion(blockList(blk),blkLimitsGC,regionType)

#ifndef FIXEDBLOCKSIZE
        deallocate(tempFlx)
        deallocate(tempArea)
        deallocate(tempGrav1d_o)
        deallocate(tempGrav1d)
        deallocate(tempDtDx)
        deallocate(tempFict)
        deallocate(regionType)
        deallocate(tempAreaLeft)
#endif

     end do !!END LOOP OVER BLOCKS

  end if
  
  ! Here we must call moveCustomRegion
  !! Add call to move custom region - SMC
!  call Grid_moveCustomRegion( blockCount, blockList, &
!       timeEndAdv, dt, dtOld)

  call Timers_stop("hy_ppm_sweep")

!  call IO_writeCheckpoint()
  
end subroutine hy_ppm_sweep


!! ---------Setting up block limits to call EOS

!!$     call Timers_start("eos gc")
!!$
!!$     eosRange = blkLimits
!!$
!!$     if(sweepDir==SWEEP_X) then
!!$        eosRange(LOW,IAXIS) = blkLimitsGC(LOW,IAXIS)
!!$        eosRange(HIGH,IAXIS) = blkLimits(LOW,IAXIS)-1
!!$     else if(sweepDir==SWEEP_Y) then
!!$        eosRange(LOW,JAXIS) = blkLimitsGC(LOW,JAXIS)
!!$        eosRange(HIGH,JAXIS) = blkLimits(LOW,JAXIS)-1
!!$     else if(sweepDir==SWEEP_Z) then
!!$        eosRange(LOW,KAXIS) = blkLimitsGC(LOW,KAXIS)
!!$        eosRange(HIGH,KAXIS) = blkLimits(LOW,KAXIS)-1
!!$     end if

!!$     call Eos_wrapped(hy_eosMode,eosRange,blockList(blk))


!!$     if(sweepDir==SWEEP_X) then
!!$        eosRange(LOW,IAXIS) = blkLimits(HIGH,IAXIS)+1
!!$        eosRange(HIGH,IAXIS) = blkLimitsGC(HIGH,IAXIS)
!!$     else if(sweepDir==SWEEP_Y) then
!!$        eosRange(LOW,JAXIS) = blkLimits(HIGH,JAXIS)+1
!!$        eosRange(HIGH,JAXIS) = blkLimitsGC(HIGH,JAXIS)
!!$     else if(sweepDir==SWEEP_Z) then
!!$        eosRange(LOW,KAXIS) = blkLimits(HIGH,KAXIS)+1
!!$        eosRange(HIGH,KAXIS) = blkLimitsGC(HIGH,KAXIS)
!!$     end if
!!$
!!$     if (hy_cvisc .eq. 0.0) then
!!$        if (NDIM == 1) then
!!$           solnData(VELY_VAR:VELZ_VAR, &
!!$                eosRange(LOW,IAXIS):eosRange(HIGH,IAXIS), &
!!$                eosRange(LOW,JAXIS):eosRange(HIGH,JAXIS), &
!!$                eosRange(LOW,KAXIS):eosRange(HIGH,KAXIS)) = 0.0
!!$        else if (NDIM == 2) then
!!$           solnData(VELZ_VAR, &
!!$                eosRange(LOW,IAXIS):eosRange(HIGH,IAXIS), &
!!$                eosRange(LOW,JAXIS):eosRange(HIGH,JAXIS), &
!!$                eosRange(LOW,KAXIS):eosRange(HIGH,KAXIS)) = 0.0
!!$        end if
!!$     end if
!!$     call Eos_wrapped(hy_eosMode,eosRange,blockList(blk))

!!$     call Timers_stop("eos gc")
