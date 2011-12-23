!!****h* source/Grid/Grid_interface
!!
!! This is the header file for the grid module that defines its
!! public interfaces.
!!***
Module Grid_interface

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
!#include "GridParticles.h"

  integer,parameter :: GRID_PDE_BND_ISOLATED  = 0
  integer,parameter :: GRID_PDE_BND_PERIODIC  = 1
  integer,parameter :: GRID_PDE_BND_DIRICHLET = 2
  integer,parameter :: GRID_PDE_BND_NEUMANN   = 3
  integer,parameter :: GRID_PDE_BND_GIVENVAL  = 4
  integer,parameter :: GRID_PDE_BND_GIVENGRAD = 5

  interface
     subroutine Grid_findExtrema(blockID,ivar,flag,extrema)
       integer,intent(IN):: blockID,ivar,flag
       real,intent(INOUT):: extrema
     end subroutine Grid_findExtrema
  end interface

  interface
     subroutine Grid_applyBCEdge(bcType,bcDir,guard,var,dataRow,face,&
          gridDataStruct, blockHandle, secondCoord, thirdCoord)
       integer,intent(IN):: bcType,bcDir,guard,var,face,gridDataStruct
       real,dimension(:),intent(INOUT)::dataRow
       integer,intent(IN),OPTIONAL:: blockHandle
       real,intent(IN),OPTIONAL :: secondCoord,thirdCoord
     end subroutine Grid_applyBCEdge
  end interface

  interface
     subroutine Grid_applyBCEdgeAllUnkVars(bcType,bcDir,guard,dataRow,face,&
          cellCenterSweepCoord, secondCoord,thirdCoord, blockHandle)
       integer,intent(IN):: bcType,bcDir,guard,face
       real,dimension(2*guard,NUNK_VARS),intent(INOUT)::dataRow
       real,intent(IN):: cellCenterSweepCoord(*), secondCoord,thirdCoord
       integer,intent(IN),OPTIONAL:: blockHandle
     end subroutine Grid_applyBCEdgeAllUnkVars
  end interface

  interface Grid_computeUserVars
     subroutine Grid_computeUserVars()
     end subroutine Grid_computeUserVars
  end interface

  interface
     subroutine Grid_computeVarNorm (level, normType, ivar, norm, leafOnly)
       implicit none
       integer, intent(IN)  :: normType, level, ivar, leafOnly
       real, intent(OUT)    :: norm
     end subroutine Grid_computeVarNorm
  end interface

  interface
     subroutine Grid_computeVarDiff(level, gr_iRefSoln, gr_iSoln, ires)
       implicit none
       integer, intent(in)          :: level, gr_iRefSoln, gr_iSoln, ires
     end subroutine Grid_computeVarDiff
  end interface

  interface Grid_conserveFluxes
     subroutine Grid_conserveFluxes( axis, level)
       integer, intent(in) :: axis, level
     end subroutine Grid_conserveFluxes
  end interface

  interface Grid_dump
     subroutine Grid_dump(var,num,blockID,gcell)
       integer, intent(IN) :: num, blockID
       integer, dimension(num), intent(IN) :: var
       logical, intent(IN) :: gcell
     end subroutine Grid_dump
  end interface


  interface
     subroutine Grid_fillGuardCells( gridDataStruct,idir,&
          minLayers,eosMode,doEos,maskSize,mask,makeMaskConsistent,&
          selectBlockType)
       integer, intent(in) :: gridDataStruct
       integer, intent(in) :: idir
       integer,optional,intent(in) :: minLayers
       integer,optional,intent(in) :: eosMode
       logical,optional,intent(IN) :: doEos
       integer,optional,intent(in) :: maskSize
       logical,optional,dimension(:), intent(IN) :: mask
       logical,optional, intent(IN) :: makeMaskConsistent
       integer,optional,intent(in) :: selectBlockType
     end subroutine Grid_fillGuardCells
  end interface

  interface Grid_finalize
     subroutine Grid_finalize()
     end subroutine Grid_finalize
  end interface

  interface Grid_getBlkBC
     subroutine Grid_getBlkBC(blockId, faces, onBoundary)
       integer, intent(in) :: blockId
       integer, dimension(2,MDIM),intent(out):: faces
       integer, optional, dimension(2,MDIM), intent(out) :: onBoundary
     end subroutine Grid_getBlkBC
  end interface

  interface Grid_getBlkBoundBox
     subroutine Grid_getBlkBoundBox(blockId,boundBox)
       integer, intent(in) :: blockId
       real, dimension(2, MDIM), intent(out) :: boundBox
     end subroutine Grid_getBlkBoundBox
  end interface

  interface Grid_getBlkCenterCoords
     subroutine Grid_getBlkCenterCoords(blockId, blockCenter)
       integer,intent(in) :: blockId
       real,dimension(MDIM),intent(out) :: blockCenter
     end subroutine Grid_getBlkCenterCoords
  end interface

  interface Grid_getBlkCornerID
     subroutine Grid_getBlkCornerID(blockId, cornerID, stride)
       integer,intent(IN)  :: blockId
       integer,dimension(MDIM), intent(OUT) :: cornerID, stride
     end subroutine Grid_getBlkCornerID
  end interface

  interface
     subroutine Grid_getBlkData(blockID, dataType, structIndex, beginCount, &
          startingPos, datablock, dataSize)
       integer, intent(in) :: blockID, structIndex, beginCount, dataType
       integer, dimension(MDIM), intent(in) :: startingPos
       integer, dimension(3), intent(in) :: dataSize
       real, dimension(datasize(1), dataSize(2), dataSize(3)),intent(out) :: datablock
     end subroutine Grid_getBlkData
  end interface

  interface
     subroutine Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC,gridDataStruct)
       integer,intent(IN)                     :: blockId
       integer,dimension(2,MDIM), intent(OUT) :: blkLimits, blkLimitsGC
       integer,optional,intent(IN)            :: gridDataStruct
     end subroutine Grid_getBlkIndexLimits
  end interface

  interface Grid_getBlkPhysicalSize
     subroutine Grid_getBlkPhysicalSize(blockId, blockSize)
       integer,intent(in) :: blockId
       real,dimension(MDIM),intent(out) :: blockSize
     end subroutine Grid_getBlkPhysicalSize
  end interface

  interface
     subroutine Grid_getBlkPtr(blockId, dataPtr,gridDataStruct)
       integer, intent(in) :: blockId
       real,dimension(:,:,:,:), pointer :: dataPtr
       integer,optional, intent(in) :: gridDataStruct
     end subroutine Grid_getBlkPtr
  end interface

  interface Grid_getBlkRefineLevel
     subroutine Grid_getBlkRefineLevel(blockID, refineLevel)
       integer,intent(in) :: blockID
       integer,intent(out) :: refineLevel
     end subroutine Grid_getBlkRefineLevel
  end interface

  interface
     subroutine Grid_getCellCoords(axis, blockID, edge, guardcell, coordinates, size)
       integer, intent(in) :: axis, blockID, edge
       integer, intent(in) :: size
       logical, intent(in) :: guardcell
       real,intent(out), dimension(size) :: coordinates
     end subroutine Grid_getCellCoords
  end interface

  interface Grid_getDeltas
     subroutine Grid_getDeltas(blockId, del)
       integer, intent(in) :: blockId
       real, dimension(MDIM), intent(out) :: del
     end subroutine Grid_getDeltas
  end interface

  interface
     subroutine Grid_getFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)
       integer, intent(IN) :: blockID
       integer, intent(IN) :: axis
       integer, intent(IN), dimension(3) :: dataSize
       real, intent(INOUT), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
       integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
       real, intent(IN), OPTIONAL :: areaLeft(:,:,:)
     end subroutine Grid_getFluxData
  end interface

  interface
     subroutine Grid_getGlobalIndexLimits(globalIndexLimits)
       integer, dimension(MDIM), intent(out) :: globalIndexLimits
     end subroutine Grid_getGlobalIndexLimits
  end interface

  interface
     subroutine Grid_getListOfBlocks(blockType, listOfBlocks,count,refinementLevel)
       integer, intent(in) :: blockType
       integer,dimension(1),intent(out) :: listOfBlocks
       integer,intent(out) :: count
       integer,intent(IN),optional :: refinementLevel
     end subroutine Grid_getListOfBlocks
  end interface

  interface
     subroutine Grid_getLocalNumBlks(numBlocks)
       integer,intent(out) :: numBlocks
     end subroutine Grid_getLocalNumBlks
  end interface

  interface Grid_getMinCellSize
     subroutine Grid_getMinCellSize(minCellSize)
       real, intent(OUT) :: minCellSize
     end subroutine Grid_getMinCellSize
  end interface

  interface
     subroutine Grid_getMinCellSizes(minCellSizes)
       real, dimension(MDIM), intent(OUT) :: minCellSizes
     end subroutine Grid_getMinCellSizes
  end interface

  interface Grid_getPlaneData
     subroutine Grid_getPlaneData(blockid, gridDataStruct, variable, beginCount, &
          plane, startingPos, datablock, dataSize)
       integer, intent(IN) :: blockid, variable, beginCount, plane, gridDataStruct
       integer, dimension(MDIM), intent(IN) :: startingPos
       integer, dimension(2), intent(IN) :: dataSize
       real, dimension(datasize(1), dataSize(2)),intent(OUT) :: datablock
     end subroutine Grid_getPlaneData
  end interface

  interface Grid_getPointData
     subroutine Grid_getPointData(blockID, gridDataStruct, variable, beginCount, &
          position, datablock)
       integer, intent(in) :: blockID, variable, beginCount, gridDataStruct
       integer, dimension(MDIM), intent(in) :: position
       real, intent(out) :: datablock
     end subroutine Grid_getPointData
  end interface

  interface Grid_getRowData
     subroutine Grid_getRowData(blockid, gridDataStruct, variable, beginCount, &
          row, startingPos, datablock, dataSize)
       integer, intent(in) :: blockid, variable, beginCount, row, gridDataStruct
       integer, dimension(MDIM), intent(in) :: startingPos
       integer, intent(in) :: dataSize
       real, dimension(datasize),intent(out) :: datablock
     end subroutine Grid_getRowData
  end interface

  interface Grid_getSingleCellCoords
     subroutine Grid_getSingleCellCoords(ind, blockId,edge, beginCount,coords)
       integer,dimension(MDIM), intent(in) :: ind
       integer, intent(in) :: blockId, edge
       integer, intent(in) :: beginCount
       real, dimension(MDIM), intent(out) :: coords
     end subroutine Grid_getSingleCellCoords
  end interface

  interface Grid_getSingleCellVol
     subroutine Grid_getSingleCellVol(blockID, beginCount, point, cellvolume)
       integer, intent(in) :: blockID, beginCount
       integer, intent(in) :: point(MDIM)
       real, intent(out)   :: cellvolume
     end subroutine Grid_getSingleCellVol
  end interface
  
  interface Grid_init
     subroutine Grid_init()
     end subroutine Grid_init
  end interface

  interface Grid_initDomain
     subroutine Grid_initDomain( restart,particlesInitialized)
       logical, intent(IN) :: restart
       logical,intent(INOUT) :: particlesInitialized
     end subroutine Grid_initDomain
  end interface

  interface Grid_markBlkDerefine
     subroutine Grid_markBlkDerefine(block,mark)
       integer, intent(IN) :: block
       logical, intent(IN) :: mark
     end subroutine Grid_markBlkDerefine
  end interface

  interface Grid_markBlkRefine
     subroutine Grid_markBlkRefine(block,mark)
       integer, intent(IN) :: block
       logical, intent(IN) :: mark
     end subroutine Grid_markBlkRefine
  end interface

  interface Grid_markRefineDerefine
     subroutine Grid_markRefineDerefine()
     end subroutine Grid_markRefineDerefine
  end interface

  interface Grid_markRefineSpecialized
     subroutine Grid_markRefineSpecialized(criterion,size,specs,lref)
       integer, intent(IN) :: criterion
       integer, intent(IN) :: size
       real,dimension(size),intent(IN) :: specs
       integer, intent(IN) ::  lref
     end subroutine Grid_markRefineSpecialized
  end interface

  interface Grid_moveParticles
     subroutine Grid_moveParticles(dataBuf, propCount, maxCount, localCount, &
       index_list, indexCount,&
       coords_in_blk)
       
       
       integer,intent(IN) :: maxCount, propCount, indexCount
       integer,intent(INOUT) :: localCount
       
       real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf
       integer, dimension(indexCount),intent(IN) :: index_list
       logical, intent(IN) :: coords_in_blk
       
     end subroutine Grid_moveParticles
  end interface

  interface
     subroutine Grid_putBlkData(blockID, gridDataStruct, variable, beginCount, &
          startingPos, datablock, dataSize)
       integer, intent(in) :: blockID, variable, beginCount, gridDataStruct
       integer, dimension(MDIM), intent(in) :: startingPos
       integer, dimension(3), intent(in) :: dataSize
       real, dimension(datasize(1), dataSize(2), dataSize(3)),intent(in) :: datablock
     end subroutine Grid_putBlkData
  end interface

  interface
     subroutine Grid_putFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)
       implicit none
       integer, intent(IN) :: blockID
       integer, intent(IN) :: axis
       integer, intent(IN), dimension(3) :: dataSize
       real, intent(IN), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
       integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
       real, intent(IN), OPTIONAL :: areaLeft(:,:,:)
     end subroutine Grid_putFluxData
  end interface

  interface Grid_putLocalNumBlks
     subroutine Grid_putLocalNumBlks(numBlocks)
       integer,intent(in) :: numBlocks
     end subroutine Grid_putLocalNumBlks
  end interface

  interface Grid_putPlaneData
     subroutine Grid_putPlaneData(blockid, gridDataStruct, variable, beginCount, &
          plane, startingPos, datablock, dataSize)
       integer, intent(in) :: blockid, variable, beginCount, plane, gridDataStruct
       integer, dimension(MDIM), intent(in) :: startingPos
       integer, dimension(2), intent(in) :: dataSize
       real, dimension(datasize(1), dataSize(2)),intent(in) :: datablock
     end subroutine Grid_putPlaneData
  end interface

  interface Grid_putPointData
     subroutine Grid_putPointData(blockid, gridDataStruct, variable, beginCount, position, datablock)
       integer, intent(in) :: blockid, variable, beginCount, gridDataStruct
       integer, dimension(MDIM), intent(in) :: position
       real, intent(in) :: datablock
     end subroutine Grid_putPointData
  end interface

  interface
     subroutine Grid_putRowData(blockid, gridDataStruct, variable, beginCount, &
          row, startingPos, datablock, dataSize)
       integer, intent(IN) :: blockid, variable, beginCount, row, gridDataStruct
       integer, dimension(MDIM), intent(IN) :: startingPos
       integer, intent(IN) :: dataSize
       real, dimension(datasize),intent(IN) :: datablock
     end subroutine Grid_putRowData
  end interface

  interface
     subroutine Grid_releaseBlkPtr(blockId, dataPtr, gridDataStruct)
       integer, intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine Grid_releaseBlkPtr
  end interface

  interface Grid_restrictAllLevels
     subroutine Grid_restrictAllLevels()
     end subroutine Grid_restrictAllLevels
  end interface

  interface
     subroutine Grid_restrictByLevels( gridDataStruct, fromLevel, toLevel, checkFinestLevel,&
          maskSize,mask)
       integer, intent(in) :: gridDataStruct
       integer, intent(in) :: fromLevel, toLevel
       logical, optional,intent(in) :: checkFinestLevel
       integer, optional,intent(in) :: maskSize
       logical,dimension(*),optional,intent(in) :: mask
     end subroutine Grid_restrictByLevels
  end interface

  interface Grid_sendOutputData
     subroutine Grid_sendOutputData()
     end subroutine Grid_sendOutputData
  end interface

  interface
     subroutine Grid_setFluxHandling(handling, status)
       implicit none
       character(len=*),intent(IN) :: handling
       integer,intent(OUT),OPTIONAL :: status
     end subroutine Grid_setFluxHandling
  end interface

  interface Grid_updateRefinement
     subroutine Grid_updateRefinement( nstep,time, gridChanged)
       integer, intent(in) :: nstep
       real, intent(in) :: time
       logical, intent(out), OPTIONAL :: gridChanged
     end subroutine Grid_updateRefinement
  end interface

  interface Grid_unitTest
     subroutine Grid_unitTest(fileUnit,perfect)

       integer, intent(in)           :: fileUnit ! Output to file
       logical, intent(inout)        :: perfect  ! Flag to indicate errors
     end subroutine Grid_unitTest
  end interface

  interface Grid_getGeometry
     subroutine Grid_getGeometry(geometry)
       integer, intent(OUT) :: geometry
     end subroutine Grid_getGeometry
  end interface


  interface
     subroutine Grid_sortParticles(dataBuf,props,localNumCount,&
          elementTypes,maxPerProc,&
          elementsPerBlk, attrib1, attrib2)

       integer,intent(INOUT) :: localNumCount
       integer,intent(IN) :: maxPerProc, props,elementTypes

       real,intent(INOUT),dimension(props,maxPerProc) :: dataBuf
       integer,intent(OUT),dimension(MAXBLOCKS,elementTypes) :: elementsPerBlk
       integer, intent(IN) :: attrib1
       integer, optional, intent(IN) :: attrib2
     end subroutine Grid_sortParticles
  end interface

  interface
     subroutine Grid_mapMeshToParticles (particles, part_props,part_blkID,&
                                         numParticles,&
                                         posAttrib,&
                                         numAttrib, attrib,&
                                         mapType,gridDataStruct)
       implicit none
       integer, INTENT(in) :: part_props, numParticles, part_blkID
       real, INTENT(inout),dimension(part_props,numParticles) :: particles
       integer, intent(IN) :: numAttrib
       integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
       integer,dimension(MDIM), intent(IN) :: posAttrib
       integer, INTENT(IN) :: mapType
       integer, optional, intent(IN) :: gridDataStruct
     end subroutine Grid_mapMeshToParticles
  end interface

  interface
     subroutine Grid_mapParticlesToMesh (particles,part_props,numParticles,&
          maxParticlesPerProc,propPart, varGrid, mode)
  
       integer,intent(IN) :: numParticles, part_props,maxParticlesPerProc
       real,dimension(part_props,maxParticlesPerProc),intent(INOUT) :: particles
       integer, INTENT(in) :: propPart, varGrid
       integer, INTENT(in), optional :: mode

     end subroutine Grid_mapParticlesToMesh
  end interface
  
  
  interface 
     subroutine Grid_solvePoisson (grav_iSoln, grav_iSrc, grav_bcTypes, &
          grav_bcValues, poisfact)
       implicit none
       integer, intent(in)    :: grav_iSoln, grav_iSrc
       integer, intent(in)    :: grav_bcTypes(6)
       real, intent(in)       :: grav_bcValues(2,6)
       real, intent(inout)    :: poisfact
     end subroutine Grid_solvePoisson
  end interface
  
  interface 
     subroutine Grid_advanceDiffusion (iVar, iSrc, iFactorB, iFactorA, bcTypes, bcValues, dt, chi, scaleFact, &
          theta, solnIsDelta, iFactorC, iFactorD, pass)       
       implicit none
       
       integer, intent(IN) :: iVar
       integer, intent(IN) :: iSrc
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: iFactorA
       real, intent(IN)    :: dt 
       real, intent(IN)    :: chi
       real, intent(IN)    :: scaleFact
       real, intent(IN)    :: theta
       logical, intent(IN) :: solnIsDelta
       integer, dimension(6),  intent(IN) :: bcTypes
       real   , dimension(2,6),intent(IN) :: bcValues
       integer, intent(IN), OPTIONAL :: pass
       integer, intent(IN), OPTIONAL :: iFactorC
       integer, intent(IN), OPTIONAL :: iFactorD   
     end subroutine Grid_advanceDiffusion
  end interface

  interface
     subroutine Grid_limitAbundance(blkLimits,solnData)
       integer, dimension(2,MDIM), INTENT(in) :: blkLimits
       real, POINTER :: solnData(:,:,:,:)
     end subroutine Grid_limitAbundance
  end interface


  interface
     subroutine Grid_renormAbundance(blockId,blkLimits,solnData)
       integer, INTENT(in) :: blockId
       integer, intent(in), dimension(2,MDIM)::blkLimits
       real,pointer :: solnData(:,:,:,:)
     end subroutine Grid_renormAbundance
  end interface


  interface
     subroutine Grid_renormMassScalars(blkLimits,solnData)
       integer, intent(in), dimension(2,MDIM)::blkLimits
       real,pointer :: solnData(:,:,:,:)
     end subroutine Grid_renormMassScalars
  end interface

  interface
     subroutine Grid_conserveField ()
     end subroutine Grid_conserveField
  end interface


  interface
     subroutine Grid_pfft(direction,inArray,outArray)
       integer, intent(IN) :: direction
       real, dimension(:),intent(IN) :: inArray
       real,dimension(:), intent(OUT) :: outArray
     end subroutine Grid_pfft
  end interface

  interface
     subroutine Grid_pfftInit( ndim, needMap,globalLen, mapSize,&
          transformType, baseDatType, jProcs, kProcs, refinementLevel)
       integer, intent(IN) :: ndim
       logical, intent(IN) :: needMap
       integer, dimension(MDIM), intent(IN) :: globalLen
       integer,dimension(MDIM),intent(OUT) ::  mapSize
       integer,dimension(MDIM),optional,intent(IN) :: transformType
       integer,dimension(0:MDIM),optional,intent(IN) :: baseDatType
       integer,optional,intent(IN) :: jProcs, kProcs
       integer,optional,intent(IN) :: refinementLevel
     end subroutine Grid_pfftInit
  end interface

  interface
     subroutine Grid_pfftFinalize()
     end subroutine Grid_pfftFinalize
  end interface

  interface
     subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
          guard,axis,face,regionData,regionSize,mask,applied,&
          blockHandle,secondDir,ThirdDir,endPoints,blkLimitsGC, idest)

       integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
       integer,dimension(REGION_DIM),intent(IN) :: regionSize
       real,dimension(regionSize(BC_DIR),&
            regionSize(SECOND_DIR),&
            regionSize(THIRD_DIR),&
            regionSize(STRUCTSIZE)),intent(INOUT)::regionData
       logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
       logical, intent(OUT) :: applied
       integer,intent(IN) :: blockHandle
       integer,intent(IN) :: secondDir,thirdDir
       integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
       integer,intent(IN),OPTIONAL:: idest
     end subroutine Grid_bcApplyToRegionSpecialized
  end interface

  interface
     subroutine Grid_bcApplyToRegion(bcType,gridDataStruct,&
          guard,axis,face,regionData,regionSize,mask,applied,&
          blockHandle,secondDir,ThirdDir,endPoints,blkLimitsGC, idest)
       integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
       integer,dimension(REGION_DIM),intent(IN) :: regionSize
       real,dimension(regionSize(BC_DIR),&
            regionSize(SECOND_DIR),&
            regionSize(THIRD_DIR),&
            regionSize(STRUCTSIZE)),intent(INOUT)::regionData
       logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
       logical, intent(OUT) :: applied
       integer,intent(IN) :: blockHandle
       integer,intent(IN) :: secondDir,thirdDir
       integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
  integer,intent(IN),OPTIONAL:: idest

     end subroutine Grid_bcApplyToRegion
  end interface

  interface
     subroutine Grid_pfftGetIndexLimits(configLimits,phaseLimits)
       integer,dimension(LOW:HIGH,MDIM),intent(OUT) :: configLimits, phaseLimits
     end subroutine Grid_pfftGetIndexLimits
  end interface


  interface
     subroutine Grid_getBlkType(blockID, blkType)
       implicit none
       integer,intent(in) :: blockID
       integer,intent(out) :: blkType
     end subroutine Grid_getBlkType
  end interface

  interface Grid_pfftMapToInput
     subroutine Grid_pfftMapToInput(gridVar, pfft_inArray)
       integer,intent(IN) :: gridVar
       real, dimension(:),intent(OUT) :: pfft_inArray
     end subroutine Grid_pfftMapToInput
     subroutine Grid_pfftMapToInput3DArr(gridVar, pfft_inArray)
       integer,intent(IN) :: gridVar
       real, dimension(:,:,:),intent(OUT) :: pfft_inArray
     end subroutine Grid_pfftMapToInput3DArr
  end interface

  interface Grid_pfftMapFromOutput
     subroutine Grid_pfftMapFromOutput(gridVar, pfft_outArray)
       integer,intent(IN) :: gridVar
       real, dimension(:),intent(IN) :: pfft_outArray
     end subroutine Grid_pfftMapFromOutput
     subroutine Grid_pfftMapFromOutput3DArr(gridVar, pfft_outArray)
       integer,intent(IN) :: gridVar
       real, dimension(:,:,:),intent(IN) :: pfft_outArray
     end subroutine Grid_pfftMapFromOutput3DArr
  end interface

  interface
     subroutine Grid_outsideBoundBox(pos,bndBox,outside,Negh)
       real,dimension(MDIM),intent(IN) :: pos
       real,dimension(LOW:HIGH,MDIM),intent(IN) :: bndBox
       logical, intent(OUT) :: outside
       integer, dimension(MDIM),intent(OUT) :: Negh
     end subroutine Grid_outsideBoundBox
  end interface

  interface
     subroutine Grid_getMaxCommonRefinement(inputComm, maxRefinement)
       implicit none
       integer, intent(IN) :: inputComm
       integer, intent(OUT) :: maxRefinement
     end subroutine Grid_getMaxCommonRefinement
  end interface

  interface

     subroutine Grid_GCPutScratch(gridDataStruct,needSetup,releaseSetup,&
          &blkCount,blkList,indCount,indList, gcCnt, putData)

       integer, intent(IN) :: gridDataStruct
       logical, intent(IN) :: needSetup, releaseSetup
       integer, optional,intent(IN) :: blkCount, indCount
       integer, optional, dimension(:), intent(IN) :: blkList
       integer, optional, dimension(:), intent(IN) :: indList
       integer, optional, dimension(NDIM),intent(IN) :: gcCnt  
       logical, optional,intent(IN) :: putData

     end subroutine Grid_GCPutScratch

  end interface

  interface

     subroutine Grid_GCTransferOneBlk(mode,gridDataStruct,blkIndex,blkData)

       logical, intent(IN) :: mode
       integer, intent(IN) :: gridDataStruct, blkIndex
       real, pointer, dimension(:,:,:,:) :: blkData
       
     end subroutine Grid_GCTransferOneBlk

  end interface

  
  interface
     
     subroutine Grid_markCustomRegion(blkID,blkLimitsGC,regionType)
       integer, intent(in) :: blkID
       integer, intent(in), dimension(LOW:HIGH,MDIM) :: blkLimitsGC
       real, intent(out), &
            DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),                   &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  regionType
     end subroutine Grid_markCustomRegion
  end interface

  interface
     subroutine Grid_getNumVars(gridStruct, nVar)  
       implicit none
       integer, intent(in) :: gridStruct
       integer, intent(out) :: nVar
     end subroutine Grid_getNumVars
  end interface

  interface
     
     subroutine Grid_fillCustomRegion(blkID,blkLimits,blkLimitsGC,regionType,sweepDir)

       integer, intent(in) :: blkID
       integer, intent(in) :: sweepDir
       integer, intent(in), dimension(LOW:HIGH,MDIM) :: blkLimitsGC, blkLimits
          real, intent(in), DIMENSION&
            (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  regionType
     end subroutine Grid_fillCustomRegion
  end interface

  interface

     subroutine Grid_updateCustomRegion(blkID,blkLimitsGC, regionType)
  
       integer, intent(in) :: blkID
       integer, intent(in), dimension(LOW:HIGH,MDIM) :: blkLimitsGC
       real, intent(inout), DIMENSION&
            (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  regionType
     end subroutine Grid_updateCustomRegion
  end interface

  interface

     subroutine Grid_moveCustomRegion(blockCount, blockList, timeEndAdv, dt, dtOld)
  
       integer, INTENT(IN) :: blockCount
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       
     end subroutine Grid_moveCustomRegion

  end interface

  interface

     subroutine Grid_addToVar(srcVar, destVar, multFactor, reset)
       integer, intent(in) :: srcVar, destVar
       real,  intent(in) :: multFactor
       logical, intent(in) :: reset
       
     end subroutine Grid_addToVar

  end interface

  interface
     subroutine Grid_primitiveToConserve(blkList,count,force)
       integer,intent(IN) :: count
       integer,dimension(count),intent(IN) :: blkList 
       logical,intent(IN) :: force
     end subroutine Grid_primitiveToConserve
  end interface

  interface
     subroutine Grid_conserveToPrimitive(blkList,count,allCells,force)
       integer,intent(IN) :: count
       integer,dimension(count),intent(IN) :: blkList 
       logical,intent(IN) :: allCells,force
     end subroutine Grid_conserveToPrimitive
  end interface

  interface

     subroutine Grid_getDomainBoundBox(boundBox)
       real,dimension(LOW:HIGH,MDIM), intent(OUT) :: boundBox
     end subroutine Grid_getDomainBoundBox
  end interface

  interface
     subroutine Grid_getDomainBC(boundary)
       integer,dimension(LOW:HIGH,MDIM), intent(OUT) :: boundary
     end subroutine Grid_getDomainBC
  end interface

  interface
     subroutine Grid_parseNonRep(strlwr, nonrep, idx)
      implicit none
      character(*), intent(in) :: strlwr
      integer, intent(out) :: nonrep, idx
     end subroutine
  end interface
  
  interface
     subroutine Grid_formatNonRep(nonrep, idx, str)
        implicit none
        integer, intent(in) :: nonrep, idx
        character(*), intent(out) :: str
     end subroutine
  end interface
  
  interface
     subroutine Grid_getVarNonRep(mapblock, var, nonrep, locidx)
       implicit none
       integer, intent(in) :: mapblock, var
       integer, intent(out) :: nonrep
       integer, intent(out), optional :: locidx
     end subroutine
  end interface
end Module Grid_interface
