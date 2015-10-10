!!****h* source/physics/Gravity/Gravity_interface
!!
!! This is the header file for the gravity module that defines its
!! public interfaces.
!!***

Module Gravity_interface
#include "constants.h"
#include "Flash.h"

  interface
     subroutine Gravity_accelAtCoords (numPoints, iCoords,jCoords,kCoords, accelDir,&
          accel, blockID, &
          potentialIndex)
       integer, intent(IN) :: accelDir, numPoints
       real, dimension(:),INTENT(in) :: iCoords,jCoords,kCoords
       real, dimension(numPoints),INTENT(OUT) :: accel
       integer, intent(IN),optional :: blockID
       integer, intent(IN),optional :: potentialIndex
     end subroutine Gravity_accelAtCoords
  end interface

  interface Gravity_accelListOfBlocks
     subroutine Gravity_accelListOfBlocks (blockCount,blockList,component, &
          accelIndex, potentialIndex)
       integer,intent(IN)                      :: blockCount
       integer,dimension(blockCount), intent(IN)     :: blockList
       integer, INTENT(in) ::  component
       integer, intent(in) :: accelIndex
       integer, intent(IN), optional :: potentialIndex
     end subroutine Gravity_accelListOfBlocks
  end interface

  interface
     subroutine Gravity_accelOneRow (pos,sweepDir,blockID, numCells, grav, &
           varIndex, extraAccelVars)
       implicit none
       integer, intent(IN) :: sweepDir,blockID,numCells
       integer, dimension(2),INTENT(in) ::pos
       real, dimension(numCells),INTENT(inout) :: grav
       integer, intent(IN), optional :: varIndex 
       integer, intent(IN),OPTIONAL      :: extraAccelVars(MDIM)
     end subroutine Gravity_accelOneRow
  end interface

  interface Gravity_accelOneBlock
     subroutine Gravity_accelOneBlock ( blockID, ngcellcomp, gvec, potentialIndex)
        integer, intent(in)             :: blockID,  ngcellcomp
        real, dimension(:,:,:,:),intent(out)  :: gvec
        integer, intent(in),optional    :: potentialIndex
     end subroutine Gravity_accelOneBlock
  end interface


  interface Gravity_computeDt
     subroutine Gravity_computeDt (blockID, dx, dy, dz, &
             blkLimits, blkLimitsGC, solnData, dt_grav, dt_minloc)
       real,intent(OUT)       ::  dt_grav
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: dx
       real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: dy
       real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: dz
#else
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: dx
       real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: dy
       real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: dz
#endif
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
       real, pointer :: solnData(:,:,:,:) 
       integer, intent(IN)    ::  blockID
       integer, intent(INOUT) :: dt_minloc(5)
     end subroutine Gravity_computeDt
  end interface

  interface Gravity_finalize
     subroutine Gravity_finalize()
     end subroutine Gravity_finalize
  end interface

  interface Gravity_init
     subroutine Gravity_init()
       
     end subroutine Gravity_init
  end interface

  interface Gravity_potentialListOfBlocks
     subroutine Gravity_potentialListOfBlocks(blockCount,blockList, &
          potentialIndex)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       integer, intent(IN), optional :: potentialIndex
     end subroutine Gravity_potentialListOfBlocks
  end interface

  interface Gravity_unitTest
     subroutine Gravity_unitTest( fileUnit, perfect)
       implicit none
       integer, intent(in) :: fileUnit
       logical, intent(out) :: perfect
     end subroutine Gravity_unitTest
  end interface

! The following subroutines and functions are needed by the tree solver
#include "FortranLangFeatures.fh"
  interface Gravity_bhAccBotNode
    subroutine Gravity_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, dimension(:,:,:) :: locCoords
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(IN) :: botnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine Gravity_bhAccBotNode
  end interface

  interface Gravity_bhAccNode
    subroutine Gravity_bhAccNode(subnode, accnode)
      real, dimension(:), intent(IN)  :: subnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine Gravity_bhAccNode
  end interface

  interface Gravity_bhBotNodeContrib
    subroutine Gravity_bhBotNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
      real, dimension(:), intent(IN) :: node
      real, intent(IN) :: ndSize
      real, dimension(MDIM+2), intent(IN) :: dr
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine Gravity_bhBotNodeContrib
  end interface


  interface Gravity_bhFillBotNode
    subroutine Gravity_bhFillBotNode(blockno, point, solnData, botnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(INOUT) :: botnode
    end subroutine Gravity_bhFillBotNode
  end interface

  interface Gravity_bhNodeContrib
    subroutine Gravity_bhNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
      real, dimension(:), intent(IN) :: node
      real, intent(IN) :: ndSize
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      real, dimension(MDIM+2), intent(IN) :: dr
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine Gravity_bhNodeContrib
  end interface

  interface Gravity_bhNormalizeNode
    subroutine Gravity_bhNormalizeNode(smr, node)
      real, dimension(MDIM), intent(IN) :: smr
      real, dimension(:), intent(INOUT) :: node
    end subroutine Gravity_bhNormalizeNode
  end interface

  interface Gravity_bhPostprocNode
    subroutine Gravity_bhPostprocNode(ndSize, node)
      real, intent(IN) :: ndSize
      real, dimension(:), intent(INOUT) :: node
    end subroutine Gravity_bhPostprocNode
  end interface

  interface Gravity_bhStartBlock
    subroutine Gravity_bhStartBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine Gravity_bhStartBlock
  end interface

  interface Gravity_bhFinalizeBlock
    subroutine Gravity_bhFinalizeBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER :: solnData
    end subroutine Gravity_bhFinalizeBlock
  end interface


  interface Gravity_bhInitFieldVar
    subroutine Gravity_bhInitFieldVar(gpotVar)
      implicit none
      integer, intent(IN) :: gpotVar
    end subroutine Gravity_bhInitFieldVar
  end interface

  interface Gravity_bhMAC
    logical function Gravity_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize2
    integer, dimension(MDIM), intent(IN) :: point
    real, dimension(MDIM+2), intent(IN) :: dr
    integer, intent(IN) :: blockno
    integer, dimension(2,MDIM)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER :: solnData
    end function Gravity_bhMAC
  end interface
  
  interface Gravity_bhPartErr
    subroutine Gravity_bhPartErr(node, ndSize, dr, perr)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize
    real, dimension(MDIM+2), intent(IN) :: dr
    real :: perr
    end subroutine Gravity_bhPartErr
  end interface
  
  interface Gravity_bhSelfContrib
    subroutine Gravity_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
      implicit none
      real, dimension(:), intent(IN) :: node
      real, dimension(MDIM), intent(IN) :: cellsize
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine Gravity_bhSelfContrib
  end interface

  interface Gravity_bhGenEwaldField
    subroutine Gravity_bhGenEwaldField()
    end subroutine Gravity_bhGenEwaldField
  end interface

  interface Gravity_bhGenEwaldFieldLevel
    subroutine Gravity_bhGenEwaldFieldLevel(nx, ny, nz, ewald_xmax, ewald_ymax, ewald_zmax, &
           & ewald_periodicity, field)
      implicit none
      integer, intent(IN) :: ewald_periodicity, nx, ny, nz
      real, intent(IN) ::  ewald_xmax, ewald_ymax, ewald_zmax
      real , intent(INOUT) :: field(-1:,-1:,-1:)
    end subroutine Gravity_bhGenEwaldFieldLevel
  end interface

  interface Gravity_bhEwald
    real function Gravity_bhEwald(x, y, z)
      implicit none
      real,intent(in)    :: x, y, z
    end function Gravity_bhEwald
  end interface

  interface Gravity_bhGetNodeStruct
    subroutine Gravity_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
      implicit none
      integer, intent(IN) :: im, ix, iy, iz
      integer, intent(INOUT) :: bnsize, nsize
    end subroutine Gravity_bhGetNodeStruct
  end interface




! end of tree solver subroutines/functions



end Module Gravity_interface
