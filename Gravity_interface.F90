!!****h* source/physics/Gravity/Gravity_interface
!!
!! This is the header file for the gravity module that defines its
!! public interfaces.
!!***

Module Gravity_interface

#include "Flash.h"
#include "constants.h"

  interface
     subroutine Gravity_sendOutputData()
       implicit none
     end subroutine Gravity_sendOutputData
  end interface

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
     subroutine Gravity_accelOneRow (pos,sweepDir,blockID, numCells, grav, ptgrav, &
          varIndex)
       integer, intent(IN) :: sweepDir,blockID,numCells
       integer, dimension(2),INTENT(in) ::pos
       real, dimension(numCells),INTENT(inout) :: grav
       real, dimension(numCells),INTENT(inout) :: ptgrav
       integer, intent(IN), optional :: varIndex 
     end subroutine Gravity_accelOneRow
  end interface

  interface Gravity_computeDt
     subroutine Gravity_computeDt (block_no, myPE, &
          blkLimits,blkLimitsGC,  &
          solnData,   &
          dt_check, dt_minloc )
       
       integer, intent(IN) :: block_no, myPE
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
       real,INTENT(INOUT)    :: dt_check
       integer,INTENT(INOUT)    :: dt_minloc(5)
       real, pointer :: solnData(:,:,:,:) 
     end subroutine Gravity_computeDt
  end interface

  interface Gravity_finalize
     subroutine Gravity_finalize()
     end subroutine Gravity_finalize
  end interface

  interface Gravity_init
     subroutine Gravity_init(myPE)
       integer, intent(IN) :: myPE
     end subroutine Gravity_init
  end interface

  interface Gravity_potentialListOfBlocks
     subroutine Gravity_potentialListOfBlocks(blockCount,blockList)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
     end subroutine Gravity_potentialListOfBlocks
  end interface

  interface Gravity_unitTest
     subroutine Gravity_unitTest(myPE, fileUnit, perfect)
       implicit none
       integer, intent(in) :: myPE, fileUnit
       logical, intent(out) :: perfect
     end subroutine Gravity_unitTest
  end interface


end Module Gravity_interface
