!!****h* source/physics/Gravity/Gravity_interface
!!
!! This is the header file for the gravity module that defines its
!! public interfaces.
!!***

Module Gravity_interface

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
          varIndex, t_in)
       integer, intent(IN) :: sweepDir,blockID,numCells
       integer, dimension(2),INTENT(in) ::pos
       real, dimension(numCells),INTENT(inout) :: grav
       integer, intent(IN), optional :: varIndex 
       real, intent(IN), optional :: t_in
     end subroutine Gravity_accelOneRow
  end interface

  interface Gravity_computeDt
     subroutine Gravity_computeDt (blockID, myPE, dt_grav, dt_minloc)
       real,intent(OUT)       ::  dt_grav
       integer, intent(IN)    ::  blockID, myPE
       integer, intent(INOUT) :: dt_minloc(5)
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
