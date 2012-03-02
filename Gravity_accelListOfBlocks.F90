!!****if* source/physics/Gravity/GravityMain/Gravity_accelListOfBlocks
!!
!! NAME
!!
!!  Gravity_accelListOfBlocks  
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelListOfBlocks(integer(IN) :: blockCount,
!!                         integer(IN)    :: blockList(blockCount),
!!                         integer(IN)    :: component,
!!                         integer(IN)    :: accelIndex,
!!                         integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!!  This version implements finite-volume differencing of a
!!  potential field.
!!
!!
!! ARGUMENTS
!!
!!   blockCount   -  The number of blocks in the list
!!   blockList    -  The list of blocks on which to calculate acceleration.
!!   component    -  The component of the acceleration to compute.
!!                   Permitted values are IAXIS, JAXIS, KAXIS.  If ALLDIR
!!                      is selected, the routine aborts messily
!!   accelIndex - variable # to store the acceleration
!!   potentialIndex -   Variable # to take as potential if present
!!
!! NOTES
!!
!!   This routine can be used as a wrapper to Gravity_accelOneRow.  Each implementation
!!   of the Gravity unit has a version of Gravity_accelOneRow, but this wrapper remains
!!   constant.
!!***



subroutine Gravity_accelListOfBlocks (blockCount, blockList, component, &
     accelIndex, potentialIndex)

  use Gravity_data, ONLY : useGravity

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Gravity_interface, ONLY : Gravity_accelOneRow
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none
  
#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN)                      :: blockCount
  integer,dimension(blockCount), intent(IN)     :: blockList
  
  integer, INTENT(in) :: component, accelIndex
  integer,intent(IN),optional :: potentialIndex
  
  
  integer       :: j, k, lb, nxzones, nyzones, nzzones,blockID
  real, pointer :: solnVec(:,:,:,:)
  
#ifdef FIXEDBLOCKSIZE
  real        :: grav(max(GRID_IHI_GC,GRID_JHI_GC, GRID_KHI_GC))
  real        :: ptgrav(max(GRID_IHI_GC,GRID_JHI_GC, GRID_KHI_GC))
#else
  real,allocatable,dimension(:) :: grav
#endif
  

  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: potVar
  integer,dimension(MDIM-1) :: pos
  !=========================================================================

  if(.not.useGravity) return

  if(present(potentialIndex)) then
     potVar=potentialIndex
  else
#ifdef GPOT_VAR
     potVar=GPOT_VAR
#else
     potVar = -1
#endif
  end if
  
  do lb = 1, blockCount
     blockID=blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     
     
     nxzones = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     nyzones = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     nzzones = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     
     
#ifndef FIXEDBLOCKSIZE
     allocate(grav(max(nxzones, nyzones, nzzones)))
#endif
     call Grid_getBlkPtr(blockID, solnVec)
     
     if (component == IAXIS) then                    ! x-component
        
        do k = 1, nzzones
           do j = 1, nyzones
              
              pos(1)=j;pos(2)=k
              call Gravity_accelOneRow (pos, SWEEP_X, blockID, nxzones, grav, ptgrav, potVar)
              solnVec(accelIndex,1:nxzones,j,k) = grav(1:nxzones)
              
           enddo
        enddo
        
        
        
     elseif (component == JAXIS) then                ! y-component
        
        do k = 1, nzzones
           do j = 1, nxzones
              
              pos(1)=j; pos(2)=k
              call Gravity_accelOneRow(pos, SWEEP_Y, blockID, nyzones, grav, ptgrav, potVar)
              solnVec(accelIndex,j,1:nyzones,k) = grav(1:nyzones)
              
           enddo
        enddo


        
     elseif (component == KAXIS) then                ! z-component
        
        do k = 1, nyzones
           do j = 1, nxzones
              
              pos(1)=j;pos(2)=k

              call Gravity_accelOneRow (pos, SWEEP_Z, blockID, nzzones, grav, ptgrav, potVar)
              solnVec(accelIndex,j,k,1:nzzones) = grav(1:nzzones)
              
           enddo
        enddo


  
     else                                        ! all components
      
        call Driver_abortFlash &
             ('[Gravity_accelListOfBlocks] ALLAXIS not supported!')
        
     endif
     
     
     
#ifndef FIXEDBLOCKSIZE
     deallocate(grav)
#endif
     call Grid_releaseBlkPtr(blockID,solnVec)
     
  enddo
  
  return
end subroutine Gravity_accelListOfBlocks


