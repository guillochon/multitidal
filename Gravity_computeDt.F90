!!****f* source/physics/Gravity/Gravity_computeDt
!!
!! NAME
!!
!!  Gravity_computeDt
!!  
!! SYNOPSIS
!!
!!  Gravity_computeDt(integer(IN)        :: blockID,
!!                    pointer            :: solnData,
!!                    real (OUT)         :: dt_grav,
!!                    integer(:)(INOUT)  :: dt_minloc(5))
!!
!! DESCRIPTION
!!
!!  Compute the timestep limiter due to the gravitational solver.
!!
!! ARGUMENTS
!!
!!  dt_grav:       Will Return the limiting timestep. Should be
!!                 set to a large value (1.D99) on input.
!!  dt_minloc(5):  An array to receive information about which
!!                 processor, block, and zone was responsible
!!                 for setting the limiting timestep.  The order
!!                 is i, j, k, b, p, where (i,j,k) = zone
!!                 indices, b = local block ID, and p = PE #.
!!                 This routine should only modify these values
!!                 if it changes dt_grav.
!!  blockID:       The local ID of the block to compute the
!!                 limiter on.
!!
!!***

subroutine Gravity_computeDt (blockID, solnData, dt_grav, dt_minloc)

!==============================================================================

#include "Flash.h"
#include "constants.h"

  use Simulation_data, ONLY : sim_gCell, grv_cfl, sim_comAccel
  use Gravity_data, ONLY : grv_meshMe
  implicit none
  
  integer, intent(IN)    ::  blockID
  real, pointer :: solnData(:,:,:,:) 
  integer, intent(INOUT) ::  dt_minloc(5)
  real,intent(OUT)       ::  dt_grav

  integer :: sizeX,sizeY,sizeZ,i,j,k,temploc(5)
  real :: delxinv, delyinv, delzinv, dt_ltemp, dt_temp
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(MDIM) :: del

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  call Grid_getCellCoords(IAXIS,blockId,CENTER,sim_gCell,xCoord,sizeX)
  call Grid_getCellCoords(JAXIS,blockId,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(KAXIS,blockId,CENTER,sim_gCell,zCoord,sizeZ)
  call Grid_getDeltas(blockID,del)

  ! Only CARTESIAN works for now
  delyinv = 1.
  delzinv = 1.

  delxinv = 1.0/del(1)
  if (NDIM > 1) &
  delyinv = 1.0/del(2)
  if (NDIM > 2) &
  delzinv = 1.0/del(3)

  dt_temp    = 0.
  temploc(:) = 0

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           dt_ltemp = sqrt(abs(solnData(SGAX_VAR,i,j,k)-sim_comAccel(1))*delxinv)
           if (NDIM > 1) dt_ltemp = max(dt_ltemp,sqrt(abs(solnData(SGAY_VAR,i,j,k)-sim_comAccel(2))*delyinv))
           if (NDIM > 2) dt_ltemp = max(dt_ltemp,sqrt(abs(solnData(SGAZ_VAR,i,j,k)-sim_comAccel(3))*delzinv))

           if (dt_ltemp > dt_temp) then
              dt_temp    = dt_ltemp
              temploc(1) = i
              temploc(2) = j
              temploc(3) = k
              temploc(4) = blockID
              temploc(5) = grv_meshMe
           endif
        enddo
     enddo
  enddo

  dt_temp = grv_cfl / dt_temp
  if (dt_temp < dt_grav) then
     dt_grav = dt_temp
     dt_minloc = temploc
  endif
  
  if(dt_grav <= 0.0) then
     print*,dt_grav
     call Driver_abortFlash("[Gravity]: Computed dt is not positive! Aborting!")
  endif

  return
end subroutine Gravity_computeDt
