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

subroutine Gravity_computeDt (blockID, &
                              dx, dy, dz, blkLimits, blkLimitsGC, &
                              solnData, dt_grav, dt_minloc)

!==============================================================================

#include "Flash.h"
#include "constants.h"

  use Simulation_data, ONLY : sim_gCell, grv_cfl, sim_comAccel
  use Gravity_data, ONLY : grv_meshMe
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getDeltas
  implicit none
  
  integer, intent(IN)    ::  blockID
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
  real, pointer :: solnData(:,:,:,:) 
  integer, intent(INOUT) ::  dt_minloc(5)
  real,intent(OUT)       ::  dt_grav

  integer :: i,j,k,temploc(5)
  real :: delxinv, delyinv, delzinv, dt_ltemp, dt_temp
  real, dimension(MDIM) :: del

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: dx
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: dy
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: dz
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: dx
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: dy
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: dz
#endif

  ! Only CARTESIAN works for now
  delyinv = 1.
  delzinv = 1.

  delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
  if (NDIM > 1) &
  delyinv = 1.0/dy(blkLimits(LOW,JAXIS))
  if (NDIM > 2) &
  delzinv = 1.0/dz(blkLimits(LOW,KAXIS))

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
