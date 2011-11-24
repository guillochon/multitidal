!!****f* source/physics/sourceTerms/Cool/Cool
!!
!! NAME
!!
!!  Cool
!!
!! SYNOPSIS
!!
!!  Cool(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt,
!!          real(IN) :: time)
!!
!!
!!
!! DESCRIPTION
!!  Apply a cooling operator on the list of blocks provided as input
!!
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!!  time : the current time
!!
!!***



subroutine Cool(blockCount,blockList,dt, time)
  use Simulation_data, ONLY : sim_coolingDensity, sim_fluffDampCoeff
  use Grid_data, ONLY: gr_smalle
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  integer :: lb, i, j, k
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(:,:,:,:),pointer :: solnData
  real,intent(IN) :: dt, time

  if (sim_coolingDensity .gt. 0.d0) then
      do lb = 1, blockCount
          call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
          call Grid_getBlkPtr(blockList(lb),solnData)
          do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
              do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                  do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                      if (solnData(DENS_VAR,i,j,k) .lt. sim_coolingDensity) then
                          solnData(EINT_VAR,i,j,k) = max(sim_fluffDampCoeff*solnData(EINT_VAR,i,j,k), gr_smalle)

                          solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                              0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                      endif
                  enddo
              enddo
          enddo
          call Grid_releaseBlkPtr(blockList(lb), solnData)
      enddo
  endif

  return
end subroutine Cool
