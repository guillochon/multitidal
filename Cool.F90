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

  use Grid_interface, ONLY : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
  use Multispecies_interface, ONLY : Multispecies_getSumInv
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Simulation_data, ONLY : sim_smallT

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"

  integer, intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt, time

  integer :: thisBlock, blockID, i, j, k
  logical :: cooledZone
  real :: sdot, sgamma, ek, ei, rho, mp, abar
  real, pointer, dimension(:,:,:,:)            :: solnData
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  double precision, dimension(SPECIES_BEGIN:SPECIES_END) :: xn

  sgamma = 3.e-22 !From Burkert 2012

  call PhysicalConstants_get("proton mass", mp)

  ! make sure that guardcells are up to date
  call Grid_fillGuardCells(CENTER, ALLDIR)

  ! loop over list of blocks passed in
  do thisBlock = 1, blockCount

     blockID = blockList(thisBlock)
     cooledZone = .FALSE.

     ! get dimensions/limits and coordinates
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Get a pointer to solution data 
     call Grid_getBlkPtr(blockID,solnData)

     ! now guaranteed that tmp, rho, etc. exist
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              cooledZone = .true.

              rho  = solnData(DENS_VAR,i,j,k)

              xn = solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)

              call Multispecies_getSumInv(A, abar, xn)
              abar = 1.d0 / abar

              sdot = sgamma*rho/(abar*mp)**2

              ! kinetic energy
              ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  & 
                   solnData(VELY_VAR,i,j,k)**2 +  & 
                   solnData(VELZ_VAR,i,j,k)**2)

              ! internal energy, add on nuclear rate*timestep
              ei = solnData(ENER_VAR,i,j,k) - ek
              ei = ei - dt*sdot
                
#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k) = ei
#endif
              solnData(ENER_VAR,i,j,k) = ei + ek
              solnData(ECOO_VAR,i,j,k) = sdot
           enddo
        enddo
     enddo

     ! we've altered the EI, let's equilabrate
     if (cooledZone) then
        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 solnData(TEMP_VAR,i,j,k) = max(solnData(TEMP_VAR,i,j,k), sim_smallT)
              enddo
           enddo
        enddo
        call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blockID)
     end if

     call Grid_releaseBlkPtr(blockID,solnData)

  end do

  return
end subroutine Cool
