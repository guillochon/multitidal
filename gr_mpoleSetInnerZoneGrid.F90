!!****if* source/Grid/GridSolvers/Multipole_experimental/gr_mpoleSetInnerZoneGrid
!!
!! NAME
!!
!!  gr_mpoleSetInnerZoneGrid
!!
!! 
!! SYNOPSIS
!!
!!  gr_mpoleSetInnerZoneGrid (integer, intent(in)    :: nR_local,
!!                            integer, intent(in)    :: nR_inner_zone,
!!                            integer, intent(in)    :: nP_inner_zone,
!!                            real,    intent(inout) :: R_inner_zone)
!!
!!
!! DESCRIPTION
!!
!!  This routine sets up the inner zone radial grid from all local inner zone
!!  radii on each processor. When exiting this routine, all processors will have
!!  a copy of the inner zone grid and its associated data.
!!
!! ARGUMENTS
!!
!!  nR_local      : number of inner zone radii on the local processor
!!  nR_inner_zone : total number of inner zone radii
!!  nP_inner_zone : total number of processors containing the inner zone
!!  R_inner_zone  : the collection of all inner zone radii
!!
!!***

subroutine gr_mpoleSetInnerZoneGrid (nR_local,      &
                                     nR_inner_zone, &
                                     nP_inner_zone, &
                                     R_inner_zone   )

  use Driver_interface,  ONLY : Driver_abortFlash
  use Grid_data,         ONLY : gr_meshMe, gr_meshComm
  use gr_mpoleData,      ONLY : dr_inner_zone,            &
                                dr_inner_zone_inv,        &
                                inner_zone_rmax,          &
                                inner_zone_qmax,          &
                                inner_zone_radii,         &
                                inner_zone_Qlower,        &
                                inner_zone_Qupper,        &
                                inner_zone_size,          &
                                inner_zone_grid,          &
                                inner_zone_grid_inv

  implicit none

#include "Flash.h"
#include "constants.h"

  include "Flash_mpi.h"

  integer, intent (in)    :: nR_local
  integer, intent (in)    :: nR_inner_zone
  integer, intent (in)    :: nP_inner_zone
  real,    intent (inout) :: R_inner_zone (1:nR_inner_zone)

  logical :: invoke_send

  integer :: dr_unit
  integer :: error
  integer :: n
  integer :: message_tag
  integer :: n_mpi_recv_calls
  integer :: nR_grid
  integer :: nR_received
  integer :: Q,Qstart
  integer :: used_space,free_space

  integer :: status (MPI_STATUS_SIZE)

  real    :: rmax,rmax_prev
!
!
!     ...Those processors which contain inner zone radii will send them to the
!        master processor, which collects them into the big array. After all inner
!        zone radii were collected on the master processor, he orders them and
!        sets up the inner zone radial grid. The inner zone radial grid is then
!        broadcast to all processors.
!
!        These are the individual steps taken:
!
!              1) Local processors -> send local inner zone radii
!              2) Master processor -> collect all local inner zone radii
!                                     into global inner zone radii array
!              3) Master processor -> express the true value global inner
!                                     zone radii in inner zone atomic
!                                     distance units
!              4) Master processor -> order global inner zone radii into
!                                     increasing order
!              5) Master processor -> determine # of inner zone grid radii
!              6) Master processor -> allocate inner zone grid
!              7) Master processor -> calculate inner zone grid radii
!              8) Master processor -> deallocate global inner zone radii array
!              9) All   processors -> broadcast # of inner zone grid radii
!             10) Local processors -> deallocate local inner zone radii array
!             11) Local processors -> allocate inner zone grid
!             12) All   processors -> broadcast inner zone grid
!
!
  message_tag = 1

  invoke_send = (nR_local > 0) .and. (gr_meshMe /= MASTER_PE)

  if (invoke_send) then

      call MPI_Send  (R_inner_zone,        &
                      nR_local,            &
                      FLASH_REAL,          &
                      MASTER_PE,           &
                      message_tag,         &
                      gr_meshComm,      &
                      error                )
  end if

  if (gr_meshMe == MASTER_PE) then

      used_space = nR_local
      free_space = nR_inner_zone - nR_local

      if (nR_local > 0) then
          n_mpi_recv_calls = nP_inner_zone - 1
      else
          n_mpi_recv_calls = nP_inner_zone
      end if

      do n = 1,n_mpi_recv_calls

         call MPI_Recv (R_inner_zone (used_space+1), &
                        free_space,                  &
                        FLASH_REAL,                  &
                        MPI_ANY_SOURCE,              &
                        MPI_ANY_TAG,                 &
                        gr_meshComm,              &
                        status,                      &
                        error                        )

         call MPI_get_Count (status,      &
                             FLASH_REAL,  &
                             nR_received, &
                             error        )

         used_space = used_space + nR_received
         free_space = free_space - nR_received
      end do 

      R_inner_zone = dr_inner_zone_inv * R_inner_zone

      call gr_mpoleHeapsort (nR_inner_zone,R_inner_zone)

      rmax_prev = ceiling (R_inner_zone (1) * inner_zone_grid_inv) * inner_zone_grid

      nR_grid = 1
      do n = 2,nR_inner_zone
         rmax = ceiling (R_inner_zone (n) * inner_zone_grid_inv) * inner_zone_grid
         if (rmax > rmax_prev) then
             nR_grid = nR_grid + 1
             rmax_prev = rmax
         end if
      end do

      if (nR_grid == 0) then
          call Driver_abortFlash ('[gr_mpoleRad3Dcartesian] ERROR: no inner zone grid radii found')
      end if

      allocate (inner_zone_radii (1:nR_grid))

      inner_zone_radii (1) = ceiling (R_inner_zone (1) * inner_zone_grid_inv) * inner_zone_grid

      nR_grid = 1
      do n = 2,nR_inner_zone
         rmax = ceiling (R_inner_zone (n) * inner_zone_grid_inv) * inner_zone_grid
         if (rmax > inner_zone_radii (nR_grid)) then
             nR_grid = nR_grid + 1
             inner_zone_radii (nR_grid) = rmax
         end if
      end do

      inner_zone_qmax = nR_grid

  end if

  call MPI_Bcast (inner_zone_qmax, &
                  1,               &
                  FLASH_INTEGER,   &
                  MASTER_PE,       &
                  gr_meshComm,  &
                  error            )

  if (gr_meshMe /= MASTER_PE) then
      allocate   (inner_zone_radii (1:inner_zone_qmax))
  end if

  call MPI_Bcast (inner_zone_radii, &
                  inner_zone_qmax,  &
                  FLASH_REAL,       &
                  MASTER_PE,        &
                  gr_meshComm,   &
                  error             )
!
!
!     ...Create the inner zone radii locator arrays to help search through
!        the inner zone grid. This is done on all processors.
!
! 
   allocate (inner_zone_Qlower (0:inner_zone_size))   ! lower index 0
 
   inner_zone_Qlower = 0
 
   Qstart = 1
   do n = 0,inner_zone_size                           ! start from 0
      do Q = Qstart,inner_zone_qmax
         dr_unit = int (inner_zone_radii (Q))
         if (dr_unit > n) then
             exit
         else if (dr_unit == n) then
             inner_zone_Qlower (n) = Q
             Qstart = Q + 1
             exit
         end if
      end do
   end do
 
   allocate (inner_zone_Qupper (0:inner_zone_size))   ! lower index 0
 
   inner_zone_Qupper = 0
 
   Qstart = inner_zone_qmax
   do n = inner_zone_size,0,-1                        ! end at 0
      do Q = Qstart,1,-1
         dr_unit = int (inner_zone_radii (Q))
         if (dr_unit < n) then
             exit
         else if (dr_unit == n) then
             inner_zone_Qupper (n) = Q
             Qstart = Q - 1
             exit
         end if
      end do
   end do
!
!
!       ...Ready!
!
!
  return
end subroutine gr_mpoleSetInnerZoneGrid
