!!****if* source/IO/IOMain/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: myPE, 
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   myPE - current processor number
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities (myPE, isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr
  use gr_mpoleData, ONLY: Xcm, Ycm, Zcm

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
  
  integer, intent(in) :: myPE
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter ::  nGlobalSum = 11         ! Number of globally-summed quantities
  real :: gsum(nGlobalSum), gvel(4) !Global summed quantities
  real :: lsum(nGlobalSum), lvel(4) !Global summed quantities

  integer :: i, j, k
  real :: dvol, ldisp, gdisp
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  gdisp = 0.
  ldisp = 0.
  gvel = 0.
  lvel = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
     
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
              lvel(1) = lvel(1) + solnData(VELX_VAR,i,j,k)
           
#endif
#ifdef VELY_VAR      

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
              lvel(2) = lvel(2) + solnData(VELY_VAR,i,j,k)
           
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
              lvel(3) = lvel(3) + solnData(VELZ_VAR,i,j,k)
#endif
              lvel(4) = lvel(4) + 1.0

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#ifdef MAGP_VAR
              ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * & 
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

              lsum(5) = lsum(5) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif

#ifdef GPOT_VAR
              ! gravitational potential energy
              lsum(8) = lsum(8) + solnData(GPOT_VAR,i,j,k)/2. * solnData(DENS_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR

#ifdef MAGP_VAR
              ! magnetic energy
              lsum(9) = lsum(9) + solnData(MAGP_VAR,i,j,k)*dvol
#endif

#ifdef ENUC_VAR
              ! nuclear energy generation
              lsum(10) = lsum(10) + solnData(DENS_VAR,i,j,k)*solnData(ENUC_VAR,i,j,k)*dvol
#endif

#ifdef ECOO_VAR
              ! neutrino cooling rate
              lsum(11) = lsum(11) + solnData(DENS_VAR,i,j,k)*solnData(ECOO_VAR,i,j,k)*dvol
#endif

           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, MPI_COMM_WORLD, error)
  call MPI_AllReduce (lvel, gvel, 4, FLASH_REAL, MPI_SUM, & 
       &                MPI_COMM_WORLD, error)
  gvel(1:3) = gvel(1:3) / gvel(4)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)

#ifdef DENS_VAR
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! dispersion energy
              ldisp = ldisp + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &   ((solnData(VELX_VAR,i,j,k) - gvel(1))**2+ & 
                   &    (solnData(VELY_VAR,i,j,k) - gvel(2))**2+ & 
                   &    (solnData(VELZ_VAR,i,j,k) - gvel(3))**2)*dvol           

#endif
#endif
#endif
#endif
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo

  call MPI_Reduce (ldisp, gdisp, 1, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, MPI_COMM_WORLD, error)

  if (MyPE == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time             ', &
                'mass              ', &
                'x-momentum        ', &
                'y-momentum        ', & 
                'z-momentum        ', &
                'E_total           ', &
                'E_kinetic         ', &
                'E_internal        ', &
                'E_gravity         ', &
                'E_nuclear         ', &
                'E_cooling         ', &
                'E_magnetic        ', &
                'E_dispersion      ', &
                'Xcm               ', &
                'Ycm               ', &
                'Zcm               '

10         format (2x,50(a20, :, 1X))

        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12) simtime, gsum, gdisp, Xcm, Ycm, Zcm     ! Write the global sums to the file.
12   format (1x, 50(es20.10, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



