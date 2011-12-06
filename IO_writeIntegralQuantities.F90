!!****if* source/IO/IOMain/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
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
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities (isFirst, simTime)

    use IO_data, ONLY : io_restart, io_statsFileName, io_globalComm
    use IO_data, ONLY : io_globalMe
    use Grid_interface, ONLY : Grid_getListOfBlocks, &
      Grid_getBlkIndexLimits, Grid_getBlkPtr, &
      Grid_releaseBlkPtr, Grid_getDeltas, Grid_getBlkBoundBox
    use Gravity_data, ONLY : grv_ptvec, grv_obvec, grv_ptmass, grv_boundvec, &
        grv_exactvec, grv_totmass, grv_thresh, grv_momacc, grv_angmomacc, grv_eneracc, grv_massacc
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get

    implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
    
    real, intent(in) :: simTime

    integer, intent(in) :: isFirst

    integer :: lb, count
    
    integer :: funit = 99
    integer :: error
    
    integer :: blockList(MAXBLOCKS)

    integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

    integer, parameter ::  nGlobalSum = 21         ! Number of globally-summed quantities

    double precision, dimension(nGlobalSum) :: gsum, lsum

    integer :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, sizeX, sizeY, sizeZ, istat
    real :: dvol, xx, yy, zz, x, y, z, v2, velx, vely, velz, ptrad, ptv2, newton
    real :: ptx, pty, ptz, ptvx, ptvy, ptvz
    real, DIMENSION(:,:,:,:), POINTER :: solnData
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    real,dimension(MDIM) :: delta
    real,dimension(LOW:HIGH,MDIM) :: bndBox
    logical :: gcell = .true.

    ! Sum quantities over all locally held leaf-node blocks.
    gsum = 0.
    lsum = 0.
    
    call Grid_getListOfBlocks(LEAF, blockList, count)
    call PhysicalConstants_get("Newton", newton)
    
    do lb = 1, count
       call Grid_getDeltas(blockList(lb),delta)
       call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
       call Grid_getBlkBoundBox(blockList(lb), bndBox)
       sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
       allocate(xCoord(sizeX),stat=istat)
       sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
       allocate(yCoord(sizeY),stat=istat)
       sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
       allocate(zCoord(sizeZ),stat=istat)

       if (NDIM == 3) call Grid_getCellCoords&
                           (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
       if (NDIM >= 2) call Grid_getCellCoords&
                           (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
       call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
       call Grid_getBlkPtr(blockList(lb),solnData)
       dvol = delta(IAXIS) * delta(JAXIS) * delta(KAXIS)
       kmax = blkLimits(HIGH,KAXIS)
       kmin = blkLimits(LOW,KAXIS)  
       jmax = blkLimits(HIGH,JAXIS)
       jmin = blkLimits(LOW,JAXIS)
       imax = blkLimits(HIGH,IAXIS)
       imin = blkLimits(LOW,IAXIS)

       ! Sum contributions from the indicated blkLimits of cells.
       do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
           if (NDIM == 3) zz = bndBox(LOW,KAXIS) + (k-kmin+0.5)*delta(KAXIS)
           do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
               if (NDIM >= 2) yy = bndBox(LOW,JAXIS) + (j-jmin+0.5)*delta(JAXIS)
               do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                
                ! mass   
#ifdef DENS_VAR
                lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
                ! momentum
                velx = solnData(VELX_VAR,i,j,k)
                lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k)*velx*dvol
#endif
#ifdef VELY_VAR      
                vely = solnData(VELY_VAR,i,j,k)
                lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k)*vely*dvol
#endif
#ifdef VELZ_VAR      
                velz = solnData(VELZ_VAR,i,j,k)
                lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k)*velz*dvol
#endif
                ! total energy
#ifdef ENER_VAR
                lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                     &                                solnData(DENS_VAR,i,j,k)*dvol
#endif

             
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
                ! kinetic energy
                v2 = velx**2.d0 + vely**2.d0 + velz**2.d0
                lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k)*dvol*v2
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
                lsum(8) = lsum(8) + solnData(GPOT_VAR,i,j,k) * solnData(DENS_VAR,i,j,k)*dvol
#endif
                xx = bndBox(LOW,IAXIS) + (i-imin+0.5d0)*delta(IAXIS)
                ptx = xx - grv_exactvec(1) - grv_ptvec(1) + grv_obvec(1)
                pty = yy - grv_exactvec(2) - grv_ptvec(2) + grv_obvec(2)
                ptz = zz - grv_exactvec(3) - grv_ptvec(3) + grv_obvec(3)
                ptrad = dsqrt(ptx**2.d0 + pty**2.d0 + ptz**2.d0)
                ptvx = velx + grv_obvec(4) - grv_ptvec(4)
                ptvy = vely + grv_obvec(5) - grv_ptvec(5)
                ptvz = velz + grv_obvec(6) - grv_ptvec(6)
                ptv2 = ptvx**2.d0 + ptvy**2.d0 + ptvz**2.d0
                x = xx - grv_exactvec(1)
                y = yy - grv_exactvec(2)
                z = zz - grv_exactvec(3)
                lsum(9)  = lsum(9)  + solnData(DENS_VAR,i,j,k)*dvol*(y*velz - z*vely)
                lsum(10) = lsum(10) + solnData(DENS_VAR,i,j,k)*dvol*(z*velx - x*velz)
                lsum(11) = lsum(11) + solnData(DENS_VAR,i,j,k)*dvol*(x*vely - y*velx)
#endif ! ifdef DENS_VAR
                ! summed quantities for only bound material
                ! momentum
                velx = solnData(VELX_VAR,i,j,k) - grv_boundvec(4)
                vely = solnData(VELY_VAR,i,j,k) - grv_boundvec(5)
                velz = solnData(VELZ_VAR,i,j,k) - grv_boundvec(6)
                v2 = 0.5d0*(velx**2.d0 + vely**2.d0 + velz**2.d0)

                if (-solnData(GPOT_VAR,i,j,k) .ge. v2) then
                    x = xx - grv_boundvec(1)
                    y = yy - grv_boundvec(2)
                    z = zz - grv_boundvec(3)
                    v2 = 0.5d0*(velx**2.d0 + vely**2.d0 + velz**2.d0)
                    lsum(12) = lsum(12) + 0.5d0*solnData(DENS_VAR,i,j,k) * dvol * v2
                    lsum(13) = lsum(13) + solnData(DENS_VAR,i,j,k) * solnData(EINT_VAR,i,j,k)*dvol
                    lsum(14) = lsum(14) + solnData(GPOT_VAR,i,j,k) * solnData(DENS_VAR,i,j,k)*dvol
                    lsum(15) = lsum(15) + solnData(DENS_VAR,i,j,k)*dvol*(y*velz - z*vely)
                    lsum(16) = lsum(16) + solnData(DENS_VAR,i,j,k)*dvol*(z*velx - x*velz)
                    lsum(17) = lsum(17) + solnData(DENS_VAR,i,j,k)*dvol*(x*vely - y*velx)
                elseif (newton*grv_ptmass/ptrad .gt. ptv2) then
                    lsum(18) = lsum(18) + solnData(DENS_VAR,i,j,k)*dvol
                    lsum(19) = lsum(19) + solnData(DENS_VAR,i,j,k)*dvol*(pty*ptvz - ptz*ptvy)
                    lsum(20) = lsum(20) + solnData(DENS_VAR,i,j,k)*dvol*(ptz*ptvx - ptx*ptvz)
                    lsum(21) = lsum(21) + solnData(DENS_VAR,i,j,k)*dvol*(ptx*ptvy - pty*ptvx)
                endif
             enddo
          enddo
       enddo
       call Grid_releaseBlkPtr(blockList(lb), solnData)

    enddo
    
    call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, MASTER_PE, io_globalComm, error)

    if (io_globalMe == MASTER_PE) then
       
       ! create the file from scratch if it is a not a restart simulation, 
       ! otherwise append to the end of the file
       if (isFirst == 0) then
          open (funit, file=trim(io_statsFileName), position='APPEND')
       else 
          if (.NOT. io_restart) then
             open (funit, file=trim(io_statsFileName)) 
             write (funit, 10)               &
                  '#time                  ', &
                  'mass                   ', &
                  'x-momentum             ', &
                  'y-momentum             ', & 
                  'z-momentum             ', &
                  'E_total                ', &
                  'E_kinetic              ', &
                  'E_internal             ', &
                  'E_gravity              ', &
                  'Ang. mom. tot, x comp. ', &
                  'y comp.                ', &
                  'z comp.                ', &
                  'E_kinetic, bound       ', &
                  'E_internal, bound      ', &
                  'E_gravity, bound       ', &
                  'Ang. mom. bnd., x comp.', &
                  'y comp.                ', &
                  'z comp.                ', &
                  'Mass bound to pt. mass ', &
                  'Ang. mom. bnd. pt., x  ', &
                  'y comp.                ', &
                  'z comp.                ', &
                  'Mass acc. pt.          ', &
                  'Mom. acc. pt., x       ', &
                  'y comp.                ', &
                  'z comp.                ', &
                  'Ang. mom. acc. pt., x  ', &
                  'y comp.                ', &
                  'z comp.                ', &
                  'Ener. acc. pt.         '

10        format (2x,58(a25, :, 1X))

          else
             open (funit, file=trim(io_statsFileName), position='APPEND')
             write (funit, 11) 
11        format('# simulation restarted')
          endif
       endif
       
       write (funit, 12) simtime, gsum, grv_massacc, grv_momacc, grv_angmomacc, grv_eneracc ! Write the global sums to the file.
12     format (1x, 58(es25.15, :, 1x))
   
       close (funit)          ! Close the file.
       
    endif
    
    call MPI_Barrier (io_globalComm, error)
    
    !=============================================================================
    
    return
end subroutine IO_writeIntegralQuantities



