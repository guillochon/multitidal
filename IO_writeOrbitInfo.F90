!!****if* source/IO/IOMain/IO_writeOrbitInfo
!!
!!
!!  NAME
!!    IO_writeOrbitInfo
!!
!!  SYNOPSIS
!!    call IO_writeOrbitInfo(integer(in) :: isFirst,
!!                           real(in)    :: simTime)
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

subroutine IO_writeOrbitInfo (isFirst, simTime)

    use IO_data, ONLY : io_restart, io_statsFileName
    use Grid_interface, ONLY : Grid_getListOfBlocks, &
      Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
      Grid_releaseBlkPtr
    use gr_isoMpoleData, ONLY: Xcm, Ycm, Zcm
    use Gravity_data, ONLY: grv_ptvec, grv_obvec, grv_ptmass, grv_bound, &
        grv_ener, grv_tot_ener, grv_boundvec, grv_exactvec, grv_peakvec, &
        grv_mpolevec
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use Grid_data, ONLY: gr_meshMe

    implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
    
    real, intent(in) :: simTime

    integer, intent(in) :: isFirst

    integer :: lb, error, blockCount
    
    integer :: funit = 95
    
    character (len=MAX_STRING_LENGTH), save :: fname 
    
    integer :: blockList(MAXBLOCKS), blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

    integer :: i, j, k
    real, DIMENSION(:,:,:,:), POINTER :: solnData

    integer :: point(MDIM)
    double precision :: r, vel, arglat, semimaj, hsq, ecc, radvel, cose, sine, &
      truanom, longperi, arguperi, newton, tinitial, raasc, inclin, &
      h1, h2, h3, p1, p2, x, y, z, vx, vy, vz
    character(len=50) :: filename

    call PhysicalConstants_get("Newton", newton)
    call Grid_getListOfBlocks(LEAF,blockList,blockCount)
    call Orbit_energy(blockCount, blockList)

    write(filename, fmt=13)

13  format("orbit.dat")

    if (gr_meshMe == MASTER_PE) then
       
       ! create the file from scratch if it is a not a restart simulation, 
       ! otherwise append to the end of the file
       if (isFirst == 0) then
          open (funit, file=trim(filename), position='APPEND')
       else 
          if (.NOT. io_restart) then
             open (funit, file=trim(filename)) 
             write (funit, 10)               &
                  '#Time                  ', &
                  'Obj. 1: X              ', &
                  'Y                      ', &
                  'Z                      ', &
                  'V_X                    ', & 
                  'V_Y                    ', &
                  'V_Z                    ', &
                  'Obj. 2: X              ', &
                  'Y                      ', &
                  'Z                      ', &
                  'V_X                    ', & 
                  'V_Y                    ', &
                  'V_Z                    ', &
                  'Self-bound: Xcm        ', &
                  'Ycm                    ', &
                  'Zcm                    ', &
                  'V_Xcm                  ', & 
                  'V_Ycm                  ', &
                  'V_Zcm                  ', &
                  'Tot. box: Xcm          ', &
                  'Ycm                    ', &
                  'Zcm                    ', &
                  'V_Xcm                  ', & 
                  'V_Ycm                  ', &
                  'V_Zcm                  ', &
                  'Peak density: X        ', &
                  'Y                      ', &
                  'Z                      ', &
                  'V_X                    ', & 
                  'V_Y                    ', &
                  'V_Z                    ', &
                  'Mpole Exp.: Xcm        ', &
                  'Ycm                    ', &
                  'Zcm                    ', &
                  'V_Xcm                  ', & 
                  'V_Ycm                  ', &
                  'V_Zcm                  ', &
                  'Radius                 ', &
                  'Velocity               ', &
                  'Arg. of latitude       ', &
                  'Semi-major axis        ', &
                  'Eccentricity           ', &
                  'Radial Velocity        ', &
                  'True Anomaly           ', &
                  'Arg. of Peri.          ', &
                  'R.A. Asc. Node         ', &
                  'Long. of Peri.         ', &
                  'Inclination            ', &
                  'Bound mass             ', &
                  'Orbital energy         ', &
                  'Orb. energy, bound     ', &
                  'Mass of point mass     '

10         format (2x,52(a25, :, 1X))

          else
             open (funit, file=trim(filename), position='APPEND')
             write (funit, 11) 
11         format('# simulation restarted')
          endif
       endif
       
       x  = grv_obvec(1) - grv_ptvec(1) + grv_boundvec(1) - grv_exactvec(1)
       y  = grv_obvec(2) - grv_ptvec(2) + grv_boundvec(2) - grv_exactvec(2)
       z  = grv_obvec(3) - grv_ptvec(3) + grv_boundvec(3) - grv_exactvec(3)
       vx = grv_obvec(4) - grv_ptvec(4) + grv_boundvec(4) - grv_exactvec(4)
       vy = grv_obvec(5) - grv_ptvec(5) + grv_boundvec(5) - grv_exactvec(5)
       vz = grv_obvec(6) - grv_ptvec(6) + grv_boundvec(6) - grv_exactvec(6)
       r = dsqrt(x**2.d0 + y**2.d0 + z**2.d0)
       h1 = y*vz - z*vy
       h2 = z*vx - x*vz
       h3 = x*vy - y*vx
       hsq = h1**2.d0 + h2**2.d0 + h3**2.d0
       raasc = datan2(h1, -h2)
       inclin = datan2(dsqrt(h1**2.d0 + h2**2.d0), h3)
       vel = dsqrt(vx**2.d0 + vy**2.d0 + vz**2.d0)
       p1 = x*dcos(raasc) - y*dsin(raasc)
       p2 = x*dcos(inclin)*dsin(raasc) + y*dcos(inclin)*dcos(raasc) - z*dsin(inclin)
       arglat = datan2(p2, p1)
       semimaj = newton*(grv_ptmass + grv_bound)*r/(2.d0*newton*(grv_ptmass + grv_bound) - r*vel**2.d0)
       ecc = dsqrt(1.d0 - hsq/newton/(grv_ptmass + grv_bound)/semimaj)
       radvel = (x*vx + y*vy + z*vz)/r
       cose = (semimaj - r)/semimaj/ecc
       sine = r*radvel/ecc/dsqrt(newton*(grv_ptmass+grv_bound)*semimaj)
       truanom = datan2(dsqrt(1.d0 - ecc**2.d0)*sine, cose - ecc)
       arguperi = arglat - truanom
       longperi = dmod(arguperi - raasc, 2.d0*PI)

       write (funit, 12) simTime, grv_ptvec, grv_obvec, grv_boundvec, grv_exactvec, grv_peakvec, grv_mpolevec, r, vel, &
           arglat, semimaj, ecc, radvel, truanom, arguperi, raasc, longperi, inclin, grv_bound, grv_tot_ener, grv_ener, grv_ptmass
12     format (1x, 52(es25.15, :, 1x))
   
       close (funit)          ! Close the file.
       
    endif
    
    call MPI_Barrier (MPI_Comm_World, error)
    
    !=============================================================================
    
    return
end subroutine IO_writeOrbitInfo
