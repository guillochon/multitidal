!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleMoments
!!
!! NAME
!!
!!  gr_mpoleMoments
!!
!! 
!! SYNOPSIS
!!
!!  gr_mpoleMoments(integer, intent(in) : idensvar)
!!
!!
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution,
!!  assuming the center of mass and the total mass have first been
!!  computed by find_center_of_mass().  On output, the Moment()
!!  array declared in mpole_common.F90 contains the moments.
!!  this subroutine does nothing if the geometry is 2-d spherical
!!
!! PARAMETERS
!!
!!  mpole_subSample  -- integer to control number of subzones in each direction
!!                           for smoothing potential calculations.
!!                      Also slows down calculation considerably when != 1
!!
!!***

subroutine gr_mpoleMoments (idensvar)

  use gr_mpoleData, ONLY : mpole_geometry,&
                         G_2DSPHERICAL,G_1DSPHERICAL,&
                         G_2DCYLINDRICAL,G_3DCARTESIAN,G_3DAXISYMMETRIC,&
                         Moment,OldMoment,Outer,Momtmp,cartfactor,&
                         leg_fact,Inner,qmax,Even,Odd,mpole_lmax,mpole_mmax,&
                         fourpi,cylfactor,Xcm,Ycm,Zcm, mpole_dumpMoments, &
                         mpole_useMatrixMPI, MomtmpMatrix, mpole_subSample
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
                             Grid_getDeltas,Grid_getBlkBoundBox,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,&
                             Grid_getBlkIndexLimits
  use Timers_interface, ONLY:  Timers_start, Timers_stop
  
  implicit none
  
#include "constants.h"
#include "Flash.h"
 include "Flash_mpi.h"


 integer, intent(in)   :: idensvar
 
 integer   :: i, j, k, l, m, error
 real      :: zonemass, dvol, delx, dely, delz, x, y, z
 
 integer   :: lb, ii, jj, kk 
 !!  integer   :: Nint = 2 !! Now replaced with RP mpole_subSample
 real      :: xx, yy, zz
 real      :: delxx, delyy, delzz
 real      :: zonedens
 
 integer :: q, imax, jmax, kmax, imin, jmin, kmin
 integer :: compressMax
 
 integer :: blockCount
 integer :: blockList(MAXBLOCKS)
 
 integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
 real,dimension(LOW:HIGH,MDIM) :: bndBox
 real,dimension(MDIM) :: delta
 real, pointer,dimension(:,:,:,:)      :: solnData
 
 
 ! this part is not used with current implementation in 2-D spherical geometry
 
 if ( mpole_geometry /= G_2DSPHERICAL ) then
    
    
    !Sum quantities over all locally held leaf blocks.
    
    OldMoment = Moment
    Moment(:,:,:,:,:) = 0.e0
    
    call Grid_getListOfBlocks(LEAF, blockList, blockCount)
    
    do lb = 1, blockCount
       
        
       call Grid_getDeltas(blockList(lb),delta)
       call Grid_getBlkBoundBox(blockList(lb),bndBox)
       !handle cases if less than 2d or 3d
       delx = delta(IAXIS)
       dely = delta(JAXIS)
       delz = delta(KAXIS)
       
       delxx = delx / real(mpole_subSample)  !half of del(IAXIS)
       delyy = dely / real(mpole_subSample)
       delzz = delz / real(mpole_subSample)
       
       !Get pointer to solution data.
       call Grid_getBlkPtr(blockList(lb), solnData)
       
       !get limits of blk.. this is compatible with blk sizes not fixed at compile time
       call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
       
       kmax = blkLimits(HIGH,KAXIS)
       kmin = blkLimits(LOW,KAXIS)  
       jmax = blkLimits(HIGH,JAXIS)
       jmin = blkLimits(LOW,JAXIS)
       imax = blkLimits(HIGH,IAXIS)
       imin = blkLimits(LOW,IAXIS)
       
       !Compute the moment contributions for this block.
       
       do k = kmin, kmax
          z = (bndBox(LOW,KAXIS) + (k-kmin)*delz - Zcm) * K3D
          
          do j = jmin, jmax
             y = (bndBox(LOW,JAXIS) + (j-jmin)*dely - Ycm) * K2D
             
             do i = imin, imax
                x = bndBox(LOW,IAXIS) + (i-imin)*delx - Xcm
                
                
     ! Note Carlo thinks this is suspicious that the zonedens depends on i,j,k
     !        rather than ii,jj,kk.  In fact it is correct, as solnData doesn't exist at ii,jj,kk
                zonedens = solnData(idensvar,i,j,k)


     ! Subsampling here occurs over the legendre polynomials/cos/sin arguments 
                do kk = 1, 1+(mpole_subSample-1)*K3D
                   zz = (z + (kk-0.5)*delzz) * K3D
                   
                   do jj = 1, 1+(mpole_subSample-1)*K2D
                      yy = (y + (jj-0.5)*delyy) * K2D
                      
                      do ii = 1, mpole_subSample
                         xx = x + (ii-0.5)*delxx

                         
                         select case (mpole_geometry)
                            
                         case (G_3DCARTESIAN)
                            ! cartfactor includes the factor of 8 implied by octant symmetry
                            ! if we are doing an octant.  The point is to treat each zone like
                            ! eight, so the total mass is computed as the mass of the full star,
                            ! not just 1 octant.
                            dvol = cartfactor * delxx * delyy * delzz
                            zonemass = zonedens * dvol
                            call gr_zoneMoments (xx, yy, zz, zonemass)

                         case (G_3DAXISYMMETRIC)
                            dvol = cartfactor * delxx * delyy * delzz
                            zonemass = zonedens * dvol
                            call gr_zoneMoments (xx, yy, zz, zonemass)
                            
                         case (G_2DCYLINDRICAL)
                            dvol = cylfactor * xx * delxx * delyy
                            zonemass = zonedens * dvol
                            call gr_zoneMoments (xx, 0., yy, zonemass)
                            
                         case (G_1DSPHERICAL)
                            dvol = fourpi/3. * ((xx+0.5*delxx)**3 - (xx-0.5*delxx)**3)
                            !  Note that the 1-D case DOES do interpolation of zonedens.  The higher-
                            !  dimensional cases do not.
                            if (i < imax) then
                               zonedens = solnData(idensvar,i,j,k) + &
                                    (solnData(idensvar,i+1,j,k) - &
                                    solnData(idensvar,i,j,k))/delx * &
                                    (xx - (x+0.5*delx))
                            else
                               zonedens = solnData(idensvar,i-1,j,k) + &
                                    (solnData(idensvar,i,j,k) - &
                                    solnData(idensvar,i-1,j,k))/delx * &
                                    (xx - (x-0.5*delx))
                            endif
                            zonemass = zonedens * dvol
                            call gr_zoneMoments (xx, 0., 0., zonemass)
                            
                         end select
                         
                      enddo
                   enddo
                enddo
                
             enddo
          enddo
       enddo
       
       call Grid_releaseBlkPtr(blockList(lb), solnData)
       
    enddo
    
    ! In zone_moments, we only added the contribution of each zone to its
    ! particular radial bin.  Each zone should contribute its inner moment to
    ! all zones with radius greater than it, and its outer moment to all zones
    ! with radius less than it.
    
    do m = 0, mpole_mmax
       do l = m, mpole_lmax
          do i = Even, Odd
             do q = 2, qmax
                Moment(q,i,Inner,l,m) = Moment(q,i,Inner,l,m) + Moment(q-1,i,Inner,l,m)
             enddo
             do q = qmax-1, 1, -1
                Moment(q,i,Outer,l,m) = Moment(q,i,Outer,l,m) + Moment(q+1,i,Outer,l,m)
             enddo
          enddo
       enddo
    enddo
    
    !==============================================================================
    
    !               Give all processors a copy of the global sums.
    if (.not. mpole_useMatrixMPI) then
        call Timers_start ("Looped Moment All Reduce")
      do m = 0, mpole_mmax
          do l = m, mpole_lmax
             do j = Inner, Outer
                do i = Even, Odd
                   call mpi_allreduce (Moment(0,i,j,l,m), Momtmp(0), qmax+1, & 
                        FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, &
                        error)
                   Moment(:,i,j,l,m) = Momtmp(:)
                enddo
             enddo
          enddo
       enddo
        call Timers_stop ("Looped Moment All Reduce")

    else
       ! Anshu and Carlo believe that the collapsed sums are valid only with mpole_mmax=0
       ! LBR's version without the loops

        call Timers_start ("Matrix Moment All Reduce")

       compressMax = (mpole_mmax+1)*(mpole_lmax+1)*2*2*(qmax+1)
       call mpi_allreduce (Moment, MomtmpMatrix, compressMax, & 
            FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, &
            error)
       Moment = MomtmpMatrix

        call Timers_stop ("Matrix Moment All Reduce")
    end if
    
    !               Normalize the moments properly.
    
    if ((mpole_lmax > 0) .and. (mpole_mmax > 0)) then
       do l = 1, mpole_lmax
          do m = 1, l
             Moment(:,:,:,l,m) = Moment(:,:,:,l,m) * Leg_fact(l,m)
          enddo
       enddo
    endif
    
 end if
 !==============================================================================

 !  dump data if requested by runtime parameter mpole_dumpMoments
 if (mpole_dumpMoments) call gr_mpoleDumpMoments()
 
 return
end subroutine gr_mpoleMoments
