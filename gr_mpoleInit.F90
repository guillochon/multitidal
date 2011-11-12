!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleInit
!!
!! NAME
!!
!!  gr_mpoleInit
!!
!! 
!! SYNOPSIS
!!
!!  gr_mpoleInit()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!! PARAMETERS
!!
!!  mpole_subSample  -- integer to control number of subzones in each direction
!!                           for smoothing potential calculations.
!!                      Also slows down calculation considerably when != 1
!!
!!
!!***

subroutine gr_mpoleInit()

  use Grid_data, ONLY : gr_myPE, gr_geometry
  use gr_mpoleData

  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_stampMessage
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get

  implicit none

#include "Flash.h"
#include "constants.h"

  real    :: Nxmax, Nymax, Nzmax
  real    :: factrl
  integer :: istat, m, l, qmid
  real    :: cxzn, cxzn_1, ds, maxrad



  character(len=256) :: str_buffer
  character(len=32)  :: int_to_str


  twopi  = 2.*PI
  fourpi = 4.*PI

  call PhysicalConstants_get("Newton", newton)

  call RuntimeParameters_get("mpole_lmax", mpole_lmax)

  call RuntimeParameters_get("quadrant", quadrant)
  call RuntimeParameters_get("octant",   octant)
  call RuntimeParameters_get("mpole_3daxisymmetric",  axisym)
  call RuntimeParameters_get("mpole_dumpMoments", mpole_dumpMoments)

  !! True if you want to collapse the looped mpi_allreduce calls
  call RuntimeParameters_get("mpole_useMatrixMPI", mpole_useMatrixMPI)

  !! Integer to control subsampling in potential and moment calculations.  Set to 1 for no subsampling
  ! This subsampling parameter was originally set to Nint6=6 in the potential calculations, and 
  !   Nint=2 in the moment calculations.  It's value is now controlled by a runtime parameter mpole_subSample
  call RuntimeParameters_get("mpole_subSample", mpole_subSample)
  call Logfile_stamp(gr_myPE,mpole_subSample,"[gr_mpoleInit] Potential and Moment subsampling set to ")
  mpole_subSampleInv = 1.0 / float(mpole_subSample)



#ifdef FLASH_GRID_PARAMESH
  call RuntimeParameters_get("lrefine_max", lrefine_max)
#else
  lrefine_max=1
#endif
  call RuntimeParameters_get("Nblockx", nBlockx)
  call RuntimeParameters_get("Nblocky", nBlocky)
  call RuntimeParameters_get("Nblockz", nBlockz)

  call RuntimeParameters_get("xmin", xmin)
  call RuntimeParameters_get("xmax", xmax)
  call RuntimeParameters_get("ymin", ymin)
  call RuntimeParameters_get("ymax", ymax)
  call RuntimeParameters_get("zmin", zmin)
  call RuntimeParameters_get("zmax", zmax)

  call RuntimeParameters_get("point_mass", point_mass)
  call RuntimeParameters_get("point_mass_rsoft", point_mass_rsoft)
  call RuntimeParameters_get("grv_outdsfac",grv_outdsfac)
  call RuntimeParameters_get("grv_outrad",grv_outrad)


  ! Check if we support the requested grid geometry.
  ! Only bounded grid geometries are supported, and
  ! support for the two which use angular coordinates
  ! (3D cylindrical and 3D spherical) are deferred for now.

  if ((NDIM == 3) .and. (gr_geometry == CARTESIAN)) then

     if (axisym) then
        mpole_geometry = G_3DAXISYMMETRIC
        mpole_mmax = 0
        call Logfile_stamp(gr_myPE, '3D axisymmetry, ignoring m > 0 moments', '[gr_mpoleInit]')
     else 
        mpole_geometry = G_3DCARTESIAN
        mpole_mmax = mpole_lmax
     endif
 
  else if (NDIM == 2) then
     if (gr_geometry == CYLINDRICAL) then
        
        mpole_geometry = G_2DCYLINDRICAL
        mpole_mmax = 0
        call Logfile_stamp(gr_myPE, '2D axisymmetry, ignoring m > 0 moments', '[gr_mpoleInit]')

     elseif (gr_geometry == SPHERICAL) then

       if (       point_mass.ne.0.e0           &
            .and. point_mass_rsoft.eq.0.e0     &
            .and. xmin.eq.0.e0 ) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: non-zero point_mass_rsoft required for point mass located at Rmin=0')
       end if

        mpole_geometry = G_2DSPHERICAL
        mpole_mmax = 0
        call Logfile_stamp(gr_myPE, '2D spherical, ignoring m > 0 moments', '[gr_mpoleInit]')
     else
        call Driver_abortFlash ('[gr_mpoleInit] ERROR: unsupported geometry')
     end if
  else if ((NDIM == 1) .and. (gr_geometry == SPHERICAL)) then
     
     mpole_geometry = G_1DSPHERICAL
     mpole_lmax = 0
     mpole_mmax = 0
     call Logfile_stamp &
          (gr_myPE, '1d spherical symmetry, ignoring l > 0 moments', '[gr_mpoleInit]')
     
  else
     
     call Driver_abortFlash ('[gr_mpoleInit] ERROR: unsupported geometry')
     
  endif
  
  ! we are only allowed to do a quadrant (i.e. enforce reflection symmetry
  ! about y=0) if we are in 2-d cylindrical or 2-d spherical coords

  if (quadrant) then
     if((mpole_geometry /= G_2DCYLINDRICAL).AND.(mpole_geometry /= G_2DSPHERICAL)) then
        call Driver_abortFlash('[gr_mpoleInit] ERROR: quadrant only allowed in 2-d angular geometry')
     endif
  end if

  ! octant is for 3-d cartesian only -- enforce that mpole_lmax = 0 for now
  if (octant) then
     if (mpole_geometry /= G_3DCARTESIAN) then
        call Driver_abortFlash('[gr_mpoleInit] ERROR: octant only allowed in 3-d Cartesian geometry')
     endif

     call Logfile_stamp(gr_myPE, '3D octant -- ignoring l > 0 moments', '[gr_mpoleInit]')
    
     mpole_lmax = 0
     mpole_mmax = 0
  endif


!               Maximum number of zones across each dimension, if each
!               dimension were to become fully refined.  Also compute
!               minimum zone spacings at the maximum refinement level.
!               Assume the domain is a Cartesian box.

  Nxmax = Nblockx * NXB * 2.e0**(lrefine_max-1)

  dxmin = (xmax - xmin) / Nxmax

  if (NDIM >= 2) then
     Nymax = Nblocky * NYB * 2.e0**(lrefine_max-1)
     dymin = (ymax - ymin) / Nymax
  else
     Nymax = 0.
     dymin = 1.
  endif

  if (NDIM == 3) then
     Nzmax = Nblockz * NZB * 2.e0**(lrefine_max-1)
     dzmin = (zmax - zmin) / Nzmax
  else 
     Nzmax = 0.
     dzmin = 1.
  endif

!               Number of radial samples to use in moment arrays.
!               Inverse of sample spacing to use for moment arrays.


  ds = (dxmin*dymin*dzmin)**(1./NDIM) / mpole_subSample / 2.d0
  !dsinv = 1.d0/ds
  if(mpole_geometry /= G_2DSPHERICAL) then
     qmax = 2 * int( sqrt(Nxmax**2 + Nymax**2 + Nzmax**2) ) * mpole_subSample
     qmax = qmax + max( 4, nint(qmax*0.01) )
     maxrad = (qmax - 1) * ds
     !qmax = 2 * qmax ! ADDED BY JFG, NEEDED FOR FORCE ON BH.
  else
     qmax = Nxmax + 1
  end if
  allocate(rad(0:qmax))
  allocate(dsinvarr(0:qmax))

  rad(0) = 0.d0
  qmax = 1
  do while (1)
     rad(qmax) = rad(qmax-1) + ds
     if (rad(qmax) .gt. grv_outrad) exit
     if (rad(qmax) .gt. maxrad) exit
     qmax = qmax + 1
  enddo
  dsinvarr(0:qmax) = 1.d0/ds
  if (rad(qmax) .lt. maxrad) then
      ds = grv_outdsfac * ds
      qmax = qmax + 1
      qmid = qmax
      do while (1)
         rad(qmax) = rad(qmax-1) + ds
         if (rad(qmax) .gt. maxrad) exit
         qmax = qmax + 1
      enddo
      dsinvarr(qmid:qmax) = 1.d0/ds
  endif

  if (gr_myPE == MASTER_PE) then
     write (str_buffer,"(A)") 'using'
     write (int_to_str, "(I7)") qmax
     write (str_buffer(len_trim(str_buffer)+1:), "(A, A)") " ", trim(adjustl(int_to_str))
     write (str_buffer(len_trim(str_buffer)+1:), "(A)") " moment array,"
     write (int_to_str, "(I7)") qmax*2*2*(mpole_lmax+1)*(mpole_mmax+1)
     write (str_buffer(len_trim(str_buffer)+1:), "(A, A)") " ", trim(adjustl(int_to_str))
     write (str_buffer(len_trim(str_buffer)+1:), "(A)") " items"

     call Logfile_stamp(gr_myPE, str_buffer, '[gr_mpoleInit]')
     write(*,*)"[gr_mpoleInit] ",str_buffer
  endif


  ! Allocate moment arrays and other data structures.

  allocate ( Moment(0:qmax,1:2,1:2,0:mpole_lmax,0:mpole_mmax), stat=istat)

  allocate ( OldMoment(0:qmax,1:2,1:2,0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: Moment() allocate failed')

  if (.not. mpole_useMatrixMPI) then
     call Logfile_stampMessage(gr_myPE,"[mpole_moments] Using looped mpi_allreduce")
     allocate( Momtmp(0:qmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: Momtmp() allocate failed')
  else
     !  Anshu and Carlo believe this is only valid for mpole_mmax=0
    call Logfile_stampMessage(gr_myPE,"[mpole_moments] Using matrix mpi_allreduce")
    if (mpole_mmax .NE. 0) then
        call Logfile_stamp &
          (gr_myPE,mpole_mmax,"[gr_mpoleInit] Cannot use mpole_useMatrixMPI=true with mmax!=0,now is ")
        call Driver_abortFlash("[gr_mpoleInit] Cannot use mpole_useMatrixMPI=true with mmax!=0")
     end if

     allocate ( MomtmpMatrix(0:qmax,1:2,1:2,0:mpole_lmax,0:mpole_mmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: MomtmpMatrix() allocate failed')
  endif

  if(mpole_geometry /= G_2DSPHERICAL) then
     allocate( costable(0:mpole_mmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: costable() allocate failed')

     allocate( sintable(0:mpole_mmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: sintable() allocate failed')
     
     allocate( rpower(0:mpole_lmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: rpower() allocate failed')
     
     allocate( Legk1(0:mpole_lmax,0:mpole_mmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: Legk1() allocate failed')
     
     allocate( rprinv(0:mpole_lmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: rprinv() allocate failed')
     
     allocate( Legk2(0:mpole_lmax,0:mpole_mmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: Legk2() allocate failed')
     
     allocate( Leg_fact(0:mpole_lmax,0:mpole_mmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: leg_fact() allocate failed')

     do m = 0, mpole_mmax
        do l = m+2, mpole_lmax
           Legk1(l,m) = float(2*l - 1) / float(l - m)
           Legk2(l,m) = float(l + m - 1) / float(l - m)
        enddo
     enddo
     
     if ((mpole_lmax > 0) .and. (mpole_mmax > 0)) then
        
        do l = 1, mpole_lmax
           factrl = 2.
           do m = 1, l
              factrl = factrl / float((l+m) * (l-m+1))
              Leg_fact(l,m) = factrl
           enddo
        enddo
        
     endif
  else
     
     maxx = NXB+2*NGUARD
     maxy = NYB+2*NGUARD

     allocate( pleg(1:maxy,0:mpole_lmax+1), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: pleg() allocate failed')

     allocate( pint(1:maxy,0:mpole_lmax+1), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: pint() allocate failed')
     
     allocate( yzn(1:maxy), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: yzn() allocate failed')
     
     allocate( xzn(1:maxx), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: xzn() allocate failed')

     allocate( r2(0:qmax), stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: r2() allocate failed')
  
     allocate(gpot(1:NXB+2*NGUARD,1:NYB+2*NGUARD,1),stat=istat)
     if (istat>0) call Driver_abortFlash ('[gr_mpoleInit] ERROR: gpot() allocate failed')


     ! Moments here are fac1 and fac2 terms; fac2 has opposite order than
     ! in original version

     Moment(:,:,:,:,:) = 0.e0
     r2(0)             = xmin**2

     !tp better accuracy?
     do m = 1,qmax
        cxzn_1 = xmin + float(m-1)*dxmin ! left global zone interface
        cxzn   = xmin + float(m)  *dxmin ! right global zone interface
        Moment(m  ,1,2,:,0) = cxzn_1/cxzn
        Moment(m-1,2,2,:,0) = Moment(m,1,2,:,0)
        r2(m)               = cxzn**2
     end do

     !tp loop index only up to mpole_lmax
     do m = 0,mpole_lmax
        Moment(:,1,2,m,0) = Moment(:,1,2,m,0)**float(m+1)
        Moment(:,2,2,m,0) = Moment(:,1,2,m,0)**float(m)
     end do
  end if

  ! If we are doing only a quadrant of the star in 2-d cylindrical coords, then
  ! we pick up an extra factor of 2 in the volume calculations because of the 
  ! symmetry about the y = 0 axis
  
  if (quadrant) then
     cylfactor = 2.e0*twopi
  else
     cylfactor = twopi
  endif
  
  if (octant) then
     cartfactor = 8.e0
  else
     cartfactor = 1.e0
  endif
  

  ! Bunch of parameters 

  fourpi_inv = 1.0/fourpi
    
  !==============================================================================
  
  ! If we are doing only a quadrant of the star in 2-d spherical coords, then
  ! we pick up an extra factor of 2 in the volume calculations because of the 
  ! symmetry about the y = 0 axis
  
  if (quadrant) then
     lstep = 2
  else
     lstep = 1
  endif


  return
end subroutine gr_mpoleInit
