!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleData
!!
!! NAME
!!
!!  gr_mpoleData
!!
!!
!! SYNOPSIS
!!
!!  use gr_mpoleData
!!
!!
!! DESCRIPTION
!!
!!  Variable declarations for the multipole Poisson solver.
!!
!!  Multipole module data structures -- these are allocated in mpole_init
!!  once the size of the mesh is computed.  These are then used to share 
!!  data across the various mpole functions
!!
!!    Moment(q,k,p,l,m)     The (l,m)-th moment (k=1, even; 2, odd;
!!                          p=1, inner; 2, outer) of the image
!!                          mass distribution as a fctn of q*dsinv
!!
!!    Momtmp(q)             Temporary array to receive radial
!!                          samples of Moment() in parallel
!!                          reduction operation
!!    MomtmpMatrix          Temporary array the size of Moment (.true. gr_useMatrixMPI)
!!
!!    mpole_useMatrixMPI     a runtime parameter to choose between the two methods
!!
!!    mpole_subSample        runtime parameter to control subsampling of potential and moment calculations
!!
!!    dsinv                 Moment array sample spacing
!!
!!    mpole_lmax            Maximum multipole moment (runtime set)
!!
!!    mpole_mmax            Maximum azimuthal moment (determined by geometry)
!!
!!    qmax                  Maximum number of radial samples
!!
!!    Mtot                  Total mass
!!
!!    X/Y/Zcm               Location of center of mass
!!
!!    mpole_geometry        Geometry type for grid
!!
!!    costable(m)           Table containing the cosine of
!!                          m * the current azimuthal angle; also
!!                          sintable(m)
!!
!!    rpower(l)             The current cell's density*rprime^l
!!
!!    rprinv(l)             The current cell's density*rprime^-(l+1)
!!
!!    Leg_fact(l,m)         Factorial normalization coefficients
!!                          for the associated Legendre function
!!    pleg(nyb,l)           legendre for spherical
!!    pint(nyb,l)
!!   
!! PARAMETERS
!!
!!  mpole_subSample  -- integer to control number of subzones in each direction
!!                           for smoothing potential calculations.
!!                      Also slows down calculation considerably when != 1
!!
!!***

module gr_mpoleData

  implicit none
  
  real, save    :: xmin, xmax, ymin, ymax, zmin, zmax
  real, save    :: dxmin, dymin, dzmin

  integer, save :: lrefine_max, nblockx, nblocky, nblockz

  real, save    :: Newton, point_mass, point_mass_rsoft, gbnd

  integer, save :: maxx, maxy, maxz, lstep

! Physical/mathematical constants -- these are filling in init_mpole

  real, save    :: pi, twopi, fourpi


  real, save     :: Mtot, Xcm, Ycm, Zcm, vXcm, vYcm, vZcm, oXcm, oYcm, oZcm, dsinv, grv_outdsfac, grv_outrad

  real, allocatable, save  :: OldMoment(:,:,:,:,:), Moment(:,:,:,:,:), Momtmp(:)
  real, allocatable, save  :: MomtmpMatrix(:,:,:,:,:)
  logical, save            :: mpole_useMatrixMPI


  real, allocatable, save  :: costable(:), sintable(:)

  real, allocatable, save  :: rpower(:), rprinv(:), dsinvarr(:), rad(:)

  real, allocatable, save  :: Legk1(:,:), Legk2(:,:), Leg_fact(:,:)

  real, allocatable, save  :: yzn(:), xzn(:), pleg(:,:),pint(:,:)

  real, allocatable, save  :: gpot(:,:,:) 

  real, allocatable, save  :: r2(:)

  integer, save            :: qmax, mpole_lmax, mpole_mmax, mpole_geometry

  integer, save            :: mpole_subSample
!  integer, save            :: Nint6  !! Was value of 6, now controlled by mpole_subSample
  real, save               :: mpole_subSampleInv  

  real, save :: fourpi_inv

! Constants used to index moment array


  integer, parameter  :: Even = 1, Odd = 2, Inner = 1, Outer = 2

! Supported geometry constants -- we do a geometry check to make sure
! that we support whatever the request geometry is in init_mpole

  integer, parameter :: G_3DCARTESIAN  = 1, G_2DCYLINDRICAL = 2, &
                        G_1DSPHERICAL  = 3, G_2DSPHERICAL   = 4, &
                        G_3DCYLINDRICAL= 5, G_3DAXISYMMETRIC = 6

! Supported boundary constants.

  integer, parameter :: MP_BND_ISOLATED = 0

! in 2-d cylindrical and 2-d spherical coordinates, we allow a single quadrant of 
! the star to be simulated if quadrant = .true.
! mpole_dumpMoments controls the output of the Moment array to a text file

  logical, save :: quadrant, octant, axisym, mpole_dumpMoments

  real, save    :: cylfactor, cartfactor

end module gr_mpoleData
