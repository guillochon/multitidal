!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoMpoleData
!!
!! NAME
!!
!!  gr_isoMpoleData
!!
!!
!! SYNOPSIS
!!
!!  use gr_isoMpoleData
!!
!!
!! DESCRIPTION
!!
!!  Variable declarations for the multipole Poisson solver.
!!
!!  Multipole module data structures -- these are allocated in init_mpole 
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
!!
!!***

module gr_isoMpoleData

  implicit none

  real,save  ::  twopi, fourpi,fourpi_inv,Nint_inv
  integer, save :: Nint = 6

  real,save  :: Mtot, Xcm, Ycm, Zcm, dsinv

  real,save, allocatable  :: Moment(:,:,:,:,:), Momtmp(:)

  real,save, allocatable  :: costable(:), sintable(:)

  real,save, allocatable  :: rpower(:), rprinv(:)

  real,save, allocatable  :: Legk1(:,:), Legk2(:,:), Leg_fact(:,:)

  integer ,save           :: qmax, mpole_lmax, mpole_mmax, mpole_geometry, mpole_newton
  real, save         :: cylfactor


  ! Constants used to index moment array

  integer, parameter :: MPOLE_EVEN = 1, MPOLE_ODD = 2, MPOLE_INNER = 1, MPOLE_OUTER = 2


  ! Supported geometry constants -- we do a geometry check to make sure
  ! that we support whatever the request geometry is in init_mpole



  ! in 2-d cylindrical coordinates, we allow a single mpole_quadrant of the star
  ! to be simulated if mpole_quadrant = .true.

  logical,save :: mpole_quadrant,mpole_octant


end module gr_isoMpoleData
