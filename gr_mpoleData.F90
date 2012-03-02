!!****if* source/Grid/GridSolvers/Multipole_experimental/gr_mpoleData
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
!!    blockCount           Holds the number of leaf blocks on a particular
!!                         processor.
!!    blockList            Holds the list of leaf blocks on a particular
!!                         processor.
!!    bndBox               Holds the bounding box for a particular
!!                         block.
!!    blkLimits            Holds the block index limits for the interior
!!                         cells of a particular block.
!!    blkLimitsGC          Holds the general block index limits including
!!                         the guard cells of a particular block (not
!!                         needed for multipole Poisson solver, but
!!                         must be declared).
!!    delta                Holds the lengths of the individual cells
!!                         inside a block.
!!    Moment_R (n,q)       The n-th moment in the q-th radial bin
!!                         corresponding to the regular solid harmonic
!!                         function of the image mass distribution.
!!                         The order in index 'n' of the moments will
!!                         be:
!!                          
!!                            1) first all cosine parts grouped together
!!                               followed by all sine parts.
!!
!!                            2) within each cosine and sine parts:
!!                               increasing in M, i.e. M = 0,1,2,3,...
!!
!!                            3) within each M increasing in L, that
!!                               is L = M,M+1,...,max_L
!!
!!                         The maximum number of moments for a particular
!!                         q-index will be: (max_L + 1) * (max_L + 1)
!!
!!    Moment_I (n,q)       Same as for the Moment_R () array but corresponding
!!                         to the irregular solid harmonic functions.
!!
!!    Scratch (n,q)        Scratch array to reduce the individual processor
!!                         Moment_R () and Moment_I () data in parallel
!!                         reduction operation.
!!
!!    R_Moment (n)         Regular solid harmonic vector used for applying
!!                         the regular solid harmonic recursion relations
!!                         at a particular point in space. Same ordering and
!!                         size as the Moment_R () array for one particular
!!                         q-index.
!!
!!    I_Moment (n)         The same as for the R_Moment () array but for
!!                         the irregular solid harmonics.
!!
!!    Denorm_R (L)         Contains the denormalization corresponding to the
!!                         angular momentum L for evaluating the potential
!!                         scalar product contribution between normalized regular
!!                         solid harmonic functions and normalized irregular
!!                         moments.
!!
!!    Denorm_I (L)         Same as Denorm_R () but between normalized irregular
!!                         solid harmonic functions and normalized regular
!!                         moments.
!!
!!    Damping_R (q)        The damping factor for the R_Moment () array
!!                         corresponding to the q-th radial bin.
!!
!!    Damping_I (q)        Same as Damping_R () but for the I_Moment () array.
!!
!!    R_Damper (n)         Damping value for the n-th regular solid harmonic
!!                         function.
!!
!!    I_Damper (n)         The same as for the R_Damper () array but for
!!                         the irregular solid harmonic function.
!!
!!    dr                   smallest radial size
!!
!!    dr_inv               Inverse of smallest radial size
!!
!!    dr_inner_zone        smallest radial size in inner zone
!!
!!    dr_inner_zone_inv    Inverse of smallest radial size in inner zone
!!
!!    max_L                Maximum multipole moment (runtime set)
!!
!!    max_M                Maximum magnetic number (runtime set)
!!
!!    max_LM               Maximum number of angular moments (determined by
!!                         geometry)
!!
!!    max_R                Maximum radius for domain
!!
!!    max_Q                Maximum number of radial bins
!!
!!    Mtot                 Total mass
!!
!!    X/Y/Z_centerofmass   Location of center of mass
!!
!!    mpole_geometry       Geometry type for grid
!!   
!!    dumpMoments          if true, produces a detailed printout of the
!!                         moments at every iteration step
!!
!!    printRadialInfo      if true, produces a detailed printout of the
!!                         radial bin structure at every iteration step
!!  
!! PARAMETERS
!!
!!
!!***

module gr_mpoleData

  implicit none
  
  logical, save :: dumpMoments
  logical, save :: G_1D, G_2D, G_3D
  logical, save :: G_CARTESIAN, G_CYLINDRICAL, G_SPHERICAL, G_POLAR
  logical, save :: ignore_inner_zone
  logical, save :: printRadialInfo
  logical, save :: symmetry_plane_2D, symmetry_axis_3D

  integer, save :: blockCount
  integer, save :: inner_zone_qmax
  integer, save :: inner_zone_size
  integer, save :: max_Q, max_L, max_2L, max_M, max_LM
  integer, save :: max_radial_zones
  integer, save :: min_radial_zone
  integer, save :: mpole_geometry
  integer, save :: n_cosine_moments
  integer, save :: outer_zone_Qshift

  real,    save :: cube_potential_factor
  real,    save :: dr, dr_inv, dr_inner_zone, dr_inner_zone_inv
  real,    save :: G_constant
  real,    save :: inner_zone_grid, inner_zone_grid_inv
  real,    save :: inner_zone_rmax
  real,    save :: max_R
  real,    save :: Mtot
  real,    save :: X_centerofmass, Y_centerofmass, Z_centerofmass
  real,    save :: xmin, xmax
  real,    save :: ymin, ymax
  real,    save :: zmin, zmax
  real,    save :: zero, one, two, three, eight,            &
                   half, third, fourth, twelfth, hundredth, &
                   twopi, halfpi, thirdpi, sixthpi,         &
                   fourpi, fourpi_inv,                      &
                   ebase,ebase_inv

  integer, allocatable, save  :: blockList         (:)
  integer, allocatable, save  :: zone_qmax         (:)
  integer, allocatable, save  :: zone_type         (:)
  integer, allocatable, save  :: inner_zone_Qlower (:)
  integer, allocatable, save  :: inner_zone_Qupper (:)
  integer, allocatable, save  :: blkLimits         (:,:)
  integer, allocatable, save  :: blkLimitsGC       (:,:)

  real,    allocatable, save  :: delta                    (:)
  real,    allocatable, save  :: NumberInv                (:)
  real,    allocatable, save  :: R_Moment                 (:)
  real,    allocatable, save  :: I_Moment                 (:)
  real,    allocatable, save  :: D_Damper                 (:)
  real,    allocatable, save  :: R_Damper                 (:)
  real,    allocatable, save  :: I_Damper                 (:)
  real,    allocatable, save  :: Denorm_D                 (:)
  real,    allocatable, save  :: Denorm_R                 (:)
  real,    allocatable, save  :: Denorm_I                 (:)
  real,    allocatable, save  :: Damping_R                (:)
  real,    allocatable, save  :: Damping_I                (:)

  integer, save :: old_max_Q
  integer, save :: old_inner_zone_qmax
  integer, save :: old_min_radial_zone
  integer, save :: old_outer_zone_Qshift
  real,    save :: old_dr, old_dr_inv, old_dr_inner_zone, old_dr_inner_zone_inv
  real,    save :: old_inner_zone_rmax
  real,    save :: old_max_R
  integer, allocatable, save  :: old_inner_zone_Qlower    (:)
  integer, allocatable, save  :: old_inner_zone_Qupper    (:)
  real,    allocatable, save  :: old_inner_zone_radii     (:)
  integer, allocatable, save  :: old_zone_qmax                (:)
  real,    allocatable, save  :: old_zone_rmax                (:)

  real,    allocatable, save  :: Old_R_Moment             (:)
  real,    allocatable, save  :: Old_I_Moment             (:)
  real,    allocatable, save  :: Old_R_Damper             (:)
  real,    allocatable, save  :: Old_I_Damper             (:)
  real,    allocatable, save  :: Old_Damping_R            (:)
  real,    allocatable, save  :: Old_Damping_I            (:)
  real,    allocatable, save  :: Old_Denorm_R             (:)
  real,    allocatable, save  :: Old_Denorm_I             (:)
  real,    allocatable, save  :: Old_Moment_R             (:,:)
  real,    allocatable, save  :: Old_Moment_I             (:,:)

  !integer, save :: tot_max_Q
  !integer, save :: tot_inner_zone_qmax
  !integer, save :: tot_min_radial_zone
  !integer, save :: tot_outer_zone_Qshift
  !real,    save :: tot_dr, tot_dr_inv, tot_dr_inner_zone, tot_dr_inner_zone_inv
  !real,    save :: tot_inner_zone_rmax
  !real,    save :: tot_max_R
  !integer, allocatable, save  :: tot_inner_zone_Qlower    (:)
  !integer, allocatable, save  :: tot_inner_zone_Qupper    (:)
  !real,    allocatable, save  :: tot_inner_zone_radii     (:)
  !integer, allocatable, save  :: tot_zone_qmax                (:)
  !real,    allocatable, save  :: tot_zone_rmax                (:)

  !real,    allocatable, save  :: Tot_R_Moment             (:)
  !real,    allocatable, save  :: Tot_I_Moment             (:)
  !real,    allocatable, save  :: Tot_R_Damper             (:)
  !real,    allocatable, save  :: Tot_I_Damper             (:)
  !real,    allocatable, save  :: Tot_Damping_R            (:)
  !real,    allocatable, save  :: Tot_Damping_I            (:)
  !real,    allocatable, save  :: Tot_Denorm_R             (:)
  !real,    allocatable, save  :: Tot_Denorm_I             (:)
  !real,    allocatable, save  :: Tot_Moment_R             (:,:)
  !real,    allocatable, save  :: Tot_Moment_I             (:,:)

  !integer, save :: old_tot_max_Q
  !integer, save :: old_tot_inner_zone_qmax
  !integer, save :: old_tot_min_radial_zone
  !integer, save :: old_tot_outer_zone_Qshift
  !real,    save :: old_tot_dr, old_tot_dr_inv, old_tot_dr_inner_zone, old_tot_dr_inner_zone_inv
  !real,    save :: old_tot_inner_zone_rmax
  !real,    save :: old_tot_max_R
  !integer, allocatable, save  :: old_tot_inner_zone_Qlower    (:)
  !integer, allocatable, save  :: old_tot_inner_zone_Qupper    (:)
  !real,    allocatable, save  :: old_tot_inner_zone_radii     (:)
  !integer, allocatable, save  :: old_tot_zone_qmax                (:)
  !real,    allocatable, save  :: old_tot_zone_rmax                (:)

  !real,    allocatable, save  :: Old_Tot_R_Moment             (:)
  !real,    allocatable, save  :: Old_Tot_I_Moment             (:)
  !real,    allocatable, save  :: Old_Tot_R_Damper             (:)
  !real,    allocatable, save  :: Old_Tot_I_Damper             (:)
  !real,    allocatable, save  :: Old_Tot_Damping_R            (:)
  !real,    allocatable, save  :: Old_Tot_Damping_I            (:)
  !real,    allocatable, save  :: Old_Tot_Denorm_R             (:)
  !real,    allocatable, save  :: Old_Tot_Denorm_I             (:)
  !real,    allocatable, save  :: Old_Tot_Moment_R             (:,:)
  !real,    allocatable, save  :: Old_Tot_Moment_I             (:,:)

  real,    allocatable, save  :: zone_scalar              (:)
  real,    allocatable, save  :: zone_exponent            (:)
  real,    allocatable, save  :: zone_lognorm             (:)
  real,    allocatable, save  :: zone_scalar_inv          (:)
  real,    allocatable, save  :: zone_exponent_inv        (:)
  real,    allocatable, save  :: zone_lognorm_inv         (:)
  real,    allocatable, save  :: zone_max_radius_fraction (:)
  real,    allocatable, save  :: zone_rmax                (:)
  real,    allocatable, save  :: inner_zone_radii         (:)
  real,    allocatable, save  :: bndBox                   (:,:)
  real,    allocatable, save  :: Moment_R                 (:,:)
  real,    allocatable, save  :: Moment_I                 (:,:)
  real,    allocatable, save  :: Scratch                  (:,:)

  integer, parameter :: ZONE_EXP  =  1, &
                        ZONE_LOG  =  2

  integer, parameter :: G_1DCARTESIAN    =  1, &
                        G_2DCARTESIAN    =  2, &
                        G_3DCARTESIAN    =  3, &
                        G_1DCYLINDRICAL  =  4, &
                        G_2DCYLINDRICAL  =  5, &
                        G_3DCYLINDRICAL  =  6, &
                        G_1DSPHERICAL    =  7, &
                        G_2DSPHERICAL    =  8, &
                        G_3DSPHERICAL    =  9, &
                        G_1DPOLAR        = 10, &
                        G_2DPOLAR        = 11, &
                        G_3DPOLAR        = 12

end module gr_mpoleData
