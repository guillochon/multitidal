!!****if* source/Grid/GridSolvers/Multipole_experimental/gr_mpoleSetOuterZoneGrid
!!
!! NAME
!!
!!  gr_mpoleSetOuterZoneGrid
!!
!! 
!! SYNOPSIS
!!
!!  gr_mpoleSetOuterZoneGrid ()
!!
!!
!! DESCRIPTION
!!
!!  This routine sets up the outer (statistical) zone radial grid.
!!
!!
!!***

subroutine gr_mpoleSetOuterZoneGrid ()

  use Grid_data,         ONLY : gr_meshMe
  use gr_mpoleData,      ONLY : ZONE_EXP,                 &
                                ZONE_LOG,                 &
                                one,zero,                 &
                                dr,dr_inv,                &
                                max_R,max_Q,              &
                                max_radial_zones,         &
                                min_radial_zone,          &
                                zone_rmax,                &
                                zone_qmax,                &
                                zone_type,                &
                                zone_scalar,              &
                                zone_lognorm,             &
                                zone_exponent,            &
                                zone_scalar_inv,          &
                                zone_lognorm_inv,         &
                                zone_exponent_inv,        &
                                zone_max_radius_fraction, &
                                inner_zone_qmax,          &
                                inner_zone_rmax,          &
                                outer_zone_Qshift
  use Gravity_data, ONLY: grv_obvec, grv_ptvec, grv_mpolevec, grv_exactvec
  use Simulation_data, ONLY: sim_startDistance
  use IO_interface, ONLY: IO_getScalar
  use Driver_data, ONLY: dr_simTime, dr_initialSimTime, dr_restart, dr_dt
  use Grid_interface, ONLY: Grid_getMinCellSize

  implicit none

#include "Flash.h"
#include "constants.h"

  integer :: outer_zone_qmin
  integer :: Q, Q_local, Q_global
  integer :: type
  integer :: zone

  real    :: r_global, r_local
  real    :: scalar, lognorm, exponent
  real    :: scl_inv, lgn_inv, exp_inv

  real, dimension(2*MDIM) :: ptvec, obvec
  real, dimension(MDIM) :: mpolevec, exactvec
  real    :: pt_dist, pt_travel_dist, mcs
!
!
!       ...Set up the characteristic values for each specified radial zone
!          according to the linear and scaling factors given for each of these
!          zones. These characteristic arrays consist of:
!
!                a) maximum global radial values for each zone
!                b) maximum global bin values for each zone 
!
!          Global here means including all the corresponding values of lower
!          zones, as opposed to local values, which are specific to each zone.
!          Set also the overall maximum number of radial bins to be expected.
!
!          The outer (statistical) zones can be of two types: exponential or
!          logarithmic. Both cases use different equations to determine the
!          radial bin boundaries:
!
!                  Exponential  =>  maximum r for local Q = s * dr * Q^e
!                  Logarithmic  =>  maximum r for local Q = s * dr * [exp^(eQ) - 1]/[exp^e - 1]
!
!          A note is in place for the logarithmic scaling. It can easily be shown
!          that:
!
!                  [exp^(eQ) - 1]/[exp^e - 1] >= Q  for all e >= 0 and Q >= 1
!
!          The equality holds only for the two cases:
!
!                                 1) limit e --> 0
!                                 2)       Q  =  1
!
!          The proof of this is by induction starting with the Q = 2 expression for
!          the ratio: exp^e + 1. Why not use a different logarithmic scaling formula
!          like: s * dr * Q * exp^(e[Q-1]) ? The reason is that this formula is
!          very hard to invert for Q.
!
!          The pair of scalar/exponential values (s,e) and to which type each zone
!          belongs to can been specified by the user through runtime parameters. 
!          To determine the maximum local Q value for each zone, the above equations
!          must be inverted:
!
!                  Exponential  =>  maximum local Q for an r = [(r/(s*dr))^(1/e)]
!                  Logarithmic  =>  maximum local Q for an r = [(1/e)*log{r*(exp^e-1)/(s*dr) + 1}]
!
!          where [] denotes the ceiling function.
!
!
  zone_rmax (0) = zero
  zone_qmax (0) = zero

  ! Must be done because Gravity_init is after this routine in Driver_initFlash.
  call Grid_getMinCellSize(mcs)

  if (dr_simTime .eq. dr_initialSimTime) then
      if (dr_restart) then
          call IO_getScalar("ptxpos", ptvec(1))
          call IO_getScalar("ptypos", ptvec(2))
          call IO_getScalar("ptzpos", ptvec(3))
          call IO_getScalar("ptxvel", ptvec(4))
          call IO_getScalar("ptyvel", ptvec(5))
          call IO_getScalar("ptzvel", ptvec(6))
          call IO_getScalar("obxpos", obvec(1))
          call IO_getScalar("obypos", obvec(2))
          call IO_getScalar("obzpos", obvec(3))
          call IO_getScalar("obxvel", obvec(4))
          call IO_getScalar("obyvel", obvec(5))
          call IO_getScalar("obzvel", obvec(6))
          call IO_getScalar("grv_ompolevec_x",  mpolevec(1))
          call IO_getScalar("grv_ompolevec_y",  mpolevec(2))
          call IO_getScalar("grv_ompolevec_z",  mpolevec(3))
          call IO_getScalar("grv_oexactvec_x",  exactvec(1))
          call IO_getScalar("grv_oexactvec_y",  exactvec(2))
          call IO_getScalar("grv_oexactvec_z",  exactvec(3))
          pt_dist = dsqrt(sum((ptvec(1:3) - obvec(1:3))**2.d0))
          pt_travel_dist = dsqrt(sum((ptvec(4:6) - obvec(4:6))**2.d0))*dr_dt*2.d0
      else
          pt_dist = sim_startDistance
          pt_travel_dist = 4.d0*mcs
      endif
  else
      pt_dist = dsqrt(sum((grv_ptvec(1:3) - grv_obvec(1:3))**2.d0))
      pt_travel_dist = max(4.d0*mcs, dsqrt(sum((grv_obvec(4:6) - grv_ptvec(4:6))**2.d0))*dr_dt*2.d0)
  endif

  do zone = 1,max_radial_zones

     type     = zone_type         (zone)
     scl_inv  = zone_scalar_inv   (zone)
     exp_inv  = zone_exponent_inv (zone)
     if (zone .eq. max_radial_zones - 2) then
         r_global = pt_dist - pt_travel_dist
     elseif (zone .eq. max_radial_zones - 1) then
         r_global = pt_dist + pt_travel_dist
     elseif (zone .eq. 1) then
         r_global = max(zone_max_radius_fraction (zone) * max_R, &
             dsqrt(sum((grv_mpolevec(1:3) - grv_exactvec(1:3))**2.d0)))
     else
         r_global = zone_max_radius_fraction (zone) * max_R
     endif
     r_local  = r_global - zone_rmax (zone-1)

     if (type == ZONE_EXP) then
         Q_local = ceiling ( (r_local * scl_inv * dr_inv) ** exp_inv )
     else if (type == ZONE_LOG) then
         lgn_inv = zone_lognorm_inv (zone)
         Q_local = ceiling ( exp_inv * log (r_local * scl_inv * dr_inv * lgn_inv + one) )
     end if

     Q_global = zone_qmax (zone-1) + Q_local

     zone_rmax (zone) = r_global
     zone_qmax (zone) = Q_global

  end do

!
!
!       ...Find the first significant outer zone. All outer zones below
!          the first significant outer zone have been 'swallowed' by the
!          size of the inner zone. If no first significant outer zone
!          is found, i.e. if the complete domain is within the inner zone,
!          the first significant outer zone will be equal to zero for
!          further processing (see below).
!
!
  min_radial_zone = 0

  do zone = max_radial_zones,1,-1
     r_global = zone_rmax (zone)
     if (r_global > inner_zone_rmax) then
         min_radial_zone = zone
     end if
  end do
!
!
!       ...Determine within the first significant outer zone, which global
!          radial bin (outer_zone_qmin) first exceeds the maximum inner zone bin.
!          This gap has to be closed when merging the inner and outer zones
!          so that at their boundary the radial bin numbers continue without
!          interuption. Note, that the gap can be positive or negative:
!
!                   +ve : outer_zone_qmin > inner_zone_qmax
!                   -ve : outer_zone_qmin < inner_zone_qmax
!
!          In either case a radial bin shift value (outer_zone_Qshift) is
!          calculated, which brings the two zones together. This shift value
!          will be added to the outer zone radial bin values later on when
!          evaluating the moment bins.
!
!
  if (min_radial_zone == 0) then

      max_Q = inner_zone_qmax

  else

      Q_local  = zone_qmax (min_radial_zone) - zone_qmax (min_radial_zone - 1)

      type     = zone_type     (min_radial_zone)
      scalar   = zone_scalar   (min_radial_zone)
      exponent = zone_exponent (min_radial_zone)

      do Q = 1,Q_local

         if (type == ZONE_EXP) then
             r_local = scalar * dr * (real (Q) ** exponent)
         else if (type == ZONE_LOG) then
             lognorm = zone_lognorm (min_radial_zone)
             r_local = scalar * dr * lognorm * (exp (exponent * real (Q)) - one)
         end if

         r_global = zone_rmax (min_radial_zone - 1) + r_local

         if (r_global > inner_zone_rmax) then
             outer_zone_qmin = zone_qmax (min_radial_zone - 1) + Q
             exit
         end if

      end do

      outer_zone_Qshift = inner_zone_qmax - outer_zone_qmin + 1

      max_Q = zone_qmax (max_radial_zones) + outer_zone_Qshift

  end if
!
!
!       ...Ready!
!
!
  return
end subroutine gr_mpoleSetOuterZoneGrid
