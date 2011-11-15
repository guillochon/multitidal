! This routine updates the position of the two tracer particles which track the
! locations of the point mass and extended object. This routine can run in one
! of two modes depending on the grv_orb3D parameter. If this parameter is set to
! false, then the code ignores the z-direction when integrating the orbits. If
! set to true, the code includes all three dimensions. To avoid very small
! values for the second derivatives in the z-direction (which is usually the
! case when the orbital angular momentum is entirely in the z-direction), the
! coordinate system is rotated 45 degrees, integrated forward in time, and then
! rotated back into the original coordinate system. In the future this rotation
! should be based on the orbit's initial angular momentum vector to ensure that
! any arbitrary initial orbit has a non-zero force in all three cartesian
! directions at all times.

subroutine Orbit_update()

    use nrtype
    use nr
    use ode_path
    use Gravity_data, ONLY: grv_ptmass, grv_mode, orb_t, orb_dt, &
        grv_ptvec, grv_obvec, grv_optvec, grv_oobvec, grv_hptvec, &
        grv_hobvec, grv_exactvec, grv_oexactvec, grv_orbTol, grv_orb3D
    use gr_mpoleData, ONLY: X_centerofmass, Y_centerofmass, Z_centerofmass
    use gr_isoMpoleData, ONLY: Xcm, Ycm, Zcm

    implicit none
    integer :: i
    double precision :: h1,hmin,x1,x2
    double precision, dimension(:), allocatable :: ystart
    double precision, dimension(3,3) :: roty, rotx, rot
    integer, dimension(2) :: ashape = (/ 3, 3 /)
    double precision, parameter :: rangle = PI/180.d0 * 45.d0
    INTERFACE
        SUBROUTINE derivs(x,y,dydx)
        IMPLICIT NONE
        double precision, INTENT(IN) :: x
        double precision, DIMENSION(:), INTENT(IN) :: y
        double precision, DIMENSION(:), INTENT(OUT) :: dydx
        END SUBROUTINE derivs
    END INTERFACE

#include "Flash.h"

    if (grv_orb3D) then
        allocate(ystart(12))

        roty = reshape((/ dcos(rangle), 0.d0,          -dsin(rangle), &
                          0.d0,         1.d0,          0.d0,          &
                          dsin(rangle), 0.d0,          dcos(rangle)   /), ashape)
        rotx = reshape((/ 1.d0,         0.d0,          0.d0,          & 
                          0.d0,         dcos(rangle),  dsin(rangle),  &
                          0.d0,         -dsin(rangle), dcos(rangle)   /), ashape)
        rot = matmul(roty, rotx)
    else
        allocate(ystart(8))
    endif
    x1=orb_t
    x2=orb_t+orb_dt

    if (grv_mode .eq. 1) then
        ! Save the old parameters
        grv_optvec = grv_ptvec
        grv_oobvec = grv_obvec
    elseif (grv_mode .eq. 3) then
        grv_ptvec = grv_optvec
        grv_obvec = grv_oobvec
    endif

    if (grv_orb3D) then
        ystart(1:3) = grv_obvec(1:3)
        ystart(4:6) = grv_ptvec(1:3)
        ystart(7:9) = grv_obvec(4:6)
        ystart(10:12) = grv_ptvec(4:6)

        ystart(1:3)   = matmul(rot,ystart(1:3))
        ystart(4:6)   = matmul(rot,ystart(4:6))
        ystart(7:9)   = matmul(rot,ystart(7:9))
        ystart(10:12) = matmul(rot,ystart(10:12))
    else
        ystart(1:2) = grv_obvec(1:2)
        ystart(3:4) = grv_ptvec(1:2)
        ystart(5:6) = grv_obvec(4:5)
        ystart(7:8) = grv_ptvec(4:5)
    endif
    h1=(x2 - x1) / 1.d4
    hmin=0.d0
    call odeint(ystart,x1,x2,grv_orbTol,h1,hmin,derivs,bsstep)
    if (grv_orb3D) then
        roty = reshape((/ dcos(rangle),  0.d0,         dsin(rangle),  &
                          0.d0,          1.d0,         0.d0,          &
                          -dsin(rangle), 0.d0,         dcos(rangle)  /), ashape)
        rotx = reshape((/ 1.d0,          0.d0,         0.d0,          & 
                          0.d0,          dcos(rangle), -dsin(rangle), &
                          0.d0,          dsin(rangle), dcos(rangle)   /), ashape)
        rot = matmul(rotx, roty)

        ystart(1:3)   = matmul(rot,ystart(1:3))
        ystart(4:6)   = matmul(rot,ystart(4:6))
        ystart(7:9)   = matmul(rot,ystart(7:9))
        ystart(10:12) = matmul(rot,ystart(10:12))

        grv_obvec(1:3) = ystart(1:3)
        grv_ptvec(1:3) = ystart(4:6)
        grv_obvec(4:6) = ystart(7:9)
        grv_ptvec(4:6) = ystart(10:12)
    else
        grv_obvec(1:2) = ystart(1:2)
        grv_ptvec(1:2) = ystart(3:4)
        grv_obvec(4:5) = ystart(5:6)
        grv_ptvec(4:5) = ystart(7:8)
    endif

    if (grv_mode .eq. 1) then
        if (grv_orb3D) then
            grv_hobvec(1:3) = ystart(1:3)
            grv_hptvec(1:3) = ystart(4:6)
            grv_hobvec(4:6) = ystart(7:9)
            grv_hptvec(4:6) = ystart(10:12)
        else
            grv_hobvec(1:2) = ystart(1:2)
            grv_hptvec(1:2) = ystart(3:4)
            grv_hobvec(4:5) = ystart(5:6)
            grv_hptvec(4:5) = ystart(7:8)
        endif
    elseif (grv_mode .eq. 3) then
        grv_obvec(1:3) = grv_obvec(1:3) - grv_oexactvec(1:3) + grv_exactvec(1:3)
    endif

    deallocate(ystart)
end subroutine Orbit_update

subroutine derivs(x,y,dydx)
    use Gravity_data, ONLY: grv_ptmass, grv_mode, orb_t, orb_dt, grv_exactvec, &
        grv_orbMinForce, grv_oexactvec, grv_totmass, grv_orb3D
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use Grid_interface, ONLY: Grid_getMinCellSize
    use gr_mpoleData, ONLY: X_centerofmass, Y_centerofmass, Z_centerofmass, &
        max_R, &
        zone_max_radius_fraction, max_radial_zones
    use gr_isoMpoleData, ONLY: Xcm, Ycm, Zcm
    use gr_mpoleInterface, ONLY: gr_mpoleGradTotPot, gr_mpoleGradTotOldPot
    use Driver_data, ONLY: dr_simTime
    use Grid_data, ONLY: gr_meshMe

    implicit none
#include "Flash.h"
#include "constants.h"

    double precision, INTENT(IN) :: x
    double precision, DIMENSION(:), INTENT(IN) :: y
    double precision, DIMENSION(:), INTENT(OUT) :: dydx
    double precision :: xcmt, ycmt
    double precision :: newton, fac
    double precision :: xcm_accel, ycm_accel, nxcm_accel, nycm_accel
    double precision, dimension(3) :: grad_pot, ptt0, dist
    double precision :: last_zone_fraction
    double precision :: max_dydx
    double precision, dimension(3,3) :: rotx, roty, rot
    integer, dimension(2) :: ashape = (/ 3, 3 /)
    double precision, parameter :: rangle = PI/180.d0 * 45.d0

    call PhysicalConstants_get("Newton", newton)

    if (grv_orb3D) then
        roty = reshape((/ dcos(rangle),  0.d0, dsin(rangle), &
                          0.d0,          1.d0, 0.d0,         &
                          -dsin(rangle), 0.d0, dcos(rangle)  /), ashape)
        rotx = reshape((/ 1.d0,          0.d0,         0.d0,          & 
                          0.d0,          dcos(rangle), -dsin(rangle), &
                          0.d0,          dsin(rangle), dcos(rangle)   /), ashape)
        rot = matmul(rotx, roty)

        dydx(1)=y(7) ! Obj. vel.
        dydx(2)=y(8)
        dydx(3)=y(9)
        dydx(4)=y(10) ! Pt. vel.
        dydx(5)=y(11)
        dydx(6)=y(12)
        dist = matmul(rot,y(4:6) - y(1:3))
    else
        dydx(1)=y(5) ! Obj. vel.
        dydx(2)=y(6)
        dydx(3)=y(7) ! Pt. vel.
        dydx(4)=y(8)
        dist = (/ y(3) - y(1), y(4) - y(2), 0.d0 /)
    endif
    if (grv_mode .eq. 3) then
        fac = (x - orb_t)/orb_dt
        dist(1:2) = dist(1:2) - (grv_exactvec(1:2) - grv_oexactvec(1:2))*fac
        if (grv_orb3D) then
            dist(3) = dist(3) - (grv_exactvec(3) - grv_oexactvec(3))*fac
        endif
    endif
    last_zone_fraction = zone_max_radius_fraction (max_radial_zones)
    if (sqrt(sum(dist**2.d0)) .gt. 0.99*max_R*last_zone_fraction) then
        print *, grv_mode, grv_exactvec, grv_oexactvec
        print *, x, orb_t, orb_dt
        print *, y, dydx
        print *, dist, sqrt(sum(dist**2.d0)), max_R*last_zone_fraction
        call Driver_abortFlash('ERROR: Point mass is beyond outermost radial zone!')
    endif
    call gr_mpoleGradTotPot(dist, grad_pot)
    !if (maxval(abs(grad_pot)) .gt. 1.d20) then
    !    print *, dist, grad_pot
    !    call Driver_abortFlash('Force too large, something bad happened in gr_mpoleGradTotPot!')
    !endif
    if (grv_orb3D) then
        roty = reshape((/ dcos(rangle), 0.d0,          -dsin(rangle), &
                          0.d0,         1.d0,          0.d0,          &
                          dsin(rangle), 0.d0,          dcos(rangle)   /), ashape)
        rotx = reshape((/ 1.d0,         0.d0,          0.d0,          & 
                          0.d0,         dcos(rangle),  dsin(rangle),  &
                          0.d0,         -dsin(rangle), dcos(rangle)   /), ashape)
        rot = matmul(roty, rotx)

        dydx(10:12) = matmul(rot,grad_pot)
    else
        dydx(7:8) = grad_pot(1:2)
    endif
    if (grv_mode .eq. 3) then
        call gr_mpoleGradTotOldPot(dist, grad_pot)
        !if (maxval(abs(grad_pot)) .gt. 1.d20) then
        !    print *, dist, grad_pot
        !    call Driver_abortFlash('Force too large, something bad happened in gr_mpoleGradTotOldPot!')
        !endif
        ptt0 = grad_pot
        if (grv_orb3D) then
            ptt0 = matmul(rot,ptt0)
            dydx(10:12) = ptt0 + (dydx(10:12) - ptt0)*fac
        else
            dydx(7:8) = ptt0(1:2) + (dydx(7:8) - ptt0(1:2))*fac
        endif
    endif

    if (grv_orb3D) then
        dydx(7:9) = -dydx(10:12)*grv_ptmass/grv_totmass
    else
        dydx(5:6) = -dydx(7:8)*grv_ptmass/grv_totmass
    endif

    !if (dydx(10) .ne. dydx(10) .or. dydx(11) .ne. dydx(11) .or. dydx(12) .ne. dydx(12)) then
    !    print *, 'dydx', dydx(10:12)
    !    print *, 'fac', fac
    !    print *, 'x, orb_t, orb_dt', x, orb_t, orb_dt
    !    call Driver_abortFlash('Force is NaN!')
    !endif

    ! Make sure none of the forces are much smaller than the maximum force
    if (grv_orb3D) then
        max_dydx = maxval(abs(dydx(7:12)))
        if (abs(dydx(7))  .lt. grv_orbMinForce*max_dydx) dydx(7)  = 0.d0
        if (abs(dydx(8))  .lt. grv_orbMinForce*max_dydx) dydx(8)  = 0.d0
        if (abs(dydx(9))  .lt. grv_orbMinForce*max_dydx) dydx(9)  = 0.d0
        if (abs(dydx(10)) .lt. grv_orbMinForce*max_dydx) dydx(10) = 0.d0
        if (abs(dydx(11)) .lt. grv_orbMinForce*max_dydx) dydx(11) = 0.d0
        if (abs(dydx(12)) .lt. grv_orbMinForce*max_dydx) dydx(12) = 0.d0
    else
        max_dydx = maxval(abs(dydx(5:8)))
        if (abs(dydx(5))  .lt. grv_orbMinForce*max_dydx) dydx(5)  = 0.d0
        if (abs(dydx(6))  .lt. grv_orbMinForce*max_dydx) dydx(6)  = 0.d0
        if (abs(dydx(7))  .lt. grv_orbMinForce*max_dydx) dydx(7)  = 0.d0
        if (abs(dydx(8))  .lt. grv_orbMinForce*max_dydx) dydx(8)  = 0.d0
    endif
end subroutine derivs
