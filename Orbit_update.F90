subroutine Orbit_update
    use nrtype
    use nr
    use ode_path
    use Gravity_data, ONLY: grv_ptmass, grv_mode, orb_t, orb_dt, &
        grv_ptvec, grv_obvec, grv_optvec, grv_oobvec
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, oXcm, oYcm, oZcm

    implicit none
    integer, parameter :: nvar=8
    integer :: i
    double precision :: eps,h1,hmin,x1,x2
    double precision, dimension(nvar) :: ystart
    INTERFACE
        SUBROUTINE derivs(x,y,dydx)
        IMPLICIT NONE
        double precision, INTENT(IN) :: x
        double precision, DIMENSION(:), INTENT(IN) :: y
        double precision, DIMENSION(:), INTENT(OUT) :: dydx
        END SUBROUTINE derivs
    END INTERFACE
    x1=orb_t
    x2=orb_t+orb_dt

    if (grv_mode .eq. 1) then
        ! Save the old parameters
        grv_optvec = grv_ptvec
        grv_oobvec = grv_obvec
    else
        grv_ptvec = grv_optvec
        grv_obvec = grv_oobvec
    endif

    ystart(1:2) = grv_obvec(1:2)
    ystart(3:4) = grv_ptvec(1:2)
    ystart(5:6) = grv_obvec(4:5)
    ystart(7:8) = grv_ptvec(4:5)
    eps=1.0d-12
    h1=(x2 - x1) / 1.d6
    hmin=0.d0
    call odeint(ystart,x1,x2,eps,h1,hmin,derivs,bsstep)
    grv_obvec(1:2) = ystart(1:2)
    grv_ptvec(1:2) = ystart(3:4)
    grv_obvec(4:5) = ystart(5:6)
    grv_ptvec(4:5) = ystart(7:8)
    if (grv_mode .eq. 2) then
        grv_obvec(1) = grv_obvec(1) + Xcm - oXcm
        grv_obvec(2) = grv_obvec(2) + Ycm - oYcm
        grv_obvec(3) = grv_obvec(3) + Zcm - oZcm
    endif
end subroutine Orbit_update

subroutine derivs(x,y,dydx)
    use Gravity_data, ONLY: grv_ptmass, grv_mode, orb_t, orb_dt
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use Grid_interface, ONLY: Grid_getMinCellSize
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, oXcm, oYcm, oZcm

    implicit none
    double precision, INTENT(IN) :: x
    double precision, DIMENSION(:), INTENT(IN) :: y
    double precision, DIMENSION(:), INTENT(OUT) :: dydx
    double precision :: r, potl, potr, xx, yy, zz, xcmt, ycmt
    double precision :: newton, dx, ptax, ptay, ptaz, fac, nptax, nptay
    double precision :: xcm_accel, ycm_accel, nxcm_accel, nycm_accel

    call PhysicalConstants_get("Newton", newton)
    call Grid_getMinCellSize(dx)
    dx = dx/2.d0

    dydx(1)=y(5) ! Obj. vel.
    dydx(2)=y(6)
    dydx(3)=y(7) ! Pt. vel.
    dydx(4)=y(8)
    xx = y(3) - y(1)
    yy = y(4) - y(2)
    if (grv_mode .eq. 2) then
        fac = (x - orb_t)/orb_dt
        xx = xx - (Xcm - oXcm)*fac
        yy = yy - (Ycm - oYcm)*fac
    endif
    r = sqrt(xx**2.d0 + yy**2.d0)
    if (xx .ne. xx .or. yy .ne. yy .or. dx .ne. dx) then
        print *, y(1:8)
        call Driver_abortFlash('failure')
    endif
    call gr_zonePotential(xx - dx, yy, 0.d0, potl)
    call gr_zonePotential(xx + dx, yy, 0.d0, potr)
    ptax = newton*(potr - potl) / dx / 2.d0
    call gr_zonePotential(xx, yy - dx, 0.d0, potl)
    call gr_zonePotential(xx, yy + dx, 0.d0, potr)
    ptay = newton*(potr - potl) / dx / 2.d0
    if (ptax .ne. ptax .or. ptay .ne. ptay) then
        print *, ptax, ptay, potr, potl, xx, yy, dx, newton, fac, x, orb_t, orb_dt
        call Driver_abortFlash('Force is NaN!')
    endif
    if (grv_mode .eq. 2) then
        call gr_zoneOldPotential(xx - dx, yy, 0.d0, potl)
        call gr_zoneOldPotential(xx + dx, yy, 0.d0, potr)
        nptax = newton*(potr - potl) / dx / 2.d0
        ptax = nptax + (ptax - nptax)*fac
        call gr_zoneOldPotential(xx, yy - dx, 0.d0, potl)
        call gr_zoneOldPotential(xx, yy + dx, 0.d0, potr)
        nptay = newton*(potr - potl) / dx / 2.d0
        ptay = nptay + (ptay - nptay)*fac
    endif
    if (ptax .eq. 0.0 .or. ptay .eq. 0.0) then
        call Driver_abortFlash('Force is zero in Orbit_update!')
    endif
    dydx(5)=newton/r**3.d0*grv_ptmass*xx! - xcm_accel
    dydx(6)=newton/r**3.d0*grv_ptmass*yy! - ycm_accel
    dydx(7)=ptax
    dydx(8)=ptay
end subroutine derivs
