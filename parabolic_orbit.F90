subroutine parabolic_orbit (t, newx, newy, newz)
    use Gravity_data, ONLY: grv_factor, peri_dist, peri_time
    implicit none

    real, intent(IN)  :: t
    real, intent(OUT) :: newx, newy, newz
    real :: gm, u, r, q, tmt, term1
    
    gm = -grv_factor
    q = peri_dist
    tmt = t - peri_time
    term1 = (sqrt(8.*q**3.+9.*gm*tmt**2.)+3.*sqrt(gm)*tmt)
    u = (term1**(2./3.)-2.*q) / (sqrt(2.*q)*term1**(1./3.))
    u = 2*atan(u)
    r = 2*q/(1+cos(u))
    newx = -r*sin(u)
    newy = r*cos(u)
    newz = 0.
end

subroutine parabolic_velocity(vel,r)
    use Gravity_data, ONLY: grv_factor, peri_dist, peri_time
    use Driver_interface, ONLY: Driver_getSimTime
    implicit none

    real, intent(OUT) :: vel, r
    real :: t, gm, u, q, tmt, term1
    
    call Driver_getSimTime(t)
    gm = -grv_factor
    q = peri_dist
    tmt = t - peri_time
    term1 = (sqrt(8.*q**3.+9.*gm*tmt**2.)+3.*sqrt(gm)*tmt)
    u = (term1**(2./3.)-2.*q) / (sqrt(2.*q)*term1**(1./3.))
    u = 2.*atan(u)
    r = 2*q/(1+cos(u))
    vel = sqrt(2.*gm/r)
end
