subroutine calc_orbit (t, m1, m2, ob1vec, ob2vec)
    use Gravity_data, ONLY: orb_ecc
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    implicit none
    double precision, intent(IN) :: t, m1, m2
    double precision, dimension(6), intent(INOUT) :: ob1vec, ob2vec

    if (orb_ecc .eq. 1.d0) then
        call parabolic_orbit(t, m1, m2, ob1vec, ob2vec)
    else
        call elliptical_orbit(t, m1, m2, ob1vec, ob2vec)
    endif
    
    return
end subroutine calc_orbit

subroutine parabolic_orbit (t, m1, m2, ob1vec, ob2vec)
    use Gravity_data, ONLY: grv_factor, peri_dist, peri_time
    implicit none

    double precision, intent(IN)  :: t, m1, m2
    double precision :: newx, newy, newz, velx, vely, velz
    double precision :: gm, u, r, q, tmt, term1, newton
    double precision, dimension(6), intent(INOUT) :: ob1vec, ob2vec
    
    call PhysicalConstants_get("Newton", newton)
    gm = newton * (m1 + m2)
    q = peri_dist
    tmt = t - peri_time
    term1 = (sqrt(8.*q**3.+9.*gm*tmt**2.)+3.*sqrt(gm)*tmt)
    u = (term1**(2./3.)-2.*q) / (sqrt(2.*q)*term1**(1./3.))
    u = 2*atan(u)
    r = 2*q/(1+cos(u))
    newx = -r*sin(u)
    newy = r*cos(u)
    newz = 0.d0

    ob1vec(1) = -newx * m2 / (m1 + m2)
    ob1vec(2) = -newy * m2 / (m1 + m2)
    ob1vec(3) = -newz * m2 / (m1 + m2)
    ob1vec(4) = -velx * m2 / (m1 + m2)
    ob1vec(5) = -vely * m2 / (m1 + m2)
    ob1vec(6) = -velz * m2 / (m1 + m2)
    ob2vec(1) = newx * m1 / (m1 + m2)
    ob2vec(2) = newy * m1 / (m1 + m2)
    ob2vec(3) = newz * m1 / (m1 + m2)
    ob2vec(4) = velx * m1 / (m1 + m2)
    ob2vec(5) = vely * m1 / (m1 + m2)
    ob2vec(6) = velz * m1 / (m1 + m2)
    return
end

subroutine elliptical_orbit (t, m1, m2, ob1vec, ob2vec)
    use Gravity_data, ONLY: grv_ptmass, peri_dist, peri_time, orb_ecc, grv_bound
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    implicit none

    double precision, parameter :: pi = 3.1415926535897932d0
    double precision, intent(IN) :: t, m1, m2
    double precision, dimension(6), intent(INOUT) :: ob1vec, ob2vec
    double precision :: gm, u, r, a, P, eccanom, ceccanom, newton
    double precision :: newx, newy, newz, velx, vely, velz

    call PhysicalConstants_get("Newton", newton)
    gm = newton*(m1 + m2)
    a = peri_dist/(1.d0 - orb_ecc)
    P = sqrt(4.d0*pi**2.d0/gm*a**3.d0)
    call keplereq(2.d0*pi/P*(t - peri_time), orb_ecc, eccanom)
    ceccanom = cos(eccanom)
    u = acos((ceccanom - orb_ecc)/(1.d0 - orb_ecc*ceccanom))
    r = peri_dist*(1.d0 + orb_ecc)/(1.d0+orb_ecc*cos(u))
    newx = -sign(r*sin(u), eccanom)
    newy = r*cos(u)
    newz = 0.d0
    velx = orb_ecc + cos(u)
    vely = sign(sin(u), eccanom)
    velx = -velx * sqrt(gm/a)/sqrt(1.d0-orb_ecc**2.d0)
    vely = -vely * sqrt(gm/a)/sqrt(1.d0-orb_ecc**2.d0)
    velz = 0.d0

    ob1vec(1) = -newx * m2 / (m1 + m2)
    ob1vec(2) = -newy * m2 / (m1 + m2)
    ob1vec(3) = -newz * m2 / (m1 + m2)
    ob1vec(4) = -velx * m2 / (m1 + m2)
    ob1vec(5) = -vely * m2 / (m1 + m2)
    ob1vec(6) = -velz * m2 / (m1 + m2)
    ob2vec(1) = newx * m1 / (m1 + m2)
    ob2vec(2) = newy * m1 / (m1 + m2)
    ob2vec(3) = newz * m1 / (m1 + m2)
    ob2vec(4) = velx * m1 / (m1 + m2)
    ob2vec(5) = vely * m1 / (m1 + m2)
    ob2vec(6) = velz * m1 / (m1 + m2)

    return
end subroutine elliptical_orbit

subroutine keplereq(m, ecc, eccanom)
    implicit none
    double precision, parameter :: pi = 3.1415926535897932d0
    double precision, parameter :: c1_6 = 1.d0/6.d0
    double precision, parameter :: c1_24 = 1.d0/24.d0
    double precision, parameter :: thresh = 1.d-10
    double precision, intent(IN) :: m, ecc
    double precision, intent(OUT) :: eccanom
    double precision :: mx, aux, alpha, beta, z, test, s0, s1, e0, se0, ce0
    double precision :: f, f1, f2, f3, f4, u1, u2, u3, u4, oldval, fe, fs, diff, mmm

    ! Range reduction of m to -pi < m <= pi
    ! ... m > pi
    mx = m

    if (mx .gt. pi) mx = mod(mx, 2.d0*pi)
    if (mx .gt. pi) mx = mx - 2.d0*pi

    ! ... m < -pi
    if (mx .le. -pi) mx = mod(mx, 2.d0*pi)
    if (mx .le. -pi) mx = mx + 2.d0*pi

    !
    ! Bail out for circular orbits...
    !
    if (ecc .eq. 0.d0) then
        eccanom = mx
        return
    endif

    aux   =  4.d0*ecc+0.5d0
    alpha = (1.d0-ecc)/aux

    beta=mx/(2.d0*aux)
    aux=sqrt(beta**2.d0+alpha**3.d0)
     
    z=beta+aux
    if (z .le. 0.d0) z = beta - aux

    test=abs(z)**(1.d0/3.d0)

    z = test
    if (z .lt. 0.d0) z = -z

    s0=z-alpha/z
    s1=s0-(0.078d0*sign(s0**5.d0,s0))/(1.d0+ecc)
    e0=mx+ecc*(3.d0*s1-4.d0*s1**3.d0)

    se0=sin(e0)
    ce0=cos(e0)

    f  = e0-ecc*se0-mx
    f1 = 1.d0-ecc*ce0
    f2 = ecc*se0
    f3 = ecc*ce0
    f4 = -f2
    u1 = -f/f1
    u2 = -f/(f1+0.5d0*f2*u1)
    u3 = -f/(f1+0.5d0*f2*u2+c1_6*f3*u2*u2)
    u4 = -f/(f1+0.5d0*f2*u3+c1_6*f3*u3*u3+c1_24*f4*u3**3.d0)

    eccanom=e0+u4

    if (eccanom .ge. 2.d0*pi) eccanom = eccanom - 2.d0*pi
    if (eccanom .lt. 0.d0) eccanom = eccanom + 2.d0*pi

    ! Now get more precise solution using Newton Raphson method
    ! for those times when the Kepler equation is not yet solved
    ! to better than 1e-10
    ! (modification J. Wilms)

    mmm=mx
    if (mmm .lt. 0.d0) mmm = mmm + 2.d0*pi
    diff=eccanom-ecc*sin(eccanom) - mmm

    if (abs(diff) .gt. 1.d-10) then
        do
            ! e-e sine-m
            fe=eccanom-ecc*sin(eccanom)-mmm
            ! f' = 1-e*cose
            fs=1.d0-ecc*cos(eccanom)
            oldval=eccanom
            eccanom=oldval-fe/fs
            if (abs(oldval-eccanom) .le. thresh) exit
        enddo 
        ! the following should be coded more intelligently ;-) 
        ! (similar to range reduction of mx...)
    endif 

    eccanom = mod(eccanom + pi, 2.d0*pi) - pi
    return
end subroutine keplereq
