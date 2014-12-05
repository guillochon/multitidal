subroutine calc_orbit (t, m1, m2, rp, tp, ecc, ob1vec, ob2vec)
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    implicit none
    double precision, intent(IN) :: t, m1, m2, rp, ecc, tp
    double precision, dimension(6), intent(INOUT) :: ob1vec, ob2vec

    if (ecc .eq. 1.d0) then
        call parabolic_orbit(t, m1, m2, rp, tp, ob1vec, ob2vec)
    else
        call elliptical_orbit(t, m1, m2, ecc, rp, tp, ob1vec, ob2vec)
    endif
    
    return
end subroutine calc_orbit

subroutine parabolic_orbit (t, m1, m2, rp, tp, ob1vec, ob2vec)
    implicit none

    double precision, intent(IN)  :: t, m1, m2, rp, tp
    double precision :: newx, newy, newz, velx, vely, velz
    double precision :: gm, u, r, tmt, term1, newton
    double precision, dimension(6), intent(INOUT) :: ob1vec, ob2vec
    
    call PhysicalConstants_get("Newton", newton)
    gm = newton * (m1 + m2)
    tmt = t - tp
    term1 = (dsqrt(8.*rp**3.+9.*gm*tmt**2.)+3.*dsqrt(gm)*tmt)
    u = (term1**(2./3.)-2.*rp) / (dsqrt(2.*rp)*term1**(1./3.))
    u = 2*datan(u)
    r = 2*rp/(1+dcos(u))
    newx = -r*dsin(u)
    newy = r*dcos(u)
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

subroutine elliptical_orbit (t, m1, m2, ecc, rp, tp, ob1vec, ob2vec)
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use Logfile_interface, ONLY: Logfile_stampMessage
    implicit none

#include "constants.h"

    double precision, intent(IN) :: t, m1, m2, ecc, rp, tp
    double precision, dimension(6), intent(INOUT) :: ob1vec, ob2vec
    double precision :: gm, u, r, a, P, eccanom, ceccanom, newton
    double precision :: newx, newy, newz, velx, vely, velz
    character(len=200) :: logstr

    call PhysicalConstants_get("Newton", newton)
    gm = newton*(m1 + m2)
    a = rp/(1.d0 - ecc)
    P = dsqrt(4.d0*PI**2.d0/gm*a**3.d0)
    call keplereq(2.d0*PI/P*(t - tp), ecc, eccanom)
    ceccanom = dcos(eccanom)
    u = dacos((ceccanom - ecc)/(1.d0 - ecc*ceccanom))
    r = rp*(1.d0 + ecc)/(1.d0+ecc*dcos(u))
    newx = -sign(r*dsin(u), eccanom)
    newy = r*dcos(u)
    newz = 0.d0
    velx = ecc + dcos(u)
    vely = sign(dsin(u), eccanom)
    velx = -velx * dsqrt(gm/a)/dsqrt(1.d0-ecc**2.d0)
    vely = -vely * dsqrt(gm/a)/dsqrt(1.d0-ecc**2.d0)
    velz = 0.d0

    write(logstr, fmt='(A30, 2ES15.8)') 'Bind ener #1, Bind ener #2:', -gm/2.d0/a, -gm/dsqrt(newx**2.d0 + newy**2.d0 + newz**2.d0) + (velx**2.d0 + vely**2.d0 + velz**2.d0)/2.d0
    call Logfile_stampMessage(logstr)

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
#include "constants.h"
    double precision, parameter :: c1_6 = 1.d0/6.d0
    double precision, parameter :: c1_24 = 1.d0/24.d0
    double precision, parameter :: thresh = 1.d-14
    double precision, intent(IN) :: m, ecc
    double precision, intent(OUT) :: eccanom
    double precision :: mx, aux, alpha, beta, z, test, s0, s1, e0, se0, ce0
    double precision :: f, f1, f2, f3, f4, u1, u2, u3, u4, oldval, fe, fs, diff, mmm

    ! Range reduction of m to -PI < m <= PI
    ! ... m > PI
    mx = m

    if (mx .gt. PI) mx = mod(mx, 2.d0*PI)
    if (mx .gt. PI) mx = mx - 2.d0*PI

    ! ... m < -PI
    if (mx .le. -PI) mx = mod(mx, 2.d0*PI)
    if (mx .le. -PI) mx = mx + 2.d0*PI

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
    aux=dsqrt(beta**2.d0+alpha**3.d0)
     
    z=beta+aux
    if (z .le. 0.d0) z = beta - aux

    test=dabs(z)**(1.d0/3.d0)

    z = test
    if (z .lt. 0.d0) z = -z

    s0=z-alpha/z
    s1=s0-(0.078d0*sign(s0**5.d0,s0))/(1.d0+ecc)
    e0=mx+ecc*(3.d0*s1-4.d0*s1**3.d0)

    se0=dsin(e0)
    ce0=dcos(e0)

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

    if (eccanom .ge. 2.d0*PI) eccanom = eccanom - 2.d0*PI
    if (eccanom .lt. 0.d0) eccanom = eccanom + 2.d0*PI

    ! Now get more precise solution using Newton Raphson method
    ! for those times when the Kepler equation is not yet solved
    ! to better than 1e-10
    ! (modification J. Wilms)

    mmm=mx
    if (mmm .lt. 0.d0) mmm = mmm + 2.d0*PI
    diff=eccanom-ecc*dsin(eccanom) - mmm

    if (dabs(diff) .gt. 1.d-10) then
        do
            ! e-e sine-m
            fe=eccanom-ecc*dsin(eccanom)-mmm
            ! f' = 1-e*cose
            fs=1.d0-ecc*dcos(eccanom)
            oldval=eccanom
            eccanom=oldval-fe/fs
            if (dabs(oldval-eccanom) .le. thresh) exit
        enddo 
        ! the following should be coded more intelligently ;-) 
        ! (similar to range reduction of mx...)
    endif 

    eccanom = mod(eccanom + PI, 2.d0*PI) - PI
    return
end subroutine keplereq
