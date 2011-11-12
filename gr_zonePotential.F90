!!****if* source/Grid/GridSolvers/Multipole/gr_zonePotential
!!
!! NAME
!!
!!  gr_zonePotential
!!
!!
!! SYNOPSIS
!!
!!  gr_zonePotential(real(in) :: xprime, 
!!                  real(IN) :: yprime,
!!                  real(IN) :: zprime,
!!                  real(OUT) :: potential)
!!
!!
!!
!! DESCRIPTION
!!
!!  A routine to calculate the value of the potential at one location using 
!!  the moments of the interior and exterior mass distributions.
!!
!!
!! ARGUMENTS
!!
!!     xprime  -   X Coordinates of the current zone, relative to
!!                                 the center of mass
!!     yprime  -   Y Coordinates of the current zone, relative to
!!                                 the center of mass
!!     zprime  -   Z Coordinates of the current zone, relative to
!!
!!   potential  -      The (estimated) value of the potential in the current 
!!                  zone (output)
!!
!!
!!
!! NOTES
!!
!!  On exit, the potential at the specified location (actually, the estimate 
!!  based on the multipole expansion) is contained in potential.
!!
!!  The associated Legendre polynomial values are calculated using a recurrence 
!!  relation given in Press et al. (2d ed) pp 247-8.  The loops over l and m 
!!  are arranged so as to most efficiently utilize this relation (ie. first 
!!  loop over m, then over l).  The cosines and sines of multiples of the 
!!  azimuthal angle are computed using the recurrence relations for the cosine 
!!  and sine functions (see Press et al. p173).
!!
!!  To obtain the zone-averaged potential, mpole_potential calls this routine 
!!  for several different points within a zone, then finds the average of the 
!!  returned values.
!!
!!***

subroutine gr_zonePotential (xprime, yprime, zprime, potential)

  use Grid_data, ONLY : gr_myPE
  use Logfile_interface, ONLY : Logfile_stamp

  use gr_mpoleData, ONLY : mpole_geometry, G_1DSPHERICAL,&
                         G_2DCYLINDRICAL,G_3DCARTESIAN,G_3DAXISYMMETRIC,&
                         mpole_mmax,mpole_lmax,costable,sintable,Moment,&
                         Legk1, Legk2, Even, Odd, Inner, Outer,&
                         rpower, rprinv, qmax, dsinvarr, rad, mpole_subsample
  use Driver_interface, ONLY: Driver_abortFlash
  use Grid_interface, ONLY: Grid_getMinCellSize
  
  implicit none


  !            Local variables:
  !
  !               rprime          The length of the current position vector
  !               qprime          The index q (for Moment()) at which the current
  !                                 cell falls
  !               lprime          The projection of rprime into the xy-plane
  !               mtemp1-7         Temporary variables
  !               l,m,h           Temporary indices
  !               trigalpha,
  !               trigbeta        Trig coefficients for sine/cosine calculations
  !               costheta        The cosine of the current polar angle
  !               Legnd1          The (m,m)th Legendre polynomial, evaluated at
  !                                 costheta
  !               Legnd2          The (m+1,m)th Legendre polynomial
  !               Legndr          The (l,m)th Legendre polynomial, with l>m+1
  
  
  real, intent(in) :: xprime, yprime, zprime
  real, intent(out) ::  potential
  
  integer :: l, m, lm1, mp1, mm1, h, qprime, qprmp1, qprmm1, i
  real    :: mtemp1, mtemp2, mtemp3, mtemp4, mtemp5, mtemp7
  real    :: trigalpha, trigbeta, costheta
  real    :: Legnd1, Legnd2, Legndr, lprime, rprime
  real    :: rscaled, f, f1, min_cell_size
  character(len=256) :: str_buffer
  
  !==========================================================================
  
  !                       Clear the potential variable.
  
  potential = 0.
  
  !                       Calculate the length of the position vector r'
  !                       and of its projection into the xy-plane (l').
  
  call Grid_getMinCellSize(min_cell_size)
  lprime = xprime**2 + yprime**2
  lprime = max((min_cell_size/mpole_subSample)**2.d0,lprime) !Added by JFG
  rprime = lprime + zprime**2
  
  lprime = sqrt(lprime)
  rprime = sqrt(rprime)
  
  qprime = -1
  do i = 0, qmax
      if (rad(i) .gt. rprime) then
          qprime = i
          f = (rprime - rad(i))*dsinvarr(i)
          f1 = 1.d0 - f
          exit
      endif
  enddo

  if (qprime .eq. -1) then
      f = 1.d0
      f1 = 0.d0
      qprime = qmax
  endif

  if (qprime .eq. qmax) then
      qprmp1 = qprime
  else
      qprmp1 = qprime + 1
  endif
  qprmm1 = qprime - 1
  
  if (qprime .gt. qmax .or. qprime .eq. -1 .or. qprime .eq. 0 .or. (xprime .eq. 0 .and. yprime .eq. 0) .or. &
      rprime .eq. 0.d0 .or. lprime .eq. 0.d0) then
     call Driver_abortFlash('arrrrrrrg!')
     write (str_buffer,*) 'gr_zonePotential: desired q = ', qprime, & 
          ' is invalid, or lprime = 0. qmax = ', qmax, ' rprime = ', rprime, ' lprime = ', lprime, &
          ' xprime = ', xprime, ' yprime = ', yprime, ' zprime = ', zprime
     call Logfile_stamp(gr_myPE, str_buffer, 'warning')
  endif
  
  !                       Compute powers of r'.
  
  mtemp1 = 1. / lprime
  mtemp2 = 1. / rprime
  
  rpower(0) = 1.
  rprinv(0) = mtemp2
  
  ! Compute some factors for higher-order moments.
  
  if (mpole_lmax > 1) then
     
     do l = 1, mpole_lmax
        lm1 = l - 1
        rpower(l) = rpower(lm1) * rprime
        rprinv(l) = rprinv(lm1) * mtemp2
     enddo
     
     ! Compute table of cosines and sines for multiples of the azimuthal angle 
     ! phi'; also compute cos(theta').
     
     costable(0) = 1.
     sintable(0) = 0.
     
     trigalpha = 1. - (xprime * mtemp1)
     trigbeta = yprime * mtemp1
     
     do m = 1, mpole_mmax
        mm1 = m - 1
        costable(m) = costable(mm1) - (trigalpha*costable(mm1) + & 
             trigbeta*sintable(mm1))
        sintable(m) = sintable(mm1) - (trigalpha*sintable(mm1) - & 
             trigbeta*costable(mm1))
     enddo
     
     costheta = zprime * mtemp2
     mtemp7 = sqrt((1.-costheta)*(1.+costheta))

  endif
  
  !                       Compute the contributions of the even and odd moments
  !                       to the potential.
  
  !                       Do (l,0) moments for l = 0...mpole_lmax.
  
  !                               Do (0,0).

  Legnd1 = 1.
  potential = potential + & 
       (f1*Moment(qprmm1,Even,Inner,0,0) + & 
       f *Moment(qprime,Even,Inner,0,0))*rprinv(0) + & 
       (f1*Moment(qprime,Even,Outer,0,0) + & 
       f *Moment(qprmp1,Even,Outer,0,0))*rpower(0)
  
  ! For 1D spherically symmetric problems only l = m = 0 contributes.

  if (mpole_geometry == G_1DSPHERICAL) return
  
  ! For 2D/3D problems, continue with l > 0.
  
  !                               Do (1,0).
  
  if (mpole_lmax >= 1) then
     Legnd2 = costheta * Legnd1
     potential = potential + & 
          (f1*Moment(qprmm1,Even,Inner,1,0) + & 
          f *Moment(qprime,Even,Inner,1,0))*rprinv(1)*Legnd2 + & 
          (f1*Moment(qprime,Even,Outer,1,0) + & 
          f *Moment(qprmp1,Even,Outer,1,0))*rpower(1)*Legnd2
  endif

  !                               Do (2,0) ... (mpole_lmax,0).
  
  do l = 2, mpole_lmax
     Legndr = costheta * Legk1(l,0) * Legnd2 - Legk2(l,0) * Legnd1
     Legnd1 = Legnd2
     Legnd2 = Legndr
     potential = potential + & 
          (f1*Moment(qprmm1,Even,Inner,l,0) + & 
          f *Moment(qprime,Even,Inner,l,0))*rprinv(l)*Legndr + & 
          (f1*Moment(qprime,Even,Outer,l,0) + & 
          f *Moment(qprmp1,Even,Outer,l,0))*rpower(l)*Legndr
  enddo
  
  ! For 2D axisymmetric problems only m = 0 contributes.
  
  if (mpole_geometry == G_2DCYLINDRICAL) return
  
  ! Also for 3D axisymmetric
  
  if (mpole_geometry == G_3DAXISYMMETRIC) return

  ! For 3D problems, continue with m > 0.
  
  !                       Do (l,m) moments for m > 0.
  
  do m = 1, mpole_lmax
     
     mp1 = m + 1
     
     !                               Do (m,m).
     
     Legnd1 = 1.
     mtemp5 = 1.
     do h = 1, m
        Legnd1 = -Legnd1 * mtemp5 * mtemp7
        mtemp5 = mtemp5 + 2.
     enddo
     mtemp3 = rprinv(m) * Legnd1
     mtemp4 = rpower(m) * Legnd1
  if (Moment(qprmm1,Even,Inner,m,m) .ne. Moment(qprmm1,Even,Inner,m,m) .or. &
      Moment(qprime,Even,Inner,m,m) .ne. Moment(qprime,Even,Inner,m,m) .or. &
      Moment(qprime,Even,Outer,m,m) .ne. Moment(qprime,Even,Outer,m,m) .or. &
      Moment(qprmp1,Even,Outer,m,m) .ne. Moment(qprmp1,Even,Outer,m,m) .or. &
      Moment(qprmm1,Odd,Inner,m,m) .ne. Moment(qprmm1,Odd,Inner,m,m) .or. &
      Moment(qprime,Odd,Inner,m,m) .ne. Moment(qprime,Odd,Inner,m,m) .or. &
      Moment(qprime,Odd,Outer,m,m) .ne. Moment(qprime,Odd,Outer,m,m) .or. &
      Moment(qprmp1,Odd,Outer,m,m) .ne. Moment(qprmp1,Odd,Outer,m,m)) then
  endif
     potential = potential + & 
          ((f1*Moment(qprmm1,Even,Inner,m,m) + & 
          f *Moment(qprime,Even,Inner,m,m))*mtemp3 + & 
          (f1*Moment(qprime,Even,Outer,m,m) + & 
          f *Moment(qprmp1,Even,Outer,m,m))*mtemp4)*costable(m) + & 
          ((f1*Moment(qprmm1,Odd,Inner,m,m) + & 
          f *Moment(qprime,Odd,Inner,m,m))*mtemp3 + & 
          (f1*Moment(qprime,Odd,Outer,m,m) + & 
          f *Moment(qprmp1,Odd,Outer,m,m))*mtemp4)*sintable(m)
     
     !                               Do (m+1,m).
     
     if (mp1 .le. mpole_lmax) then
        Legnd2 = costheta * (2*m+1) * Legnd1
        mtemp3 = rprinv(mp1) * Legnd2
        mtemp4 = rpower(mp1) * Legnd2
        potential = potential + & 
             ((f1*Moment(qprmm1,Even,Inner,mp1,m) + & 
             f *Moment(qprime,Even,Inner,mp1,m))*mtemp3 + & 
             (f1*Moment(qprime,Even,Outer,mp1,m) + & 
             f *Moment(qprmp1,Even,Outer,mp1,m))*mtemp4)*costable(m) + & 
             ((f1*Moment(qprmm1,Odd,Inner,mp1,m) + & 
             f *Moment(qprime,Odd,Inner,mp1,m))*mtemp3 + & 
             (f1*Moment(qprime,Odd,Outer,mp1,m) + & 
             f *Moment(qprmp1,Odd,Outer,mp1,m))*mtemp4)*sintable(m)
     endif
     
     !                               Do (m+2,m) ... (mpole_lmax,m).
     
     do l = m+2, mpole_lmax
        Legndr = costheta * Legk1(l,m) * Legnd2 - Legk2(l,m) * Legnd1
        Legnd1 = Legnd2
        Legnd2 = Legndr
        mtemp3 = rprinv(l) * Legndr
        mtemp4 = rpower(l) * Legndr
        potential = potential + & 
             ((f1*Moment(qprmm1,Even,Inner,l,m) + & 
             f *Moment(qprime,Even,Inner,l,m))*mtemp3 + & 
             (f1*Moment(qprime,Even,Outer,l,m) + & 
             f *Moment(qprmp1,Even,Outer,l,m))*mtemp4)*costable(m) + & 
             ((f1*Moment(qprmm1,Odd,Inner,l,m) + & 
             f *Moment(qprime,Odd,Inner,l,m))*mtemp3 + & 
             (f1*Moment(qprime,Odd,Outer,l,m) + & 
             f *Moment(qprmp1,Odd,Outer,l,m))*mtemp4)*sintable(m)
     enddo
     
  enddo
  
  !===============================================================================
  
  if (potential .ne. potential) then
      print *, qprime, qprmm1, qprmp1, qmax, rprime, lprime, f, f1, rad(i), dsinvarr(i)
      call Driver_abortFlash('potential is NaN in zonePotential!')
  endif

  return
end subroutine gr_zonePotential
