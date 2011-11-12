      subroutine polytr(order,mtot,rhoc,polyk,mu,mode,
     1                  x,y,yp,radius,rho,mass,prss,ebind,rhom,
     2                  ztemp,zbeta,exact,xsurf,ypsurf,n,iend,ipos)
      implicit none
      include 'const.dek'

c..this routine computes the physical properties of a polytropic star.
c..see chandrasekhars "stellar structure"


c..input :
c..order = order of the polytrope; gamma = 1 + 1/order
c..mode  = 1;  total mass and central density are given. 
c..            the polytropic constant is then computed.
c..      = 2;  total mass and polytropic constant are given. 
c..            the central density is then computed.
c..      = 3;  central density and polytropic constant given.
c..            the total mass is computed.
c..mu    = molecular weight
c..n     = dimension of the output arrays


c..output:
c..x = dimensionless radial coordinate
c..y = dimensionless density coordinate
c..yp = dimensionless derivative dx/dy
c..radius = radial coordinate (solar units)
c..density = mass density in g/cm**3
c..mass    = mass interior to the radial coordinate
c..prss    = pressure in erg/cm**3
c..ebind   = binding energy in erg/gr
c..rhom    = ratio of mean desnity to central density
c..ztemp   = temperature in K
c..zbeta   = ratio of gas to toal pressure
c..exact   = exact solution of n = 0, 1, or 5
c..xsurf   = x value where y is zero, dimensionless surface of the star
c..ypsurf  = dx/dy at the dimensionless suurface of the star
c..iend    = total number of solution points
c..ipos    = number of solution points with y > 0


c..declare the pass
      integer          n,iend,ipos,mode
      double precision order,mtot,rhoc,polyk,mu,
     1                 x(n),y(n),yp(n),radius(n),rho(n),mass(n),prss(n),
     2                 ebind(n),rhom(n),ztemp(n),zbeta(n),exact(n),
     3                 xsurf,ypsurf


c..locals variables
      external          lanemd,rkqc,tlane,cross
      logical           succes
      integer           i,k,mm,xdim,ydim,nok,nbad,iprint,iat
      parameter         (xdim=600, ydim=2)

      double precision  xrk(xdim),yrk(ydim,xdim),bc(ydim),stptry,
     1                  stpmin,stpmax,tol,odescal,lo,hi,zbrent,start,
     2                  sstop,xx,rscale,secday,f1,f2,f3


c..communication with routine lanemd and routine tlane
      double precision  oord,ppres,den,mmu,beta
      common   /poly1/  oord,ppres,den,mmu,beta



c..communication with routine cross
      integer           jpmax,jp
      parameter         (jpmax=20)
      double precision  xa(jpmax),ya(jpmax),ypa(jpmax),yps
      common   /zsurf/  xa,ya,ypa,yps,jp



c..various constants
      double precision  twoth,fpi,fpig
      parameter        (twoth  = 2.0d0/3.0d0,
     1                  fpi    = 4.0d0 * pi, 
     2                  fpig   = fpi * g)




c..initialize
       oord  = order
       mmu   = mu
       lo    = 1.0d2
       hi    = 1.0d9


c..set the initial conditions 
c..nominally x = 0 bc(1) = 0.0d0 and bc(2) = 1.0d0
c..start away from the x=0 singularity and move the boundary conditions as well

      start   = 1.0d-6
      call le_series(order,start,bc(2),bc(1))


c..set the numerical control parameters
c..integration for order=0 is so good that it makes ugly plots
c..so enforce a maximum step size 

      stptry  = 1.0d-8
      stpmin  = 1.0d-12
      stpmax  = 1.0e24
      if (order .eq. 0.0) stpmax = 0.05
      sstop   = 50.0d0
      tol     = 1.0d-8
      odescal = 1.0d0
      iprint  = 0


c..integrate to get the dimensionless solution

      call podeint(start,stptry,stpmin,sstop,bc,
     1            tol,stpmax,xdim,
     2            xrk,yrk,xdim,ydim,xdim,ydim,
     3            nok,nbad,k,odescal,iprint,
     4            lanemd,rkqc)


c..set the total number of points in the solution
      iend = min(k,n)


c..set the number of points with y > 0
      do i=1,iend
       if (yrk(2,i) .le. 0.0) goto 21
       ipos = i
      end do
 21   continue


c..transfer the solution in the output arrays
      do i=1,iend
       x(i)  = xrk(i)
       y(i)  = yrk(2,i)
       yp(i) = yrk(1,i)
      end do


c..find the zero crossing value
c..jp is the number of points to use in the polynomial fit

      if (ipos .lt. iend) then
       jp = 6
       iat = max(1,min(ipos - jp/2 + 1,iend - jp + 1))
       do mm=1,jp
        xa(mm)  = x(iat  + (mm-1))
        ya(mm)  = y(iat  + (mm-1))
        ypa(mm) = yp(iat + (mm-1))
       enddo
       xsurf  = zbrent(cross,xa(1),xa(jp),tol)
       ypsurf = yps
      else
       xsurf  = x(ipos)
       ypsurf = yp(ipos)
      end if





c..the rest of this routine figures out the physical, dimensional solution

c..for order < 1
c..serious problems of computing the physical properties for order < 1.
c..cases of real physical interest always have order > 1. 
c..see chandra pg 106, last part of section 8 for an amusing commentary.
c..in this case, only return the dimensionless solution x and y. 

      if (order .le. 1.0) then
       do i=1,ipos
        if (order .eq. 0.0) then
         exact(i)  = 1.0d0 - x(i)*x(i)/6.0d0
        else if (order .eq. 1.0) then
         exact(i)  = dsin(x(i))/x(i)
        end if
       enddo
       !return
      end if



c..for order > 1
      f1 = -fpi * xsurf * xsurf * ypsurf
      f2 = (order + 1.0d0) /fpig
      f3 = (3.0d0 - order)/(2.0d0 * order) 


c..total mass and central density given. compute the polytropic k
      if (mode .eq. 1) then
       polyk = 1.0d0/f2 * (mtot*msol/(f1*rhoc**f3))**twoth


c..total mass and polytropic k given. compute the central density
      else if (mode .eq. 2) then
       if (f3 .eq. 0.0) then 
        rhoc = 1.0d0
       else
        rhoc = ( mtot*msol/(f1 * (f2*polyk)**1.5d0 ) )**(1.0d0/f3)
       endif


c..central density and polytropic k given. compute the total mass.
      else if (mode .eq. 3) then
       mtot = f1/msol * (f2*polyk)**1.5d0 * rhoc**f3  

c..a bad mode
      else 
       stop 'unknown mode in routine polytr'
      end if



c..the constant radial scale factor
      rscale = dsqrt( f2*polyk * rhoc**(1.0d0/order - 1.0d0) )


c..start the output loop
      do i=1,ipos

c..radius
       radius(i) = rscale * x(i)

c..density; copy the value to common block
       rho(i) = rhoc * y(i)**order
       den    = rho(i)

c..the mass interior to the radius
       mass(i) = fpi * rhoc * rscale**3 * (-x(i)*x(i)*yp(i))

c..the total pressure
       prss(i)= polyk * rho(i)**(1.0d0 + 1.0d0/order)
       ppres  = prss(i)

c..the gravitational binding energy (take care of infinity at n=5)
       ebind(i)  = -1.0e30
       if (order .ne. 5.0  .and.  radius(i) .ne. 0.0) 
     1  ebind(i) = (3.0d0*mass(i)*mass(i)*g)/((order - 5.0d0)*radius(i))

c..the central to mean density ratio;if n=5 output the exact solution
       rhom(i)   = 1.0d0
       if (yp(i) .ne. 0.0)  rhom(i) = -x(i) / (3.0d0 * yp(i))
       if (order .eq. 5.0) exact(i) = 1.0d0/dsqrt(1.0d0+x(i)*x(i)/3.0d0)

c..the temperature is given by a root find on the gas + radiation pressures
c..be sure the solution is bracketed

       call zbrac(tlane,lo,hi,succes)
       if (.not.succes) then
        write(6,35) 'den=',den,' pres=',ppres,' lo=',lo,' hi=',hi
35      format(1x,4(a,1pe11.3))
        write(6,*) 'cannot bracket temperature, but carrying on'
        ztemp(i) = 0.0d0
        zbeta(i) = 0.0d0
        lo       = 1.0d2
        hi       = 1.0d9

c..reset the search limits for the next trip
       else
        ztemp(i)= zbrent(tlane,lo,hi,tol)
        lo      = 0.1d0 * ztemp(i)
        hi      = 2.0d0 * ztemp(i)
        zbeta(i) = beta
       end if

c..back for another point or bail
      enddo
      return
      end







      subroutine lanemd(x,y,dydx)
      implicit none

c..this routine evaluates the lane-emden equation of order n

c..declare
      double precision x,y(*),dydx(*),rz
      double complex   z


c..communication with other routines
      double precision  oord,ppres,den,mmu,beta
      common   /poly1/  oord,ppres,den,mmu,beta


c..allow for a negative y, although not physical, it exists
c..mathematically, and is useful for determining surface quantities

      z = dcmplx(y(2),0.0d0) 
      rz = abs(z**oord)

      dydx(1) = -2.0d0 * y(1)/x - rz
      dydx(2) = y(1)
      return
      end





      double precision function cross(xsurf)
      implicit none

c..this routine is used by a root finder to find the x coordinate 
c..where the y coordinate goes to zero. that is, find the surface of the star.

c..declare the pass
      double precision xsurf

c..local variables
      double precision ysurf,dy

c..communication with routine cross
      integer           jpmax,jp
      parameter         (jpmax=20)
      double precision  xa(jpmax),ya(jpmax),ypa(jpmax),yps
      common   /zsurf/  xa,ya,ypa,yps,jp


      call polint(xa,ya,jp,xsurf,ysurf,dy)
      call polint(xa,ypa,jp,xsurf,yps,dy)


c..we seek the zero of this function
      cross = ysurf
      return
      end 






      double precision function tlane(t)
      implicit none
      include 'const.dek'

c..given the total pressure, density, and mean molecular weight
c..this routine is used by a root finder to find the temperature 
c..at which the sum of the ideal gas and radiation pressures is 
c..equal to the total pressure given by the lane-emden solution.

c..declare
      double precision t

c..communication with other routines
      double precision  oord,ppres,den,mmu,beta
      common   /poly1/  oord,ppres,den,mmu,beta


c..ratio of gas to total pressure
      beta  = (avo*kerg*den*t)/(mmu*ppres)  

c..we seek the zero of this function
      !write(6,*) 'asol=',asol,' ppres=',ppres,' t=',t,' beta=',beta
      tlane = (asol*t*t*t*t)/(3.0d0*ppres) + beta - 1.0d0
      !write(6,*) 'tlane=',tlane
      return
      end 







      subroutine le_series(xn,z,f,fp)
      implicit none

c..given the order xn and a point z of the lane-emden equation,
c..this routine returns value and derivative of the function
c..through and an order 14 seriees expansion.

c..mathematica expansion from david reiss
c..http://www.scientificarts.com/laneemden/Links/laneemden_lnk_1.html

c..declare the pass
      double precision xn,z,f,fp

c..local variables
      double precision xn2,xn3,xn4,p6,p8,p10,p12,
     1                 z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12

      double precision c2,c4,c6,c8,c10,c12
      parameter       (c2  = 1.0d0/6.0d0,
     1                 c4  = 1.0d0/120.0d0,
     2                 c6  = 1.0d0/15120.0d0,
     3                 c8  = 1.0d0/3265920.0d0,
     4                 c10 = 1.0d0/1796256000.0d0,
     5                 c12 = 1.0d0/840647808000.0d0)



c..here we go
      xn2 = xn * xn
      xn3 = xn * xn2
      xn4 = xn * xn3

      p6  = xn*(8.0d0*xn - 5.0d0)
      p8  = xn*(122.0d0*xn2 - 183.0d0*xn + 70.0d0)
      p10 = xn*(5032.0d0*xn3 - 12642.0d0*xn2 + 10805.0d0*xn - 3150.0d0)
      p12 = xn*(183616.0d0*xn4 - 663166.0d0*xn3 + 915935.0d0*xn2 
     1                         - 574850.0d0*xn + 138600.0d0)

      z2  = z * z
      z3  = z * z2
      z4  = z * z3
      z5  = z * z4
      z6  = z * z5
      z7  = z * z6
      z8  = z * z7
      z9  = z * z8
      z10 = z * z9
      z11 = z * z10
      z12 = z * z11


c..the series expansion and its derivative with respect to z
      f  =  1.0d0 - c2*z2 + c4*z4*xn - c6*z6*p6
     1      + c8*z8*p8 - c10*z10*p10 + c12*z12*p12

      fp = -2.0d0*c2*z + 4.0d0*c4*z3*xn - 6.0d0*z5*p6
     1     + 8.0d0*c8*z7*p8 - 10.0d0*c10*z9*p10 + 12.0d0*c12*z11*p12

      return
      end








      subroutine podeint(start,stptry,stpmin,stopp,bc,
     1                   eps,stpmax,kmax, 
     2                   xrk,yrk,xphys,yphys,xlogi,ylogi,
     3                   nok,nbad,kount,odescal,iprint,
     4                   derivs,steper)  
      implicit none


c..basic ode integrator from numerical recipes.
c..special hooks added for polytropes.

c..declare  
      external         derivs,steper
      integer          nok,nbad,nmax,nstpmax,kmax,kount,xphys,
     1                 yphys,xlogi,ylogi,iprint,i,j,nstp
      parameter        (nmax = 20, nstpmax=10000)  
      double precision bc(yphys),stptry,stpmin,eps,stopp,start, 
     1                 yscal(nmax),y(nmax),dydx(nmax),  
     2                 stpmax,xrk(xphys),yrk(yphys,xphys),odescal, 
     3                 x,xsav,h,hdid,hnext,zero,one,tiny,ttiny
      parameter        (zero=0.0, one=1.0, tiny=1.0e-30, ttiny=1.0e-15)


c..here are the format statements for printouts as we integrate
100   format(1x,i4,1p10e12.4)



c..initialize   
      if (ylogi .gt. yphys) stop 'ylogi > yphys in routine odeint'
      if (yphys .gt. nmax)  stop 'yphys > nmax in routine odeint'
      x     = start   
      h     = sign(stptry,stopp-start) 
      nok   = 0 
      nbad  = 0
      kount = 0   


c..store the first step 
      do i=1,ylogi
       y(i) = bc(i)  
      enddo
      xsav = x


c..take at most nstpmax steps
      do nstp=1,nstpmax
       call derivs(x,y,dydx)


c..scaling vector used to monitor accuracy  
       do i=1,ylogi
        yscal(i) = abs(y(i)) + abs(h * dydx(i)) + tiny
       enddo


c..store intermediate results   
       if (kmax .gt. 0) then
        if ( kount .lt. (kmax-1) ) then  
         kount         = kount+1  
         xrk(kount)    = x   
         do i=1,ylogi 
          yrk(i,kount) = y(i)
         enddo
         if (iprint .eq. 1) then
          write(6,100) kount,xrk(kount),(yrk(j,kount), j=1,ylogi)
         end if
        end if
        xsav=x 
       end if



c..if the step can overshoot the stop point then cut it
       if ((x+h-stopp)*(x+h-start).gt.zero) h = stopp - x  


c..do an integration step
       call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext,derivs)   
       if (hdid.eq.h) then
        nok = nok+1   
       else 
        nbad = nbad+1 
       end if


c..bail if the solution gets "too" negative, we only need a negative y values
c..points to determine the zero crossing, the radius of the star

       if (y(2) .lt. -0.1) return


c..this is the normal exit point, save the final step   
       if (nstp.eq.nstpmax .or. (x-stopp)*(stopp-start) .ge. zero) then
        do i=1,ylogi  
         bc(i) = y(i) 
        enddo
        if (kmax.ne.0) then   
         kount         = kount+1  
         xrk(kount)    = x   
         do i=1,ylogi 
          yrk(i,kount) = y(i) 
         enddo
         if (iprint .eq. 1) then
           write(6,100) kount,xrk(kount),(yrk(j,kount), j=1,ylogi)
         end if
        end if
        return  
       end if


c..set the step size for the next iteration; stay above stpmin
       h = min(hnext,stpmax)
       if (abs(hnext).lt.stpmin) stop 'hnext < stpmin in odeint'


c..back for another iteration or death
      enddo
      write(6,*) '> than nstpmax steps required in odeint' 
      return
      end







      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)  
      implicit none

c..fifth order, step doubling, runge-kutta ode integrator with monitering of 
c..local truncation errors. input are the vector y of length n, which has a 
c..known the derivative dydx at the point x, the step size to be attempted 
c..htry, the required accuracy eps, and the vector yscal against which the 
c..error is to be scaled.  on output, y and x are replaced by their new values,
c..hdid is the step size that was actually accomplished, and hnext is the 
c..estimated next step size. derivs is a user supplied routine that computes 
c..the right hand side of the first order system of odes. plug into odeint.

c..declare  
      external         derivs   
      integer          n,nmax,i 
      parameter        (nmax = 2000)  
      double precision x,htry,eps,hdid,hnext,y(n),dydx(n),yscal(n), 
     1                 ytemp(nmax),ysav(nmax),dysav(nmax),fcor,safety,
     2                 errcon,pgrow,pshrnk,xsav,h,hh,errmax 
      parameter        (fcor=1.0d0/15.0d0, pgrow = -0.2d0, 
     1                  pshrnk = -0.25d0,  safety=0.9d0,  errcon=6.0e-4)


c..note errcon = (4/safety)**(1/pgrow)  
c..nmax is the maximum number of differential equations

c..save the initial values 
      h      = htry
      xsav   =  x
      do i=1,n   
       ysav(i)  = y(i)
       dysav(i) = dydx(i)
      enddo


c..take two half steps  
1     hh = 0.5d0*h  
      call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)   
      x  = xsav + hh 
      call derivs(x,ytemp,dydx) 
      call rk4(ytemp,dydx,n,x,hh,y,derivs)  
      x  = xsav + h  
      if (x .eq. xsav) stop 'stepsize not significant in rkqc' 


c..now take the large step  
      call rk4(ysav,dysav,n,xsav,h,ytemp,derivs)


c..ytemp is the error estimate  
      errmax = 0.0d0
      do i=1,n   
       ytemp(i) = y(i) - ytemp(i)  
       errmax   = max(errmax,abs(ytemp(i)/yscal(i)))   
      enddo
      errmax     = errmax/eps 


c..truncation error too big, reduce the step size and try again 
      if (errmax .gt. 1.0) then
       h = safety * h * (errmax**pshrnk) 
       go to  1  


c..truncation within limits, compute the size of the next step  
      else  
       hdid = h  
       if (errmax.gt.errcon) then
        hnext = safety * h * (errmax**pgrow)
       else
        hnext = 4.0d0 * h
       end if   
      end if 


c..mop up the fifth order truncation error  
      do i=1,n   
       y(i) = y(i) + ytemp(i)*fcor 
      enddo
      return
      end





       subroutine rk4(y,dydx,n,x,h,yout,derivs) 
       implicit none

c..given values for the variables y(1:n) and their derivatives dydx(1:n) known
c..at x, use the fourth order runge-kutta method to advance the solution over
c..an interval h and return the incremented variables in yout(1:n) (which need
c..not be a distinct array from y). one supplies the routine derivs which 
c..evaluates the right hand side of the ode's.   

c..declare  
       external          derivs 
       integer           n,nmax,i   
       parameter         (nmax = 2000)  
       double precision  x,h,y(n),dydx(n),yout(n),
     1                   yt(nmax),dyt(nmax),dym(nmax),
     2                   hh,h6,xh   


c..initialize the step sizes and weightings 
       hh = h*0.5d0
       h6 = h/6.0d0 
       xh = x + hh  

c..the first step   
       do i=1,n  
        yt(i) = y(i) + hh*dydx(i)  
       enddo


c..the second step  
       call derivs(xh,yt,dyt)   
       do i=1,n  
        yt(i) = y(i) + hh*dyt(i)   
       enddo


c..the third step   
       call derivs(xh,yt,dym)   
       do i=1,n  
        yt(i)  = y(i) + h*dym(i)
        dym(i) = dyt(i) + dym(i)   
       enddo


c..the fourth step and accumulate the increments with the proper weights
       call derivs(x+h,yt,dyt)  
       do i=1,n  
        yout(i) = y(i) + h6*(dydx(i) +dyt(i) + 2.0d0*dym(i)) 
       enddo
       return   
       end  





      subroutine zbrac(func,x1,x2,succes)   
      implicit none

c..given a function func and an initial guessed range x1 to x2, the routine 
c..expands the range geometrically until a root is bracketed by the returned
c..values x1 and x2 (in which case succes returns as .true.) or until the range
c..becomes unacceptably large (in which case succes returns as .false.).    
c..success  guaranteed for a function which has opposite sign for sufficiently
c..large and small arguments.   

c..declare 
      external          func
      logical           succes  
      integer           ntry,j  
      parameter         (ntry=500)
      double precision  func,x1,x2,factor,f1,f2 
      parameter         (factor = 1.6d0)
c.. 
      if (x1 .eq. x2) stop ' x1 = x2 in routine zbrac'
      f1 = func(x1) 
      f2 = func(x2) 
      succes = .true.   
      do j=1,ntry
       if (f1*f2 .lt. 0.0) return   
       if (abs(f1) .lt. abs(f2)) then   
        x1 = x1 + factor * (x1-x2)
        f1 = func(x1)   
       else 
        x2 = x2 + factor * (x2-x1)
        f2 = func(x2)   
       end if   
      enddo
      succes = .false.  
      return
      end   





      double precision function zbrent(func,x1,x2,tol)  
      implicit none

c..using brent's method this routine finds the root of a function func between
c..the limits x1 and x2. the root is when accuracy is less than tol 

c..declare
      external          func
      integer           itmax,iter  
      parameter         (itmax=100)  
      double precision  func,x1,x2,tol,a,b,c,d,e,fa,
     1                  fb,fc,xm,tol1,p,q,r,s,eps   
      parameter         (eps=3.0d-15)  

c..note: eps the the machine floating point precision


c..initialize
      a  = x1
      b  = x2
      fa = func(a)  
      fb = func(b)  
      if ( (fa .gt. 0.0  .and. fb .gt. 0.0)  .or.
     1     (fa .lt. 0.0  .and. fb .lt. 0.0)       ) then
       write(6,100) x1,fa,x2,fb
100    format(1x,' x1=',1pe11.3,' f(x1)=',1pe11.3,/,
     1        1x,' x2=',1pe11.3,' f(x2)=',1pe11.3)
       stop 'root not bracketed in routine zbrent'   
      end if
      c  = b
      fc = fb   


c..start the iteration loop
      do iter =1,itmax   

c..rename a,b,c and adjusting bound interval d  
       if ( (fb .gt. 0.0  .and. fc .gt. 0.0)  .or.
     1      (fb .lt. 0.0  .and. fc .lt. 0.0)      ) then
        c  = a   
        fc = fa 
        d  = b-a 
        e  = d   
       end if   
       if (abs(fc) .lt. abs(fb)) then   
        a  = b   
        b  = c   
        c  = a   
        fa = fb 
        fb = fc 
        fc = fa 
       end if   
       tol1 = 2.0d0 * eps * abs(b) + 0.5d0 * tol
       xm   = 0.5d0 * (c-b) 


c..convergence check
       if (abs(xm) .le. tol1 .or. fb .eq. 0.0) then 
        zbrent = b  
        return  
       end if   


c..attempt quadratic interpolation  
       if (abs(e) .ge. tol1 .and. abs(fa) .gt. abs(fb)) then
        s = fb/fa   
        if (a .eq. c) then  
         p = 2.0d0 * xm * s   
         q = 1.0d0 - s 
        else
         q = fa/fc  
         r = fb/fc  
         p = s * (2.0d0 * xm * q *(q-r) - (b-a)*(r - 1.0d0))  
         q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
        end if  


c..check if in bounds   
        if (p .gt. 0.0) q = -q   
        p = abs(p)  

c..accept interpolation 
        if (2.0d0*p .lt. min(3.0d0*xm*q - abs(tol1*q),abs(e*q))) then   
         e = d  
         d = p/q

c..or bisection 
        else
         d = xm 
         e = d  
        end if  

c..bounds decreasing to slowly use bisection
       else 
        d = xm  
        e = d   
       end if   

c..move best guess to a 
       a  = b
       fa = fb  
       if (abs(d) .gt. tol1) then   
        b = b + d   
       else 
        b = b + sign(tol1,xm)   
       end if   
       fb = func(b) 

      enddo
      stop 'too many iterations in routine zbrent'  
      end   





      subroutine polint(xa,ya,n,x,y,dy)
      implicit none

c..given arrays xa and ya of length n and a value x, this routine returns a 
c..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
c..such that ya = p(xa) ya then the returned value is y = p(x) 


c..declare the pass
      integer          n
      double precision xa(n),ya(n),x,y,dy


c..local variables
      integer          nmax,ns,i,m
      parameter        (nmax=10)
      double precision c(nmax),d(nmax),dif,dift,ho,hp,w,den


c..find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

c..first guess for y
      y = ya(ns)

c..for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

c..after each column is completed, decide which correction c or d, to add
c..to the accumulating value of y, that is, which path to take in the table
c..by forking up or down. ns is updated as we go to keep track of where we
c..are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end


