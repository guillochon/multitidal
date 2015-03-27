      subroutine fss002(grid,ngrid,u,du)
      implicit none
      save

!this routine computes second order accurate first derivatives 
!on an arbitrary grid

!input:
!ngrid       = number of points in the grid, 
!grid(ngrid) = array of independent values 
!u(ngrid)    = function values at the grid points

!output:
!du(ngrid)   = first derivative values at the grid points

!declare the pass
      integer          ngrid
      double precision grid(ngrid),u(ngrid),du(ngrid)


!local variables
      integer          m,n
      parameter        (m=2, n=3)
      integer          i
      double precision coef(n)


!first point; one sided
      call fdcoef(m,n,grid(1),grid(1),coef)
      du(1) = coef(1)*u(1) + coef(2)*u(2) + coef(3)*u(3)


!middle of the grid; central differences
      do i=2,ngrid-1
       call fdcoef(m,n,grid(i),grid(i-1),coef)
       du(i) = coef(1)*u(i-1) + coef(2)*u(i) + coef(3)*u(i+1)
      enddo

!last point; one sided
      call fdcoef(m,n,grid(ngrid),grid(ngrid-2),coef)
      du(ngrid) = coef(1)*u(ngrid-2)+coef(2)*u(ngrid-1)+coef(3)*u(ngrid)
      return
      end

      subroutine fss004(grid,ngrid,u,du)
      implicit none
      save

!this routine computes fourth-order accurate first derivatives 
!on an arbitrary grid

!input:
!ngrid       = number of points in the grid, 
!grid(ngrid) = array of independent values 
!u(ngrid)    = function values at the grid points

!output:
!du(ngrid)   = first derivative values at the grid points

!declare the pass
      integer          ngrid
      double precision grid(ngrid),u(ngrid),du(ngrid)


!local variables
      integer          m,n
      parameter        (m=2, n=5)
      integer          i
      double precision coef(n)


!first point; one sided
      do i=1,2
       call fdcoef(m,n,grid(i),grid(i),coef)
       du(i) = dot_product(coef, u(i:i+4))
      enddo


!middle of the grid; central differences
      do i=3,ngrid-2
       call fdcoef(m,n,grid(i),grid(i-2),coef)
       du(i) = dot_product(coef, u(i-2:i+2))
      enddo

!last point; one sided
      do i=ngrid-1,ngrid
       call fdcoef(m,n,grid(i),grid(i-4),coef)
       du(i) = dot_product(coef, u(i-4:i))
      enddo

      return
      end

      subroutine fdcoef(mord,nord,x0,grid,coef)
      implicit none
      save

!this routine implements simple recursions for calculating the weights
!of finite difference formulas for any order of derivative and any order
!of accuracy on one-dimensional grids with arbitrary spacing.

!from bengt fornberg's article
!generation of finite difference formulas on arbitrary spaced grids. 
!math. comp., 51(184):699-706, 1988. 


!input:
!mord       = the order of the derivative 
!nord       = order of accuracy n
!x0         = point at which to evaluate the coefficients
!grid(nord) = array containing the grid

!output:
!coef(nord) = coefficients of the finite difference formula.


!declare the pass
      integer          mord,nord
      double precision x0,grid(nord),coef(nord)


!local variables
      integer          nu,nn,mm,nmmin,mmax,nmax
      parameter        (mmax=8, nmax=10)
      double precision weight(mmax,nmax,nmax),c1,c2,c3,c4,pfac


!zero the weight array
      do nu=1,nord
       do nn=1,nord
        do mm=1,mord
         weight(mm,nn,nu) = 0.0d0
        enddo
       enddo
      enddo

      weight(1,1,1) = 1.0d0
      c1            = 1.0d0
      nmmin         = min(nord,mord)
      do nn = 2,nord
       c2 = 1.0d0
       do nu=1,nn-1
        c3 = grid(nn) - grid(nu)
        c2 = c2 * c3
        c4 = 1.0d0/c3
        pfac = grid(nn) - x0
        weight(1,nn,nu) = c4 * ( pfac * weight(1,nn-1,nu) )
        do mm=2,nmmin
         weight(mm,nn,nu) = c4 * ( pfac * weight(mm,nn-1,nu) &
                           - dfloat(mm-1)*weight(mm-1,nn-1,nu) )
        enddo
       enddo
       pfac = (grid(nn-1) - x0)
       weight(1,nn,nn) = c1/c2 * (-pfac*weight(1,nn-1,nn-1))
       c4 = c1/c2
       do mm=2,nmmin
        weight(mm,nn,nn) = c4 * (dfloat(mm-1)*weight(mm-1,nn-1,nn-1) - &
                                 pfac*weight(mm,nn-1,nn-1))
       enddo
       c1 = c2
      enddo

!load the coefficients
      do nu = 1,nord
       coef(nu) = weight(mord,nord,nu)
      enddo
      return
      end
