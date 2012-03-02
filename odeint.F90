MODULE ode_path
    USE nrtype
    INTEGER(I4B) :: nok,nbad,kount
    REAL(DP), DIMENSION(:), POINTER :: xp
    REAL(DP), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path

SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,bsstep)
    USE nrtype
    USE nrutil, ONLY : nrerror,reallocate
    USE ode_path
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
    REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
    INTERFACE
        SUBROUTINE derivs(x,y,dydx)
        USE nrtype
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(IN) :: y
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
        END SUBROUTINE derivs
!BL
        SUBROUTINE bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
        USE nrtype
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
        REAL(DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
        REAL(DP), INTENT(INOUT) :: x
        REAL(DP), INTENT(IN) :: htry,eps
        REAL(DP), INTENT(OUT) :: hdid,hnext
        INTERFACE
            SUBROUTINE derivs(x,y,dydx)
            USE nrtype
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP), DIMENSION(:), INTENT(IN) :: y
            REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
            END SUBROUTINE derivs
        END INTERFACE
        END SUBROUTINE bsstep
    END INTERFACE
    REAL(DP), PARAMETER :: TINY=1.0e-30_dp
    INTEGER(I4B), PARAMETER :: MAXSTP=10000000
    INTEGER(I4B) :: nstp
    REAL(DP) :: h,hdid,hnext,x
    REAL(DP), DIMENSION(size(ystart)) :: dydx,y,yscal
    x=x1
    h=sign(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    y(:)=ystart(:)
    do nstp=1,MAXSTP
        call derivs(x,y,dydx)
        if (x .ne. x .or. any(y .ne. y) .or. any(dydx .ne. dydx)) then
            print *, x, y, dydx
            call Driver_abortFlash('NaN encountered!')
        endif
        yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
        if ((x+h-x2)*(x+h-x1) > 0.d0) h=x2-x
        call bsstep(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
        if (hdid == h) then
            nok=nok+1
        else
            nbad=nbad+1
        end if
        if ((x-x2)*(x2-x1) >= 0.d0) then
            ystart(:)=y(:)
            RETURN
        end if
        if (abs(hnext) < hmin)&
            call nrerror('stepsize smaller than minimum in odeint')
        h=hnext
    end do
    call nrerror('too many steps in odeint')
END SUBROUTINE odeint
