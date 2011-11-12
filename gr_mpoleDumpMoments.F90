!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleDumpMoments
!!
!! NAME
!!  gr_mpoleDumpMoments
!!
!! SYNOPSIS
!!  gr_mpoleDumpMoments()
!!
!! DESCRIPTION
!!
!!  Utility routine to output the Moment array to a text file.  Writes out in 
!!    one of two HARDWIRED formats:
!!           OUTPUT = 1 sorts by Even/Odd/Inner/Outer
!!           OUTPUT = 2 sorts by q  (Cal's request, current default)
!!  The information is written out to a file named basenm_momentDump.txt, where
!!    basenm is the runtime parameter.  The file is appended at each time.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!  gr_mpoleDumpMoment [BOOLEAN] -- set to true to have this routine called
!!  
!! NOTES
!!
!!  At the end of each time, the distinctive phrase "Chakka Khan" is inserted
!!    to make post-processing easier.  Ask Calhoun why.
!!
!!***


subroutine gr_mpoleDumpMoments()
  use Grid_data, ONLY : gr_myPE
  use gr_mpoleData, ONLY : Moment, mpole_lmax, mpole_mmax, qmax, rad
  use Driver_interface, ONLY : Driver_getSimTime
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
   
  ! arguments none

  ! local variables
  integer ::  posBlank
  real    :: time
  integer  :: fileUnit, k,p,q,l,m
  character (len=MAX_STRING_LENGTH), save :: string, baseName
  character (len=MAX_STRING_LENGTH), save :: fileName
  integer     :: nStep, nStepStart
  logical  :: firstCall = .true.
  character(len=5)  :: kString(2), pString(2)
  integer    :: OUTPUT = 2     ! HARDWIRED -- change if different format required
  real       :: radialDist

  kString(1) = "ODD  "
  kString(2) = "EVEN "
  pString(1) = "INNER"
  pString(2) = "OUTER"
  fileUnit = 236
  

  ! only write out on first processor
  if (gr_myPE /= MASTER_PE) return

  ! Open up the file
  if (firstCall) then
     call RuntimeParameters_get("basenm",baseName)
     posBlank = index(baseName,' ')
     fileName = baseName(:posBlank-1) // 'dumpMoments.txt'

     open (fileUnit, file=fileName)
     write(fileUnit,90)

     firstCall = .false.
  else
     open (fileUnit, file=fileName, position='APPEND')
  end if
  
  call Driver_getSimTime(time)
  write(fileUnit,95)time, qmax,mpole_lmax,mpole_mmax

   if (OUTPUT .EQ. 1) then

      do k = 2, 1, -1
         do p = 1, 2, 1
            write(fileUnit,*)   ! blank line
            write(fileUnit,100)kString(k),pString(p)
            write(fileUnit,110)
            do q = 0, qmax
               write(fileUnit,200)q,((Moment(q,k,p,l,m),l=0,mpole_lmax),m=0,mpole_mmax)
            end do
         end do
      end do
   else if (OUTPUT .EQ. 2) then

      write(fileUnit,500)
      do q = 0, qmax
         radialDist = rad(q)
         do p = 1, 2, 1
            do l = 0, mpole_lmax
               do m = 0, mpole_mmax
                  do k = 1, 2, 1
                     write(fileUnit,550)q,radialDist,pString(p),l,m,kString(k),Moment(q,k,p,l,m)
                  end do
               end do
            end do
         end do
      end do

   endif

   write (fileUnit,96)
   write (fileUnit,*)  !! blank line to end
   close (fileUnit)
  
90  format("Moment dump")
95  format(T5,"time = ",G12.4,5X,"qmax = ",I5,5X,"mpole_lmax=",I5,5X,"mpole_mmax=",I5)
96  format(" Chakka Khan Chakka Khan I feel for you this is the end of this timestep")

100 format("k= ",A5," p= ",A5)
110 format(T5,"q",T15," Moment(q,k,p,l,m), l=0,mpole_lmax, m=0,mpole_mmax)")
200 format(I5,((T10,9G10.4)))

500 format(T3,"q",T5,"rad",T20,"p",T25,"l",T30,"m",T35,"k",T42,"Moment")
550 format(I4,T5,G10.4,1X,A5,1X,I4,1X,I4,1X,A5,1X,G12.6)
  
  

end subroutine gr_mpoleDumpMoments


