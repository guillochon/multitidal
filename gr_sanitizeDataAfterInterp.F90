!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_sanitizeDataAfterInterp
!!
!! NAME
!!
!!  gr_sanitizeDataAfterInterp
!!
!!
!! SYNOPSIS
!!
!!  gr_sanitizeDataAfterInterp(integer(in) :: blkList(count),
!!                             integer(in) :: count,
!!                             integer(in) :: layers(MDIM))
!!
!!
!!
!! DESCRIPTION
!!
!!  Given a list of blocks of data, loop over all of the blocks and
!!  check whether solution data in certain variables lie in a reasonable
!!  range of values.
!!
!!  Energies (ENER_VAR and EINT_VAR) are expected to be .ge. gr_smalle,
!!  and the density (DENS_VAR) is expected to be .ge. gr_smallrho,
!!  where gr_smalle and gr_smallrho are lower bounds coming from the
!!  runtime parameters smalle and smlrho, respectively.
!!
!!  For data that do satisfy these expectations, warning messages are
!!  generated, but the offending data is not modified.
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   layers - number of guardcell layers to be included in the check 
!!
!! NOTES
!!
!!  The checks are based on gr_conserveToPrimitive, which is called to
!!  convert solution data back from conserved form when using the old
!!  (convertToConsvdForMeshCalls) way of ensuring that the mesh
!!  handles data interpolation in conserved form.
!!
!!  This is meant to be called where gr_conserveToPrimitive used to be
!!  called when using the new (convertToConsvdInMeshInterp) way of
!!  ensuring that the mesh handles data interpolation in conserved
!!  form.
!!
!! SEE ALSO
!!
!!  gr_conserveToPrimitive
!!
!!
!! BUGS
!!
!!  This routine accesses the global variable storage 
!!  array unk directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

#define DEBUG_CONSCONV

subroutine gr_sanitizeDataAfterInterp(blkList,count, info, layers)

  use Grid_data, ONLY : gr_smallrho,gr_smalle
  use Grid_interface, ONLY : Grid_getMyPE
  use Logfile_interface, ONLY : Logfile_stamp
  use physicaldata, ONLY:unk, gcell_on_cc
  use tree, ONLY:nodetype
  use paramesh_dimensions, ONLY: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd, kl_bndi

  implicit none
#include "constants.h"
#undef REAL_FORMAT
#define REAL_FORMAT "(1PG23.16)"
#include "Flash.h"

  integer,intent(IN) :: count
  integer, dimension(count), intent(IN) :: blkList
  character(len=*), intent(IN) :: info
  integer,dimension(MDIM), intent(IN):: layers
  integer :: n, block

  integer :: myPe, i,j
  integer :: iskip, jskip, kskip
  integer :: il,iu,jl,ju,kl,ku
  character(len=32), dimension(4,2) :: block_buff
  character(len=32)                 :: number_to_str
  
#ifndef HIDE_SANITIZE

111 format (a,a,a1,(1x,a18,'=',a),(1x,a2,'=',a5),(1x,a5,'=',a),(1x,a4,'=',a))
112 format (i3,1x,16(1x,1PG8.2))

  iskip = NGUARD - layers(IAXIS)
  jskip = (NGUARD - layers(JAXIS)) * K2D
  kskip = (NGUARD - layers(KAXIS)) * K3D
  il = il_bnd + iskip
  iu = iu_bnd - iskip
  jl = jl_bnd + jskip
  ju = ju_bnd - jskip
  kl = kl_bnd + kskip
  ku = ku_bnd - kskip

  call Grid_getMyPE(myPe)

  do n = 1,count
     block=blkList(n)

#ifdef DENS_VAR
     if (gcell_on_cc(DENS_VAR)) then
        ! small limits -- in case the interpolants are not monotonic
        if (any(unk(DENS_VAR,il:iu,jl:ju,kl:ku,block) .LT. gr_smallrho)) then
           write (block_buff(1,1), '(a18)') 'min. unk(DENS_VAR)'
           !        write (number_to_str, '('//REAL_FORMAT//',a1)') minval(unk(DENS_VAR,il:iu,jl:ju,kl:ku,block)), ','
           write (number_to_str, '(G30.22)') minval(unk(DENS_VAR,il:iu,jl:ju,kl:ku,block))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') myPe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') block
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'type'
           write (number_to_str, '(i7)') nodetype(block)
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

           call Logfile_stamp(myPe, block_buff, 4, 2, 'WARNING '//info)
           print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)
#ifdef DEBUG_CONSCONV
           ! For 2D, this prints a slice at the lowest k index that is interior - KW
           do j=ju_bnd,jl_bnd,-1
              print 112, j, (unk(DENS_VAR,i,j,kl_bndi,block), i=il_bnd,iu_bnd)
           end do
#endif
        end if
     end if
#endif

#ifdef ENER_VAR               
     if (gcell_on_cc(ENER_VAR)) then
        ! energy
        if (any(unk(ENER_VAR,il:iu,jl:ju,kl:ku,block) .LT. gr_smalle*0.999999999)) then
           write (block_buff(1,1), '(a)') 'min. unk(ENER_VAR)'
           write (number_to_str, '('//REAL_FORMAT//')') minval(unk(ENER_VAR,il:iu,jl:ju,kl:ku,block))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') myPe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') block
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'type'
           write (number_to_str, '(i7)') nodetype(block)
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

           call Logfile_stamp(myPe, block_buff, 4, 2, 'WARNING '//info)
           print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)
#ifdef DEBUG_CONSCONV
           do j=ju_bnd,jl_bnd,-1
              print 112, j, (unk(ENER_VAR,i,j,kl_bndi,block), i=il_bnd,iu_bnd) 
           end do
#endif
        end if
     end if
#endif
#ifdef EINT_VAR
     if (gcell_on_cc(EINT_VAR)) then
        if (any(unk(EINT_VAR,il:iu,jl:ju,kl:ku,block) .LT. gr_smalle*0.999999999)) then
           write (block_buff(1,1), '(a)') 'min. unk(EINT_VAR)'
           write (number_to_str, '('//REAL_FORMAT//')') minval(unk(EINT_VAR,il:iu,jl:ju,kl:ku,block))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') myPe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') block
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'type'
           write (number_to_str, '(i7)') nodetype(block)
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

           call Logfile_stamp(myPe, block_buff, 4, 2, 'WARNING '//info)
           print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)
#ifdef DEBUG_CONSCONV
           do j=ju_bnd,jl_bnd,-1
              print 112, j, (unk(EINT_VAR,i,j,kl_bndi,block), i=il_bnd,iu_bnd) 
           end do
#endif
        end if
     end if
#endif

  end do
#endif

  return 
end subroutine gr_sanitizeDataAfterInterp
        
