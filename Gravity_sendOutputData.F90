!!****f* source/Gravity/Gravity_sendOutputData
!!
!! NAME
!!  Gravity_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Gravity_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Gravity unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Gravity_sendOutputData()

    use IO_interface, ONLY : IO_setScalar
    use gr_mpoleData, ONLY: Mtot, X_centerofmass, Y_centerofmass, Z_centerofmass 
    use gr_isoMpoleData, ONLY: isoMtot => Mtot, Xcm, Ycm, Zcm
    use Gravity_data
    implicit none

#include "Flash.h"

#ifdef FLASH_MPOLE
    call IO_setScalar("Mtot", Mtot)
    call IO_setScalar("X_centerofmass", X_centerofmass)
    call IO_setScalar("Y_centerofmass", Y_centerofmass)
    call IO_setScalar("Z_centerofmass", Z_centerofmass)
#else
    call IO_setScalar("Mtot", isoMtot)
    call IO_setScalar("X_centerofmass", Xcm)
    call IO_setScalar("Y_centerofmass", Ycm)
    call IO_setScalar("Z_centerofmass", Zcm)
#endif
    call IO_setScalar("ptxpos", grv_ptvec(1))
    call IO_setScalar("ptypos", grv_ptvec(2))
    call IO_setScalar("ptzpos", grv_ptvec(3))
    call IO_setScalar("ptxvel", grv_ptvec(4))
    call IO_setScalar("ptyvel", grv_ptvec(5))
    call IO_setScalar("ptzvel", grv_ptvec(6))
    call IO_setScalar("obxpos", grv_obvec(1))
    call IO_setScalar("obypos", grv_obvec(2))
    call IO_setScalar("obzpos", grv_obvec(3))
    call IO_setScalar("obxvel", grv_obvec(4))
    call IO_setScalar("obyvel", grv_obvec(5))
    call IO_setScalar("obzvel", grv_obvec(6))
    call IO_setScalar("optxpos", grv_optvec(1))
    call IO_setScalar("optypos", grv_optvec(2))
    call IO_setScalar("optzpos", grv_optvec(3))
    call IO_setScalar("optxvel", grv_optvec(4))
    call IO_setScalar("optyvel", grv_optvec(5))
    call IO_setScalar("optzvel", grv_optvec(6))
    call IO_setScalar("oobxpos", grv_oobvec(1))
    call IO_setScalar("oobypos", grv_oobvec(2))
    call IO_setScalar("oobzpos", grv_oobvec(3))
    call IO_setScalar("oobxvel", grv_oobvec(4))
    call IO_setScalar("oobyvel", grv_oobvec(5))
    call IO_setScalar("oobzvel", grv_oobvec(6))
    call IO_setScalar("grv_oexactvec_x",  grv_oexactvec(1))
    call IO_setScalar("grv_oexactvec_y",  grv_oexactvec(2))
    call IO_setScalar("grv_oexactvec_z",  grv_oexactvec(3))
    call IO_setScalar("grv_oexactvec_vx", grv_oexactvec(4))
    call IO_setScalar("grv_oexactvec_vy", grv_oexactvec(5))
    call IO_setScalar("grv_oexactvec_vz", grv_oexactvec(6))
    call IO_setScalar("grv_ompolevec_x",  grv_ompolevec(1))
    call IO_setScalar("grv_ompolevec_y",  grv_ompolevec(2))
    call IO_setScalar("grv_ompolevec_z",  grv_ompolevec(3))
    call IO_setScalar("grv_ompolevec_vx", grv_ompolevec(4))
    call IO_setScalar("grv_ompolevec_vy", grv_ompolevec(5))
    call IO_setScalar("grv_ompolevec_vz", grv_ompolevec(6))
    call IO_setScalar("bndxpos", grv_boundvec(1))
    call IO_setScalar("bndypos", grv_boundvec(2))
    call IO_setScalar("bndzpos", grv_boundvec(3))
    call IO_setScalar("bndxvel", grv_boundvec(4))
    call IO_setScalar("bndyvel", grv_boundvec(5))
    call IO_setScalar("bndzvel", grv_boundvec(6))
    call IO_setScalar("dynrefinemax", grv_dynRefineMax)
    call IO_setScalar("ptmass", grv_ptmass)
    call IO_setScalar("optmass", grv_optmass)
    call IO_setScalar("ototmass", grv_ototmass)
    call IO_setScalar("momacc_x", grv_momacc(1))
    call IO_setScalar("momacc_y", grv_momacc(2))
    call IO_setScalar("momacc_z", grv_momacc(3))
    call IO_setScalar("angmomacc_x", grv_angmomacc(1))
    call IO_setScalar("angmomacc_y", grv_angmomacc(2))
    call IO_setScalar("angmomacc_z", grv_angmomacc(3))
    call IO_setScalar("eneracc", grv_eneracc)
    call IO_setScalar("massacc", grv_massacc)
end subroutine Gravity_sendOutputData

