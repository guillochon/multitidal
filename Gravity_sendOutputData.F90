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
    use gr_mpoleData, ONLY : Mtot, Xcm, Ycm, Zcm, oXcm, oYcm, oZcm
    use Gravity_data, ONLY: grv_ptvec, grv_obvec, grv_optvec, grv_oobvec, grv_dynRefineMax
    implicit none

    call IO_setScalar("Mtot", Mtot)
    call IO_setScalar("Xcm", Xcm)
    call IO_setScalar("Ycm", Ycm)
    call IO_setScalar("Zcm", Zcm)
    call IO_setScalar("oXcm", oXcm)
    call IO_setScalar("oYcm", oYcm)
    call IO_setScalar("oZcm", oZcm)
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
    call IO_setScalar("dynrefinemax", grv_dynRefineMax)
end subroutine Gravity_sendOutputData

