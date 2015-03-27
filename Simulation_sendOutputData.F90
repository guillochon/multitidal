!!****f* source/Simulation/Simulation_sendOutputData
!!
!! NAME
!!  Simulation_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Simulation_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Simulation unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Simulation_sendOutputData()
    use IO_interface, ONLY :  IO_setScalar
    use Simulation_data, ONLY : sim_fixedPartTag, sim_mpoleVX, sim_mpoleVY, sim_mpoleVZ

    implicit none

    call IO_setScalar("fixedparttag", sim_fixedPartTag)
    call IO_setScalar("sim_mpolevx", sim_mpoleVX)
    call IO_setScalar("sim_mpolevy", sim_mpoleVY)
    call IO_setScalar("sim_mpolevz", sim_mpoleVZ)
end subroutine Simulation_sendOutputData

