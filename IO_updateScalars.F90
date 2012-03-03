!!****if* source/IO/IOMain/IO_updateScalars
!!
!! NAME
!!
!!  IO_updateScalars
!!
!!
!! SYNOPSIS
!!
!!  IO_updateScalars()
!!
!!
!!
!! DESCRIPTION
!!
!!  It calls all of the Unit output routines that allow a user to write out
!!  scalar data.  The purpose is to send each unit's scalars to a scalar list
!!  to be output together in a checkpoint file.
!!
!! 
!! ARGUMENTS
!!
!!
!! NOTES
!!
!!
!!
!! SIDE EFFECTS
!!
!!
!!***

subroutine IO_updateScalars()
  use Hydro_interface, ONLY : Hydro_sendOutputData
  use Gravity_interface, ONLY : Gravity_sendOutputData
  use Simulation_interface, ONLY : Simulation_sendOutputData
  use Cosmology_interface, ONLY : Cosmology_sendOutputData
  use Driver_interface, ONLY : Driver_sendOutputData
  use Particles_interface, ONLY : Particles_sendOutputData
  use Grid_interface, ONLY : Grid_sendOutputData
  use IO_interface, ONLY : IO_sendOutputData
implicit none

  call Grid_sendOutputData()

  call Driver_sendOutputData()

  call IO_sendOutputData()

  call Hydro_sendOutputData()

  call Particles_sendOutputData()

  !!call IO_outputScalars()
  
  call Cosmology_sendOutputData()

  call Gravity_sendOutputData()

  !! add in other calls to output
  call Simulation_sendOutputData()

  
end subroutine IO_updateScalars
