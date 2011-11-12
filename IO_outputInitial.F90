!!****if* source/IO/IOMain/IO_outputInitial
!!
!! NAME
!!
!!  IO_outputInitial
!!
!!
!! SYNOPSIS
!!
!!
!!  IO_outputInitial(integer(in) :: myPE, 
!!                   integer(in) :: numProcs,
!!                   integer(in) :: nbegin,
!!                   real(in) :: initialSimTime  
!!                  
!!
!!
!! DESCRIPTION
!!
!!
!!  This routine is called before the main timestep loop.  It outputs the 
!!  initial data to a checkpoint file and plotfile, and particle plotfiles
!!
!!  If particles are not included a stub (empty) routine will be called.
!!
!!
!! ARGUMENTS
!!
!!  myPE - current processor
!!  numProcs - number of processors running the simulation
!!  nbegin - initial step of simulation
!!  initialSimTime - initial simulation time
!!
!!
!!***


subroutine IO_outputInitial(myPE, numProcs, nbegin, initialSimTime)

  use IO_data, ONLY : io_integralFreq, io_memoryStatFreq, &
       io_redshift, io_justCheckpointed, io_restart, &
       io_alwaysRestrictCheckpoint, io_alwaysComputeUserVars
  use Grid_interface, ONLY : Grid_restrictAllLevels, &
    Grid_computeUserVars
  use IO_interface, ONLY : IO_writeIntegralQuantities, &
    IO_writeCheckpoint, IO_writePlotfile, IO_writeParticles, IO_writeOrbitInfo

  implicit none

#include "Flash_mpi.h"

  integer, intent(in) :: myPE, numProcs, nbegin
  real, intent(in) :: initialSimTime
  logical :: forcePlotfile
  
  forcePlotfile = .false.

  !------------------------------------------------------------------------------
  ! Dump out memory usage statistics if we are monitoring them,
  ! BEFORE opening files for output for the first time.
  !------------------------------------------------------------------------------
  if (io_memoryStatFreq > 0) call io_memoryReport(myPE, numProcs)

  !write the diagnostic quantities for the .dat file
  if(io_integralFreq > 0) then
     call IO_writeIntegralQuantities(myPE, 1, initialSimTime)
     call IO_writeOrbitInfo(myPE, 1, initialSimTime)
  end if


  !Ensure valid data throughout grid and ancestor blocks
  if(.not. io_alwaysRestrictCheckpoint) call Grid_restrictAllLevels()
  if(.not. io_alwaysComputeUserVars) call Grid_computeUserVars()
  
  if(.not. io_restart) then
     call IO_writeCheckpoint(myPE, numProcs)
     io_justCheckpointed = .true.
  else
     io_justCheckpointed = .false.
  end if

  if( io_restart) forcePlotfile = .true.
  call IO_writePlotfile(myPE, numProcs, forcePlotfile)

  call IO_writeParticles(myPE, numProcs, .false.)

  !------------------------------------------------------------------------------
  ! Dump out memory usage statistics again if we are monitoring them,
  ! AFTER having written the initial output files.
  !------------------------------------------------------------------------------
  if (io_memoryStatFreq > 0) call io_memoryReport(myPE, numProcs)

end subroutine IO_outputInitial
