!!****if* source/IO/IOMain/IO_writePlotfile
!!
!!
!!
!!
!! NAME
!!
!!  IO_writePlotfile
!!
!!
!! SYNOPSIS
!!
!!  IO_writePlotfile(logical(in), optional :: forced)
!!
!!
!!
!! DESCRIPTION
!!
!!  This is a generic call to write the important simulation data to a
!!  plotfile file.  A plotfile file writes a few different types of
!!  data to a file, first the physical data like, temperature, pressure, density
!!  etc.  for all cells on the grid.
!!  Secondly, in order to recreate the simulation from a plotfile file a
!!  number of other single quantities are needed as well.  We call these
!!  scalar values which include simTime, dt, nstep, globalNumBlocks etc.
!!  We also store descriptive strings that describe the simulation run.
!!
!!  The same IO_writePlotfile routine is called regardless of the type of
!!  file being written, (such as hdf5 parallel, hdf5 serial or pnetcdf)
!!  IO_writePlotfile prepares the Grid_ioData (like getting the
!!  globalNumBlocks) and collects the scalars wanting to be stored
!!  from each unit. 
!!  IO_writePlotfile then calls three methods, io_initFile, io_writeData, 
!!  and io_closeFile.  Each of these routines _is_ specific
!!  to the type of io library used and have their own implementation.  
!!  In addition, io_writeData has its own
!!  implementation for io library and type of grid (UG, Paramesh, or other)
!!
!!  Since plotfiles are used for visualization purposes, and not for 
!!  restarting a run, to keep plotfile sizes manageable, data is 
!!  written out in single precision
!!
!!  In FLASH IO_writePlotfile is called from IO_output (or IO_outputInitial or
!!  IO_outputFinal) IO_output checks whether it is time to output a plotfile.
!!  The runtime parameters that control writing plotfiles are tplot, the 
!!  simulation time between plotfiles and nplot, the number of timesteps 
!!  between plotfiles.
!!  
!!
!!  We have put IO_writePlotfile in the API because a user may want to write
!!  a plotfile at another time or for another reason without having to go through
!!  IO_output.  For most flash users IO_writePlotfile will only ever be
!!  called through IO_output and possibly IO_outputFinal and IO_outputInitial.
!!
!! ARGUMENTS
!! 
!!  forced - should this be considered a "forced" plotfile
!!
!! NOTES
!!  
!!  We have added functionality to compress plotfile data further.  Byte packing
!!  is currently only available with hdf5 and the uniform grid.  To use
!!  byte packing set the runtime parameter bytePack to .true. in your flash.par
!!
!!  For those familiar with FLASH2, breaking up the plotfile routine into
!!  these four different methods is a change.  Because FLASH3 now supports
!!  different grid packages and we are committed to supporting both
!!  hdf5 and parallel netCDF having each grid and io library writing its
!!  own plotfile file proved to be a lot of code duplication.  We believe
!!  that while dividing up the plotfile routines created more files it 
!!  will in the end be easier to maintain.
!!
!!***


subroutine IO_writePlotfile( forced)

  use IO_data, ONLY : io_plotFileNumber, io_unklabels, &
       io_doublePrecision, io_nPlotVars, io_forcedPlotFileNumber, &
       io_ignoreForcedPlot, io_flashRelease, io_globalMe, io_wrotePlot, &
       io_oldPlotFileName
  use Logfile_interface, ONLY : Logfile_stampMessage, Logfile_stamp
  use Grid_interface, ONLY : Grid_computeUserVars, &
    Grid_restrictAllLevels
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use IO_interface, ONLY : IO_updateScalars, IO_writeParticles, IO_writeUserArray
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
  include "Flash_mpi.h"

  logical, optional, intent(in) :: forced
  
  integer :: fileID
  logical :: fileIsForced
  character (len=MAX_STRING_LENGTH) :: filename
  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff
  logical :: restrictNeeded


  if (present(forced)) then
     fileIsForced = (forced .and.(.not. io_ignoreForcedPlot))
  else
     fileIsForced = .false.
  end if

  ! If FLASH version string has not been initialized yet...
  if (io_flashRelease(1:1) == ' ') call io_prepareSimInfo()
  
  ! Ensure the Grid is in a consistient state for plotfile output.
  call io_restrictBeforeWrite( restrictNeeded)
  if (restrictNeeded .eqv. .true.) then
     call Grid_restrictAllLevels()
  end if
  call Grid_computeUserVars()
  
  
  io_doublePrecision = .false.

  if((io_nPlotVars == 0) .and. (io_globalMe == MASTER_PE)) then
     print *, "WARNING: you have called IO_writePlotfile but no plot_vars are defined."
     print *, 'Put the vars you want in the plotfile in your flash.par (plot_var_1 = "dens")'
     call Logfile_stampMessage("WARNING you have called IO_writePlotfile but no plot_vars are defined.")
     call Logfile_stampMessage('put the vars you want in the plotfile in your flash.par (plot_var_1 = "dens")')
  end if

  
  call IO_updateScalars()

  call Timers_start("writePlotfile")
  
  !-----------------------------------------------------------------------------
  ! open the file
  !-----------------------------------------------------------------------------
  IO_TIMERS_START("create file")
  if(fileIsForced) then
     call io_initFile( io_forcedPlotFileNumber, fileID, filename, '_plt_cnt_', fileIsForced)
  else
     call io_initFile( io_plotFileNumber, fileID, filename, '_plt_cnt_', fileIsForced)
  end if
  IO_TIMERS_STOP("create file")

  if (io_globalMe == MASTER_PE) then
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "plotfile"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "[IO_writePlotfile] open")
  end if



  call Grid_computeUserVars()

  call io_writeData( fileID)

  !JFG
  call IO_writeParticles( .true.)

  call IO_writeUserArray()
  !End JFG

  !----------------------------------------------------------------------
  ! close the file
  !----------------------------------------------------------------------
  IO_TIMERS_START("close file")
  call io_closeFile( fileID)
  IO_TIMERS_STOP("close file")


  !increment the checkpoint number unless it is a dump checkpoint file
  !DUMP_IOFILE_NUM typically 9999 or some large number, located in constants.h 
  if((.not.fileIsForced) .and. io_plotFileNumber /= DUMP_IOFILE_NUM) then
     io_plotFileNumber = io_plotFileNumber + 1
  else
     io_forcedPlotFileNumber = io_forcedPlotFileNumber + 1
  end if

  call Timers_stop("writePlotfile")

  if (io_globalMe == MASTER_PE) then
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "plotfile"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     print *, '*** Wrote plotfile to ', trim(filename),  ' ****' 
     call Logfile_stamp( strBuff, 2, 2, "[IO_writePlotfile] close")
  end if

  io_wrotePlot = .true.
  io_oldPlotFileName = filename

  return
end subroutine IO_writePlotfile




