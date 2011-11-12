!!****if* source/IO/IOMain/IO_output
!!
!! NAME
!!
!!  IO_output
!!
!!
!! SYNOPSIS
!!
!!
!!  IO_output(integer(in)           :: myPE, 
!!            integer(in)           :: numProcs, 
!!            real(in)              :: simTime, 
!!            real(in)              :: dt,
!!            integer(in)           :: nstep,
!!            integer(in)           :: nbegin,
!!            logical(out)          :: endRun,
!!            integer(in), optional :: outputType)
!!
!!
!! DESCRIPTION
!!
!!  This routine is called after every sweep to check if it is time to
!!  output.  Checkpoints and plotfiles and .dat files
!!  are handled through this function.
!!  
!!  
!!  A checkpoint can be triggered in a few different ways.  First, enough
!!  wall clock time has elapsed since the code has started running that
!!  we want to force a checkpoint now.  This is controlled by the
!!  wall_clock_checkpoint runtime parameter.  This is always timed from
!!  the start of execution of the code, and not since the last checkpoint,
!!  so we can ensure that we get a checkpoint dumped right before the
!!  queue window closes.  Checkpoints are also produced after a given
!!  amount of simulation time has elapsed -- this is controlled by the
!!  checkpointFileIntervalTime runtime parameter.  Finally, a checkpoint can be forced in
!!  to be produced every checkpointFileIntervalStep timesteps or by creating .dump_checkpoint
!!  file in the execution directory.
!!
!!  Plotfiles are produced equally spaced in simulation time, tplot time
!!  units apart. They can also be produced on demand by temporarily creating
!!  .dump_plotfile execution control file.
!!  
!!
!!  After every sweep, IO_writeIntegralQuantities is called to compute some global
!!  quantities and write them to the flash.dat file.
!!
!!  A checkpoint is given a name according to the basename (specified
!!  via the basenm runtime parameter) and the filetype (pnetcdf,
!!  HDF5).  Each separate checkpoint is given a unique number suffix.  The
!!  value to start with (or restart from if this simulation was a restart
!!  from a previous checkpoint) is set by the checkpointFileNumber runtime parameter.
!!  Or if it is a restart, checkpointFileNumber can be saved as a scalar in the checkpoint
!!  file.
!!
!!  Execution will be aborted immediately if .kill file is present
!!  in the execution directory. If .dump_restart is present then a
!!  checkpoint file is saved before aborting execution.
!!
!!  For more detail about runtime parameters controlling IO output
!!  look at the IO/common Config file or the setup_params text file
!!  written to the object directory.
!!
!! ARGUMENTS
!!
!!  myPE - current processor
!!  numProcs - number of processors running the simulation
!!  simTime - simulation time
!!  dt - timestep
!!  nstep - current time step number
!!  nbegin - beginning time step number, if starting from scratch this is 0,
!!      it could be different in case of restart
!!  endRun - will be set to .TRUE. on return if existence of a .dump_restart or
!!      a .kill file was detected; .FALSE. otherwise.
!!  outputType - an integer that denotes the type of output files you expect to 
!!               generate with this call. If this argument is omitted IO_output
!!               will test to see if it is time to output every kind of file.
!!
!! NOTES
!!
!!  Variables that start with "io_", like io_checkpointFileNumber and io_checkpointFileIntervalTime
!!  are located in the IO_data fortran 
!!  module.  The "io_" is meant to indicated that this variable has
!!  IO Unit scope.  Variable without the "io_" in front are local
!!  variables to this subroutine.
!!
!!  The outputType argument has a series of constants declared in constants.h which
!!  allow you to state which filed you wish to have output (povided that the normal 
!!  output conditions are met) by this call.  The options are:
!!
!!  CHECKPOINT_FILE_ONLY 
!!  PLOTFILE_ONLY 
!!  PRTICLE_FILE_ONLY 
!!  CHECKPOINT_AND_PLOTFILE 
!!  CHECKPOINT_AND_PARTICLEFILE 
!!  PLOTFILE_AND_PARTICLEFILE 
!!  ALL_FILES 
!!
!!  If the optional outputType argument is omitted, ALL_FILES is assumed.
!!
!! SIDE EFFECTS
!!
!!  For visualization purposes, the data is restricted up the entire tree,
!!  so the data is valid on all levels.  This allows multi-resolution vis
!!  techniques to be applied.  Immediately after checkpointing or writing
!!  a plotfile is the only time the data is guaranteed to be valid on all
!!  levels.
!!
!! SEE ALSO
!!   IO_writeCheckpoint, IO_writeIntegralQuantities, IO_writeParticles,
!!   IO_writePlotfile, IO_writeUserArray, IO_setScalar
!!
!!***


subroutine IO_output(myPE, numProcs, simTime, dt, nstep, nbegin, endRun, outputType)

  use IO_data, ONLY : io_nextCheckPointStep, io_nextCheckpointTime, &
       io_nextPlotFileStep, io_nextPlotfileTime, &
       io_checkpointFileNumber, io_plotFileNumber, &
       io_checkpointFileIntervalTime, io_checkpointFileIntervalStep, &
       io_plotfileIntervalTime, io_plotfileIntervalStep, io_integralFreq, &
       io_redshift, io_lastWallClockCheckpoint, io_wallClockCheckpoint, io_rollingCheckpoint, &
       io_justCheckpointed, io_memoryStatFreq, &
       io_nextCheckpointZ, io_nextPlotFileZ, io_checkpointFileIntervalZ, io_plotfileIntervalZ, &
       io_alwaysComputeUserVars, io_alwaysRestrictCheckpoint
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use Driver_interface, ONLY : Driver_abortFlash,Driver_finalizeFlash
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Grid_interface, ONLY : Grid_restrictAllLevels, &
    Grid_computeUserVars
  use IO_interface, ONLY : IO_writeIntegralQuantities, &
    IO_writeCheckpoint, IO_writePlotfile, IO_writeParticles, IO_writeOrbitInfo
  use io_ptInterface, ONLY : io_ptCorrectNextPartTime, io_ptResetNextFile

  implicit none


#include "constants.h"
#include "Flash_mpi.h"

  integer, intent(in) :: myPE, numProcs, nbegin, nstep
  real, intent(in) :: simTime, dt
  integer, intent(in), optional :: outputType
  logical :: outputCheckpoint
  logical :: outputPlotfile
  logical :: outputParticleFile

  integer :: ii, ierr
  
  character(len=MAX_STRING_LENGTH) :: outstr
  
  real :: dtCheckpoint, wtime, sav_nextPlotfileTime, sav_nextParticleFileTime
  real :: newCPUSeconds
  real,save :: lastSimTime = HUGE(simTime)*0.5
  logical :: forceCheckpoint, filteredparticles
  logical   :: filter
  
  logical :: dumpPlotfileExist, dumpRestartExist
  logical :: dumpCheckpointExist, killExist
  logical , intent(OUT) :: endRun 
  logical :: justRestricted

  real :: currentRedshift !current cosmological redshift value
  real :: saveNextRedshift


  

  if(present(outputType)) then
     
     outputCheckpoint = IAND(outputType, CHECKPOINT_FILE_ONLY) /= 0
     outputPlotfile = IAND(outputType, PLOTFILE_ONLY) /= 0
     outputParticleFile = IAND(outputType, PARTICLE_FILE_ONLY) /= 0

  else

     outputCheckpoint = .true.
     outputPlotfile = .true.
     outputParticleFile = .true.
     
  end if

  !==============================================================================
  
  !only set this to true if we want to terminate the run.
  endRun = .false.
 
  !set to true once we have called Grid_restrictAllLevels
  justRestricted = .false.
  

  !have to get the current redshift
  call Cosmology_getRedshift(currentRedshift)

  !------------------------------------------------------------------------------
  ! Dump out memory usage statistics if we are monitoring them
  !------------------------------------------------------------------------------
  if ( (io_memoryStatFreq > 0) .AND. (mod(nstep,io_memoryStatFreq) == 0) ) call io_memoryReport(mype, numprocs)
  
  
  !------------------------------------------------------------------------------
  ! Write values of integral quantities -- these are stored in flash.dat
  !------------------------------------------------------------------------------
  
  call Timers_start("diagnostics")
  if(io_integralFreq > 0) then
     if (mod (nstep, io_integralFreq) == 0) then
        !This should happen only once per timestep
        if (abs(simTime - lastSimTime) .gt. TINY(simTime)) then
            call IO_writeIntegralQuantities(myPE, 0, simTime)
            call IO_writeOrbitInfo(myPE, 0, simTime)
        endif
        lastSimTime = simTime
     end if
  end if
  call Timers_stop("diagnostics")
  
  
  io_justCheckpointed = .false.
  !------------------------------------------------------------------------------
  ! dump out a checkpoint file if it is time (either enough timesteps
  ! have elapsed, or the simulation time since the last checkpoint is
  ! > checkpointFileIntervalTime, or the simulation redshift has changed by checkpointFileIntervalRedshift,  but not
  ! if it's the first step of a run
  !------------------------------------------------------------------------------
  
  
  if (((nstep == io_nextCheckpointStep) .or. &
       (io_checkpointFileIntervalTime > 0.e0  .and. simTime >= io_nextCheckpointTime) .or.& 
       (io_checkpointFileIntervalZ < HUGE(1.) .and. currentRedshift <= io_nextCheckpointZ) &
       .and. nstep /= nbegin) .and. outputCheckpoint) then 
     
     
     ! make all the levels have valid data -- normally we only evolve on
     ! the leaf blocks, but for visualization purposes, it is nice to be
     ! able to look at other nodes in the tree
     call Timers_start("checkpointing")
     if(.not. justRestricted) then
        if(.not. io_alwaysRestrictCheckpoint) call Grid_restrictAllLevels()
        justRestricted = .true.
     end if
    
     if(.not. io_alwaysComputeUserVars) call Grid_computeUserVars()
     
     
     ! Compute the nextPlotFileTime to be stored in the checkpoint file;
     ! this computation is repeated below.
     sav_nextPlotfileTime = io_nextPlotfileTime !save current value - KW
     if ( io_plotFileIntervalTime > 0.e0 ) then
        if (simTime >= io_nextPlotFileTime) then
           io_nextPlotFileTime = io_nextPlotFileTime + &
                (int((simTime-io_nextPlotFileTime)/io_plotFileIntervalTime) + 1)*io_plotFileIntervalTime
        end if
     end if
     
     
     if ( io_checkpointFileIntervalTime > 0.e0 ) then 
        
        do while (simTime >= io_nextCheckpointTime)
           io_nextCheckpointTime = io_nextCheckpointTime + io_checkpointFileIntervalTime
        end do
        
     end if


     !Adjust cosmological redshift
     if(io_checkpointFileIntervalZ < HUGE(1.)) then
        do while(currentRedshift <= io_nextCheckpointZ)
           io_nextCheckpointZ = io_nextCheckpointZ - io_checkpointFileIntervalZ
        end do
     end if
     
     call io_ptCorrectNextPartTime(simTime, sav_nextParticleFileTime)
     
     ! compute the checkpoint number, in case we are rolling
     io_checkpointFileNumber = mod(io_checkpointFileNumber,io_rollingCheckpoint)
     
     call IO_writeCheckpoint(myPE, numProcs)
     
     io_nextCheckpointStep = nstep + io_checkpointFileIntervalStep
     
     io_justCheckpointed = .true.
     
     io_nextPlotfileTime = sav_nextPlotfileTime !restore original value - KW
     
     call io_ptResetNextFile(sav_nextParticleFileTime)
     
     call Timers_stop("checkpointing")
     
  end if

  
  
  !------------------------------------------------------------------------------
  ! Compute the time since the last wall clock triggered checkpoint.  If
  ! too much time has elapsed (as specified by the runtime parameter 
  ! wallClockCheckpoint, then checkpoint, and store the current wall
  ! clock time.  This allows us to checkpoint just before a queue window
  ! closes, maximizing queue time. Also, prevent from writing the same
  ! checkpoint file twice.
  !------------------------------------------------------------------------------

  if ( (io_wallClockCheckpoint > 0.e0) .and. outputCheckpoint) then

     if (myPE == MASTER_PE) then

        dtCheckpoint = MPI_Wtime() - io_lastWallClockCheckpoint
        
        if (dtCheckpoint > io_wallClockCheckpoint) then
     
           if (myPE == MASTER_PE .and. .not.io_justCheckpointed)  & 
                write(*,*) 'checkpointing due to elapsed wall clock time from last checkpoint '

           forceCheckpoint = .true.
        else
           forceCheckpoint = .false.
        endif
     endif




     ! broadcast whether we are forcing a checkpoint to all the other processors
     call MPI_Bcast(forceCheckpoint, 1, MPI_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  
     if (forceCheckpoint) then
        
        if ( .not.io_justCheckpointed) then
        
           call Timers_start("checkpointing")
           if(.not. justRestricted) then
              if(.not. io_alwaysRestrictCheckpoint) call Grid_restrictAllLevels()
              justRestricted = .true.
           end if
           
           if(io_alwaysComputeUserVars) call Grid_computeUserVars()


           ! Compute the nextPlotFileTime to be stored in the checkpoint file;
           ! this computation is repeated below.
           sav_nextPlotfileTime = io_nextPlotfileTime !save current value - KW
           if ( io_plotFileIntervalTime > 0.e0 ) then
              if (simTime >= io_nextPlotFileTime) then
                 io_nextPlotFileTime = io_nextPlotFileTime + &
                      (int((simTime-io_nextPlotFileTime)/io_plotFileIntervalTime) + 1)*io_plotFileIntervalTime
              end if
           end if
           
           call io_ptCorrectNextPartTime(simTime, sav_nextParticleFileTime)
           io_checkpointFileNumber = mod(io_checkpointFileNumber,io_rollingCheckpoint)
           call IO_writeCheckpoint(myPE, numProcs)
           io_justCheckpointed = .true.
           
           io_nextPlotfileTime = sav_nextPlotfileTime !restore original value - KW
           call io_ptResetNextFile(sav_nextParticleFileTime)
           call Timers_stop("checkpointing")
        
         end if
 
           ! store the wallclock time of this checkpoint in the runtime scratch database
           io_lastWallClockCheckpoint = MPI_Wtime()
        
     endif
     
  end if



  !------------------------------------------------------------------------------
  ! check if it is time to dump out a plot file, if so, go for it.
  !------------------------------------------------------------------------------
  

     
  if (((nstep == io_nextPlotFileStep) .OR. &
       (io_plotFileIntervalTime > 0.e0 .AND. simTime >= io_nextPlotFileTime) .OR. &
       (io_plotFileIntervalZ < HUGE(1.) .AND. currentRedshift <= io_nextPlotFileZ)) .and. outputPlotfile) then
     
     call Timers_start("plotfile")
     if(.not. justRestricted) then
      !  call Grid_restrictAllLevels()
        justRestricted = .true.
     end if

     !call Grid_computeUserVars()


     if ( io_plotFileIntervalTime > 0.e0 ) then

        if (simTime >= io_nextPlotFileTime) then

           io_nextPlotFileTime = io_nextPlotFileTime + &
                (int((simTime-io_nextPlotFileTime)/io_plotFileIntervalTime) + 1)*io_plotFileIntervalTime

        end if
     end if

     if(io_plotFileIntervalZ < HUGE(1.)) then
        if(currentRedshift <= io_nextPlotFIleZ) then
           io_nextPlotFIleZ = io_nextPlotFileZ + &
                (int((currentRedshift+io_nextPlotFileZ)/io_plotFIleIntervalZ) + 1)*io_plotFileIntervalZ

        end if
     end if

     call IO_writePlotfile(myPE, numProcs)

     call Timers_stop("plotfile")
     io_nextPlotFileStep = nstep + io_plotFileIntervalStep


  end if





  !call IO_writeParticles.  If particles are not included in a simulation then
  !a stub will be called here
  !checks to see if it is time to write particles are done internally in the
  !IO_writeParticles routine
 
    
  !Logic is being moved out of IO_writeParticles.  When it is called, it had better be used!
  if(outputParticleFile) call IO_writeParticles(myPE, numProcs, .false.)

  

  !------------------------------------------------------------------------------
  ! if the file ".dump_plotfile" exists, write plot file #9999 and continue
  !------------------------------------------------------------------------------
  
  ! only the MASTER_PE should look for this file, since some systems have disk 
  ! hanging off each node, that is not linked into a single namespace
  
  if (myPE == MASTER_PE) then
     inquire(file = ".dump_plotfile", EXIST=dumpPlotfileExist)
  end if
  
  ! broadcast this result to all processors now
  
  call MPI_Bcast(dumpPlotfileExist, 1, MPI_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  
  if ( outputPlotfile .and. dumpPlotfileExist ) then
     
     if ( myPE == MASTER_PE ) then
        write(*,*) '[OUTPUT] .dump_plotfile file found: writing plotfile.'
     end if
     
     call Timers_start("plotfile")
     if(.not. justRestricted) then
        call Grid_restrictAllLevels()
        justRestricted = .true.
     end if
     
     call Grid_computeUserVars()
     !Note that this plotfile is considered "forced"
     call IO_writePlotfile(myPE, numProcs, .true.)
     call Timers_stop("plotfile")
     
  end if



  !------------------------------------------------------------------------------
  ! if the file ".dump_checkpoint" exists, dump checkpoint file 9999 and continue
  !------------------------------------------------------------------------------
  
  ! only the MASTER_PE should look for this file, since some systems have disk 
  ! hanging off each node, that is not linked into a single namespace

  if ( myPE == MASTER_PE) then
     inquire(file = ".dump_checkpoint", EXIST=dumpCheckpointExist)
  end if

  ! broadcast this result to all processors now

  call MPI_Bcast(dumpCheckpointExist, 1, MPI_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  if ( outputCheckpoint .and. dumpCheckpointExist ) then

     if ( myPE == MASTER_PE ) then
        write(*,*) '[IO_OUTPUT] .dump_checkpoint file found: writing checkpoint file.'
     end if
 
     call Timers_start("checkpointing")
     if(.not. justRestricted) then
        if (.not. io_alwaysRestrictCheckpoint) call Grid_restrictAllLevels()
        justRestricted = .true.
     end if

     if(.not. io_alwaysComputeUserVars) call Grid_computeUserVars()

     call IO_writeCheckpoint(myPE, numProcs)
     call Timers_stop("checkpointing")

  end if

  !------------------------------------------------------------------------------
  ! if the file ".dump_restart" exists, it means that the machine is going down,
  ! and we want to dump an output file right now and abort
  !------------------------------------------------------------------------------
  
  ! only the MASTER_PE should look for this file, since some systems have disk 
  ! hanging off each node, that is not linked into a single namespace

  if (  myPE == MASTER_PE) then
     inquire(file = ".dump_restart", EXIST=dumpRestartExist)
  end if

  ! broadcast this result to all processors now

  call MPI_Bcast(dumpRestartExist, 1, MPI_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  if ( outputCheckpoint .and. dumpRestartExist ) then

     if ( myPE == MASTER_PE ) then
        write(*,*) '[OUTPUT]  .dump_restart file found: writing checkpoint and aborting.'
     end if
    
     ! tell main loop to exit
     endRun = .true.
     call Timers_start("checkpointing")
     if(.not. justRestricted) then
        if(.not. io_alwaysRestrictCheckpoint) call Grid_restrictAllLevels()
        justRestricted = .true.
     end if

     if(.not. io_alwaysComputeUserVars) call Grid_computeUserVars()

     io_checkpointFileNumber = mod(io_checkpointFileNumber,io_rollingCheckpoint)          
     call IO_writeCheckpoint(myPE, numProcs)
     call Timers_stop("checkpointing")

     call Logfile_stamp(myPE, '.dump_restart file found, exiting' , '[IO_OUTPUT]')

     
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     write(*,*) '[IO_output] .dump_restart file found: exiting.'

     !DO NOT CALL Driver_finalizeFlash here, will cause problems! -PR


  end if

  !------------------------------------------------------------------------------
  ! if the file ".kill" file exists, abort execution immediately
  !------------------------------------------------------------------------------
  
  ! only the MASTER_PE should look for this file, since some systems have disk 
  ! hanging off each node, that is not linked into a single namespace

  if ( myPE == MASTER_PE ) then
     inquire(file = ".kill", EXIST=killExist)
  end if
  
  ! broadcast this result to all processors now
  
  call MPI_Bcast(killExist, 1, MPI_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  
  if ( killExist ) then
     
    ! tell main loop to exit
    endRun = .true.


     call Logfile_stamp(myPE, '.kill file found, exiting' , '[IO_OUTPUT]')



     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     write(*,*) '[IO_output] .kill file found: exiting.'

     !DO NOT CALL Driver_finalizeFlash here, will cause problems! -PR


  end if


  return
end subroutine IO_output
