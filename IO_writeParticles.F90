!!****if* source/IO/IOParticles/IO_writeParticles
!!
!! NAME
!!
!! IO_writeParticles
!!
!!
!! SYNOPSIS
!!
!! IO_writeParticles(logical(in) :: particlesToCheckpoint)
!! 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!     This routine writes out the particle data.  This is a general routine
!!     which will then call io_writeParticleData which is a routine specific
!!     to either hdf5 or parallel netcdf, or a plain old fortran write.
!!
!!     Particle data is written in two places.  First, particle data must be
!!     included in the checkpoint (restart) files in order to capture the 
!!     entire state of the simulation and to be able to restart.  Often, however,
!!     particle data needs to be written very frequently and is written to its 
!!     own particle plotfile without any other mesh data. It is possible to 
!!     separate the particles from the 
!!     checkpoint files because particles are not associated with the mesh data.
!!
!!     Particles are written in double precision to a checkpoint file and
!!     single precision in particle plotfiles  
!!
!!     The functionality of writing the particle data to a checkpoint file or
!!     a particle plotfile is the same so to eliminate code duplication we 
!!     have added an argument to the IO_writeParticles interface to indicate
!!     if we are writing particles to a checkpoint file or to a particle plotfile.
!!     (see below)
!! ARGUMENTS
!!
!!
!!      particlesToCheckpoint -   logical value - .true. if particles are
!!                               written to a checkpoint. .false. if particles
!!                               are written to a particle plot file
!!    
!!
!! NOTES
!!   To control particle output there are a few different parameters
!!   
!!
!!***

subroutine IO_writeParticles( particlesToCheckpoint)

  use IO_data, ONLY : io_chkptFileID, io_justCheckpointed, io_globalMe
  use Driver_interface, ONLY : Driver_getSimTime, Driver_getNStep
  use Logfile_interface, ONLY : Logfile_stamp
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY : Particles_getLocalNum, &
    Particles_updateAttributes, Particles_manageLost, Particles_sinkSyncWithParticles
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use IOParticles_data, ONLY : io_particleFileNumber, &
       io_particleFileIntervalTime, io_nextParticleFileTime, &
       io_particleFileIntervalStep, io_nextParticleFileStep, &
       io_dumpParticleFileExist, useParticles, &
       io_particleFileIntervalZ, io_nextParticleFileZ
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use IO_interface, ONLY : IO_updateScalars  
  use io_ptInterface, ONLY : io_ptWriteParticleData
  use Particles_data, ONLY : pt_numLocal
  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"
#include "Particles.h"
  logical, intent(in) :: particlesToCheckpoint

  integer :: i , ierr

  character (len=MAX_STRING_LENGTH) :: filename
  integer :: pptFileID

 ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff
  integer           :: particleOffset, localNumParticles
  integer           :: globalNumParticles
  character (len=OUTPUT_PROP_LENGTH) :: partAttributeLabels(NPART_PROPS)

  real          :: lsimTime   !local sim time variable
  real, save    :: oldsimTime = 0.0  !initial value does not matter
  integer       :: lsimGen    !local sim generation variable
  integer, save :: oldsimGen = -1 !different from any valid generation, to trigger update at startup.

  logical, save :: firstCall = .true.
  logical restart
  logical :: forceParticleFile

  real :: currentRedshift
  integer :: nstep



  !if we turn off particles return!
    
  logical :: gatheredUseParticles
  call MPI_ALLREDUCE(useParticles, gatheredUseParticles, 1, MPI_LOGICAL, &
                     MPI_LOR, MPI_COMM_WORLD, ierr)
  if(.not. gatheredUseParticles) then
     return
  end if

  call Driver_getNStep(nstep)

  !don't bother if you don't want output
  if ((.NOT. particlesToCheckpoint) &
       .AND. (io_particleFileIntervalStep == 0) &
       .AND. (abs(io_particleFileIntervalTime) .lt. TINY(1.0)) &
       .AND. (abs(io_particleFileIntervalZ) - HUGE(1.) .lt. TINY(1.))) then
    return
  end if

  forceParticleFile = .false.

  !collect data from all the units for outputting
  call IO_updateScalars()
  call Driver_getSimTime(lsimTime, lsimGen)


  !Are we at the top of a restart?
  restart = firstCall
  firstCall = .false. 

  if (.NOT. particlesToCheckpoint) then
  !------------------------------------------------------------------------------
  ! check to see if ".dump_particle_file" exists
  !------------------------------------------------------------------------------
  
     
  ! only the MASTER_PE should look for this file, since some systems have disk
  ! hanging off each node, that is not linked into a single namespace
  
     if ( io_globalMe == MASTER_PE ) then
        inquire(file = ".dump_particle_file", EXIST=io_dumpParticleFileExist)
     end if
  
  ! broadcast this result to all processors now
  
     call MPI_Bcast(io_dumpParticleFileExist, 1, MPI_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  
     if ( io_dumpParticleFileExist ) then
        forceParticleFile = .true.
        if ( io_globalMe == MASTER_PE ) then
           write(*,*) '[IO_writeParticles] .dump_particle_file file found: writing particle file.'
        end if
          
     endif

  end if
  ! pull in sink particles
  call Particles_sinkSyncWithParticles(sink_to_part=.true.)


  !map the particle property names from an int to a string
  do i=1, NPART_PROPS
     call Simulation_mapIntToStr(i, partAttributeLabels(i),MAPBLOCK_PART)
  end do

  call Particles_manageLost(PART_EXPAND)

  call Particles_getLocalNum(1, localNumParticles)
#ifdef DEBUG_PARTICLES
  print*,'in IO the number of particles ',localNumParticles
#endif
  !get particle offset for each processor
  call io_getParticleOffset( localNumParticles, globalNumParticles, particleOffset)  


  if(globalNumParticles <= 0 .and. io_globalMe == MASTER_PE) then
     print *, "WARNING: globalNumParticles = 0!!!"
  end if

  !check to see if we are writing particles to the checkpoint file
  !or their own particle plotfile
  if(particlesToCheckpoint) then
     if(oldsimTime /= lsimTime .OR. oldsimGen /= lsimGen) then
        call Particles_updateAttributes()
        if (io_globalMe == MASTER_PE) then
           call Logfile_stamp( 'done called Particles_updateAttributes()', "[IO_writeParticles]")
        end if
        oldsimTime=lsimTime
        oldsimGen =lsimGen
     end if

     call io_ptWriteParticleData( io_chkptFileID, globalNumParticles, &
          localNumParticles, particleOffset, partAttributeLabels, particlesToCheckpoint)

  else
     !we are writing a particle file, but first see if it is time to output
     !Note this is here because particles are not included in every simulation
     !and variables like io_particleFileIntervalTime were not being initialized and caused errors


     !get the sim time and store it in a local variable
     call Driver_getSimTime(lsimTime)
     call Cosmology_getRedshift(currentRedshift)
     
     if ((io_particleFileIntervalTime > 0.e0 .AND. lsimTime >= io_nextParticleFileTime) &
          .OR. io_justCheckpointed .OR. io_dumpParticleFileExist .OR. &
          (io_nextParticleFileStep == nstep)  .OR. &
          (io_particleFileIntervalZ < HUGE(1.) .AND. currentRedshift <= io_nextParticleFileZ)) then

        
        call Timers_start("writeParticles")

        ! open the file
        call io_initFile(io_particleFileNumber, pptFileID, filename, "_part_", forceParticleFile)
        
        
        if (io_globalMe == MASTER_PE) then
           write (strBuff(1,1), "(A)") "type"
           write (strBuff(1,2), "(A)") "particles"
           write (strBuff(2,1), "(A)") "name"
           write (strBuff(2,2), "(A)") trim(filename)
           call Logfile_stamp( strBuff, 2, 2, "[IO_writeParticles] open")
        end if
        
        
        if(oldsimTime /= lsimTime .OR. oldsimGen /= lsimGen) then
           call Particles_updateAttributes()
           if (io_globalMe == MASTER_PE) then
              call Logfile_stamp( 'done called Particles_updateAttributes()', "[IO_writeParticles]")
           end if
           oldsimTime=lsimTime
           oldsimGen =lsimGen
        end if

        call io_ptWriteParticleData( pptFileID, globalNumParticles, &
             localNumParticles, particleOffset, partAttributeLabels, particlesToCheckpoint)
        
        call io_closeFile( pptFileID)
        
        !increment the particle plotfile number
        !only increment if we are writing a single particle dataset to a file
        io_particleFileNumber = io_particleFileNumber + 1
        
        call Timers_stop("writeParticles")

        if (lsimTime >= io_nextParticleFileTime) then
           io_nextParticleFileTime = io_nextParticleFileTime + &
                (int((lsimTime-io_nextParticleFileTime)/io_particleFileIntervalTime) + 1)*io_particleFileIntervalTime
        
        end if

        if(currentRedshift <= io_nextParticleFileZ) then
           io_nextParticleFileZ = io_nextParticleFileZ - &
                (int((currentRedshift - io_nextParticleFileZ)/io_particleFileIntervalZ) + 1) * io_particleFileIntervalZ
        end if
        
        if(io_nextParticleFileStep == nstep) then
           io_nextParticleFileStep = io_nextParticleFileStep + io_particleFileIntervalStep
        end if
        
        if (io_globalMe == MASTER_PE) then
           write (strBuff(1,1), "(A)") "type"
           write (strBuff(1,2), "(A)") "particles"
           write (strBuff(2,1), "(A)") "name"
           write (strBuff(2,2), "(A)") trim(filename)
           call Logfile_stamp( strBuff, 2, 2, "[IO_writeParticles] close")
           print *, '*** Wrote particle file to ', trim(filename),  ' ****' 
           write (strBuff(1,1), "(A)") "file"
           write (strBuff(1,2), "(A)") trim(filename)
        end if
     end if
  end if
  call Particles_manageLost(PART_COLLAPSE)

  ! detach sink particles [CF off-domain fix]
  call Particles_sinkSyncWithParticles(sink_to_part=.false.)

  return

end subroutine IO_writeParticles
