!!****if* source/Driver/DriverMain/Split/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This routine implements the Strang splitting scheme for time
!! advancement. A single step in the this driver 
!! includes two sweeps, the first one in order XYZ, and
!! the second one in order ZYX. This driver works with directionally
!! split operators only. The routine also controls the regridding of
!! the mesh if necessary and the simulation output.
!!
!!  
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_myPE or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif


subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_myPE, dr_numProcs, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_simGeneration, dr_fSweepDir, dr_rSweepDir,&
       dr_nstep, dr_dtOld, dr_dtNew, dr_restart, dr_elapsedWCTime, &
       dr_redshiftInitial, dr_redshiftFinal, dr_redshift, dr_redshiftOld, &
       dr_useRedshift
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close, Logfile_stampMessage
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
      Grid_getListOfBlocks, Grid_updateRefinement
  use Hydro_interface, ONLY : Hydro
  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks
  use IO_interface, ONLY :IO_output,IO_outputFinal
  use Cosmology_interface, ONLY : Cosmology_redshiftHydro, &
      Cosmology_solveFriedmannEqn, Cosmology_getRedshift
  use gr_mpoleData, ONLY: Xcm, Ycm, Zcm
  use Gravity_data, ONLY: grv_mode, orb_t, orb_dt, grv_ptvec, grv_obvec
  use Simulation_data, ONLY: sim_tRelax
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(4,2) :: strBuff
  character(len=15) :: numToStr
  character(len=200) :: logstr
  
  logical :: gridChanged
  logical :: endRun !Should we end our run on this iteration?

  double precision :: tinitial
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: lb


  endRun = .false.
  call Logfile_stamp(dr_myPE, 'Entering evolution loop' , '[Driver_evolveFlash]')


  call Timers_start("evolution")

!!******************************************************************************
!! Start of Evolution Loop
!!******************************************************************************

  do dr_nstep = dr_nbegin, dr_nend


     !!Step forward in time. See bottom of loop for time step calculation.
     
     call Grid_getLocalNumBlks(localNumBlocks)
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)
     if (dr_myPE == MASTER_PE) then

        write (numToStr(1:), '(I10)') dr_nstep
        write (strBuff(1,1), "(A)") "n"
        write (strBuff(1,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_simTime
        write (strBuff(2,1), "(A)") "t"
        write (strBuff(2,2), "(A)") trim(adjustl(numToStr))
        
        if (.not. dr_useRedshift) then

           write (numToStr(1:), "(1PE12.6)") dr_dt
           write (strBuff(3,1), "(A)") "dt"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           call Logfile_stamp(dr_myPE, strBuff(1:3,:), 3, 2, "step")

        else

           write (numToStr(1:), "(F8.3)") dr_redshift
           write (strBuff(3,1), "(A)") "z"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           write (numToStr(1:), "(1PE12.6)") dr_dt
           write (strBuff(4,1), "(A)") "dt"
           write (strBuff(4,2), "(A)") trim(adjustl(NumToStr))
           
           call Logfile_stamp(dr_myPE, strBuff, 4, 2, "step")

        endif

     end if

     !!--------------------------------------------------------------------
     !!- Start Physics Sequence
     !!--------------------------------------------------------------------
#ifdef DEBUG_DRIVER
     print*, 'going into Hydro/MHD'
#endif

     call Timers_start("cosmology")
     call Cosmology_solveFriedmannEqn(dr_simTime, dr_dt)
     call Timers_stop("cosmology")

     call RuntimeParameters_get('tinitial',tinitial)

     if (dr_simTime .gt. tinitial + sim_tRelax) then
         grv_mode = 1
         if (dr_simTime .lt. tinitial + sim_tRelax) then
             orb_t = tinitial + sim_tRelax
             orb_dt = dr_simTime + dr_dt - tinitial - sim_tRelax
         else
             orb_t = dr_simTime
             orb_dt = dr_dt
         endif
         call Orbit_update
     endif

     dr_simTime = dr_simTime + dr_dt

     dr_simGeneration = 0

     call Timers_start("hydro")
#ifdef DEBUG_DRIVER
     print*,'going into hydro'
#endif
     call Hydro(dr_myPE, dr_numProcs, blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, dr_fSweepDir)
     call Timers_stop("hydro")

     
#ifdef DEBUG_DRIVER
     print*, 'return from Hydro/MHD timestep'
#endif

     call Timers_start("sourceTerms")
     call Driver_sourceTerms(blockCount, blockList, dr_dt)
     call Timers_stop("sourceTerms")
#ifdef DEBUG_DRIVER
     print*,'done source terms'
     print*, 'return from Drivers_sourceTerms '
#endif
     call Timers_start("Particles_advance")
     call Particles_advance(dr_dtOld, dr_dt)
#ifdef DEBUG_DRIVER
     print*, 'return from Particles_advance '
#endif
     call Timers_stop("Particles_advance")
   
     call Gravity_potentialListOfBlocks(blockCount,blockList)
#ifdef DEBUG_DRIVER
     print*, 'return from Gravity_potential '
#endif

     call Bound_mass(blockCount, blockList)
     if (dr_simTime - dr_dt .gt. tinitial + sim_tRelax) then
         grv_mode = 2
         if (dr_simTime - dr_dt .lt. tinitial + sim_tRelax) then
             orb_t = tinitial + sim_tRelax
             orb_dt = dr_simTime - tinitial - sim_tRelax
         else
             orb_t = dr_simTime - dr_dt
             orb_dt = dr_dt
         endif
         call Orbit_update
     endif

     call Timers_start("cosmology")
     call Cosmology_redshiftHydro(dr_myPE, dr_numprocs, blockCount, blockList)
     call Timers_stop("cosmology")

     call Timers_start("compute dt")
     call Driver_computeDt(dr_myPE, dr_numProcs,  &
                         dr_nbegin, dr_nstep, &
                         dr_simTime, dr_dtOld, dr_dtNew)
     call Timers_stop("compute dt")
     dr_dt = dr_dtNew                                        ! store new

     dr_dtOld = dr_dt   ! This will be the same as dr_dt for the second half-step

!!******************************************************************************
!!Second "half-step" of the evolution loop
!!******************************************************************************

     call Timers_start("cosmology")
     call Cosmology_solveFriedmannEqn(dr_simTime, dr_dt)
     call Timers_stop("cosmology")

     if (dr_simTime .gt. tinitial + sim_tRelax) then
         grv_mode = 1
         if (dr_simTime .lt. tinitial + sim_tRelax) then
             orb_t = tinitial + sim_tRelax
             orb_dt = dr_simTime + dr_dt - tinitial - sim_tRelax
         else
             orb_t = dr_simTime
             orb_dt = dr_dt
         endif
         call Orbit_update
     endif

     dr_simTime = dr_simTime + dr_dt

     dr_simGeneration = 0

     call Timers_start("hydro")
     call Hydro(dr_myPE, dr_numProcs, blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, dr_rSweepDir)
     call Timers_stop("hydro")
  
     call Timers_start("sourceTerms")
     call Driver_sourceTerms(blockCount, blockList, dr_dt)
     call Timers_stop("sourceTerms")

     call Timers_start("Particles_advance")
     call Particles_advance(dr_dt, dr_dt)
     call Timers_stop("Particles_advance")
     
     call Gravity_potentialListOfBlocks(blockCount,blockList)

     call Bound_mass(blockCount, blockList)
     if (dr_simTime - dr_dt .gt. tinitial + sim_tRelax) then
         grv_mode = 2
         if (dr_simTime - dr_dt .lt. tinitial + sim_tRelax) then
             orb_t = tinitial + sim_tRelax
             orb_dt = dr_simTime - tinitial - sim_tRelax
         else
             orb_t = dr_simTime - dr_dt
             orb_dt = dr_dt
         endif
         call Orbit_update
     endif

     call Timers_start("cosmology")
     call Cosmology_redshiftHydro(dr_myPE, dr_numprocs, blockCount, blockList)
     call Timers_stop("cosmology")

     !--------------------------------------------------------------------
     !- End Physics Sequence -- Start Simulation Bookkeeping
     !--------------------------------------------------------------------

     !output a plotfile before the grid changes
     call Timers_start("IO_output")
     call IO_output(dr_myPE, dr_numProcs, dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun, PLOTFILE_AND_PARTICLEFILE)
     call Timers_stop("IO_output")


     !!if (itemp_limit) .eq. 1) call Hydro_timstepPrecompute()

     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement(dr_myPE, dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")

     if (gridChanged) dr_simGeneration = dr_simGeneration + 1

     dr_dtOld = dr_dt                     ! backup needed old 
     ! calculate new
     
     call Timers_start("compute dt")
     call Driver_computeDt(dr_myPE, dr_numProcs,  &
                         dr_nbegin, dr_nstep, &
                         dr_simTime, dr_dtOld, dr_dtNew)
     call Timers_stop("compute dt")
     dr_dt = dr_dtNew                                        ! store new
     

     !!-----------------------------------------------------------------
     !! Output for current step in evolution
     !!-----------------------------------------------------------------

     call Timers_start("IO_output")
     call IO_output(dr_myPE, dr_numProcs, dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun, CHECKPOINT_FILE_ONLY)
     call Timers_stop("IO_output")

!!*****************************************************************************
!!  Evolution Loop -- check termination conditions
!!*****************************************************************************


     !Exit if a .dump_restart or .kill was found during the last step
     if(endRun) exit

     !call Particles_dump(dr_myPE, blockCount, blockList, dr_nstep, dr_simTime, dr_dt)


     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift  
     !!        (also called redshiftFinal) 
     !!  (iii) the wall clock time is greater than the maximum 
     !!        (wall_clock_time_max)

     !!Update redshift from Driver's POV.  Need this for exit condition. -PR
     !!old redshift needed for accurate restarts.
     dr_redshiftOld = dr_redshift
     call Cosmology_getRedshift(dr_redshift)
     
     if (dr_simTime >= dr_tmax) then
        if(dr_myPE == MASTER_PE) then
           print *, "exiting: reached max SimTime"
        end if
        exit
     end if
     
     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_elapsedWCTime >  dr_wallClockTimeLimit) then
        if(dr_myPE == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        end if
        exit
     end if

     if (dr_redshift < dr_redshiftfinal .and. dr_useRedshift) then
        if(dr_mype == MASTER_PE) then
           print *, "exiting: reached redshiftfinal"
        end if
        exit
     end if

  enddo
  !The value of dr_nstep after the loop is (dr_nend + 1) if the loop iterated for
  !the maximum number of times.  However, we need to retain the value that
  !dr_nstep had during the last loop iteration, otherwise the number for nstep
  !that will be stored in a final checkpoint file will be wrong.
  dr_nstep = min(dr_nstep,dr_nend)

!!******************************************************************************
!! End of Evolution Loop
!!******************************************************************************



  call Timers_stop("evolution")

  call Logfile_stamp(dr_myPE, 'Exiting evolution loop' , '[Driver_evolveFlash]')

  !if a file termination, this may already be done.
  if(.NOT.endRun) call IO_outputFinal(dr_myPE, dr_numProcs)

  call Timers_getSummary(dr_myPE, max(0,dr_nstep-dr_nbegin+1))


  call Logfile_stamp(dr_myPE, "FLASH run complete.", "LOGFILE_END")

  call Logfile_close(dr_myPE)


  return
  
end subroutine Driver_evolveFlash



