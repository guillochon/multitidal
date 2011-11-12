!!****if* source/Driver/DriverMain/Driver_init
!!
!! NAME
!!  Driver_init
!!
!! SYNOPSIS
!!  Driver_init(integer (IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  Perform the Driver unit initializations.
!!  Gets runtime parameters from the flash.par, or from checkpoint file
!!  if run is a restart.
!!  
!!  Initializes the simulation time, initial dt, timestep, whether 
!!  from scratch or restart.
!!  Also, verifies that the initial dt passes CFL criteria.
!!
!! ARGUMENTS
!!
!!  myPE - current processor
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used by the basic Driver unit.
!!   Your specific implementation may have more runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You may have overwritten these values with the flash.par values
!!   for your specific run.  
!!
!!    dtinit [REAL]
!!        Initial timestep
!!    dtmax [REAL]
!!        Maximum timestep
!!    dtmin [REAL]
!!        Minimum timestep
!!    nbegin [INTEGER]
!!        First timestep
!!    nend [INTEGER]
!!        Maximum number of timesteps to take
!!    restart [BOOLEAN]
!!        Is this a restart run?
!!    tinitial [REAL]
!!        Initial simulation time
!!    tmax [REAL]
!!        Maximum simulation time
!!    wall_clock_time_limit [REAL]
!!        Total wall clock time limit (seconds)
!!    tstep_change_factor {REAL]
!!        factor to allow multiplicative increase in dt until it
!!        hits the CFL condition. This lets initial dt be
!!        very conservative initially, but increase rapidly to find the
!!        the optimum value.
!!    zInitial  [REAL]
!!        The initial redshift value.  < 0 if redshift is not being used.
!!    zFinal    [REAL]
!!        The final redshift value.  If reached, the simulation ends.
!!
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_myPE or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_init(myPE)
  use Driver_data
  use Driver_interface, ONLY:  Driver_putTimeStamp, Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use IO_interface, ONLY : IO_getScalar  
  use Cosmology_interface, ONLY : Cosmology_redshiftToTime

  implicit none       

#include "Flash.h"
#include "constants.h" 
  
  integer, intent(in) :: myPE
  character(len=40), save                :: dateStr
  integer, save :: sweepOrder


  ! get the parameters needed by Driver
  call RuntimeParameters_get("eachProcWritesOwnAbortLog", dr_eachProcWritesOwnAbortLog)
  call RuntimeParameters_get("dr_abortPause", dr_abortPause)
  !! hand dr_abortPause out to C routines to avoid architecture-dependent code
  call driver_abortflashc_set_pause(dr_abortPause)

  call RuntimeParameters_get("nend", dr_nend)

  call RuntimeParameters_get("restart", dr_restart)
  call RuntimeParameters_get("tmax", dr_tmax)
  call RuntimeParameters_get("tinitial",dr_initialSimTime)

  call RuntimeParameters_get("dtinit",dr_dtInit)
  call RuntimeParameters_get("dtmin",dr_dtMin)
  call RuntimeParameters_get("dtmax",dr_dtMax)
  call RuntimeParameters_get("tstep_change_factor", dr_tstepChangeFactor)
  call RuntimeParameters_get("wall_clock_time_limit",dr_wallClockTimeLimit)
  call RuntimeParameters_get("zInitial", dr_redshiftInitial)
  call RuntimeParameters_get("zFinal", dr_redshiftFinal)
  call RuntimeParameters_get("sweepOrder",sweepOrder)

  !Check to see if redshift is being used.
  dr_useRedshift = (dr_redshiftInitial >= 0.0)
  
!!  This is not yet in use. This parameter controls timestep based
!!  upon the temperature.
!!  call RuntimeParameters_get("temp_factor", dr_tempFactor)

!! initialize a timeStamp for the run; needed in case Logfile is not included
  call current_date_time(dateStr)
  call Driver_putTimeStamp(dateStr)

  if(sweepOrder == 123) then
     dr_fSweepDir=SWEEP_XYZ
     dr_rSWeepDir=SWEEP_ZYX
  elseif(sweepOrder == 132) then
     dr_fSweepDir=SWEEP_XZY
     dr_rSWeepDir=SWEEP_YZX
  elseif(sweepOrder == 213) then
     dr_fSweepDir=SWEEP_YXZ
     dr_rSWeepDir=SWEEP_ZXY
  elseif(sweepOrder == 231) then
     dr_fSweepDir=SWEEP_YZX
     dr_rSWeepDir=SWEEP_XZY
  elseif(sweepOrder == 312) then
     dr_fSweepDir=SWEEP_ZXY
     dr_rSWeepDir=SWEEP_YXZ
  elseif(sweepOrder == 321) then
     dr_fSweepDir=SWEEP_ZYX
     dr_rSWeepDir=SWEEP_XYZ
  else
     call Driver_abortFlash("Driver_init : invalid SweepOrder for Hydro")
  end if
  
  
  if (dr_restart) then
  
     !get values from the scalar list from the previous run
     call IO_getScalar("nstep", dr_nstep)
     call IO_getScalar("time", dr_simTime)
     call IO_getScalar("dt", dr_dt)
     call IO_getScalar("dtOld", dr_dtOld)
     if(dr_useRedshift) call IO_getScalar("redshift", dr_redshift)

     dr_nbegin = dr_nstep
     dr_initialSimTime = dr_simTime

  else if (.not. dr_restart) then
     call RuntimeParameters_get("nbegin", dr_nbegin)
     dr_simTime = dr_initialSimTime
     dr_nstep = dr_nbegin
     if(dr_useRedshift) then
        call Cosmology_redshiftToTime(dr_redshiftInitial, dr_simTime)
        dr_redshift = dr_redshiftInitial
        dr_initialSimTime = dr_simTime
     else
        dr_redshift = 0.0
     endif

  end if

  dr_simGeneration = 0
  dr_particlesInitialized=.false.
  return
end subroutine Driver_init
