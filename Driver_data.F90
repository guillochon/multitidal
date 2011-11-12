!!****ih* source/Driver/DriverMain/Driver_data
!!
!! NAME
!!  Driver_data
!!
!! SYNOPSIS
!!
!!  use Driver_data
!!
!! DESCRIPTION 
!!  Driver_data is a fortran module that holds variables with
!!  the Driver Unit scope.  All variables located in Driver_data are
!!  accessible to subroutines in the Driver unit only.  (Note this 
!!  is a convention, there is nothing in the fortran language that
!!  prevents another routine from using Driver_data.  It is the FLASH3
!!  architecture that tries to enforce these rules.)
!!  
!!  All variables located in the Driver_data fortran module start with 
!! "dr_".  This is to indicate where they come from, and to make it easier
!! on the developer and user to see which variables are local variables and 
!! which belong to the Driver unit.
!!
!!
!! NOTES
!!  For those familiar with FLASH2:  The fortran data modules in FLASH3 
!!  (like this Driver_data) take the place of the database in FLASH2.
!!  In FLASH2 all the data was centrally located in the database.  In 
!!  FLASH3 each unit stores its own data.  When another unit needs to 
!!  get data from a different unit, accessor functions are used.  For example
!!  simulation time is stored in the Driver owned variable dr_simTime.  The
!!  IO unit queries the Driver unit with the accessor function
!!  Driver_getSimTime() to get the simTime and to check to see if a checkpoint
!!  needs to be written.
!! 
!!***

module Driver_data

#include "constants.h"


  integer, save :: dr_nBegin, dr_nStep, dr_nEnd
  real , save ::  dr_tmax      
  real , save ::  dr_dtInit, dr_dtMin, dr_dtMax
  real, save :: dr_tstepChangeFactor, dr_tempFactor
  real, save :: dr_wallClockTimeLimit
  real, save :: dr_redshiftMin, dr_redshift, dr_redshiftOld
  logical , save ::  dr_restart,dr_particlesInitialized
  logical , save ::  dr_eachProcWritesOwnAbortLog = .FALSE.
  integer, save :: dr_abortPause = 0
  
  character(len=40), save ::  dr_timeStamp

  !! Driver variables
  integer, save :: dr_myPE, dr_numProcs  
  integer, save :: dr_simGeneration
  real, save :: dr_dtOld, dr_dtNew, dr_simTime
  real, save :: dr_dt, dr_initialSimTime

  !! wall clock time variables
  real, save :: dr_elapsedWCTime, dr_initialWCTime
  real, save :: dr_redshiftInitial, dr_redshiftFinal
  logical, save :: dr_useRedshift  
  integer, save :: dr_fSweepDir, dr_rSweepDir

end module Driver_data

