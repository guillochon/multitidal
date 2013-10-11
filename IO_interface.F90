!!****h* source/IO/IO_interface
!!
!! This is the header file for the IO module that defines its
!! public interfaces.
!!
!!***

Module IO_interface
#ifndef IO_H
#define IO_H


  ! Added by JFG
  interface IO_writeOrbitInfo
     subroutine IO_writeOrbitInfo (isfirst, simTime)
       integer, intent(in) :: isfirst
       real, intent(in) :: simTime
     end subroutine IO_writeOrbitInfo
  end interface
  ! End JFG

  interface IO_setScalar
     subroutine IO_setScalarReal (name, value)
       character(len=*),intent(in)                :: name
       real,intent(in)                            :: value
     end subroutine IO_setScalarReal
     
     subroutine IO_setScalarInt(name, value)
       character(len=*),intent(in)                :: name
       integer,intent(in)                         :: value
     end subroutine IO_setScalarInt
     
     subroutine IO_setScalarStr (name, value)
       character(len=*),intent(in)                :: name, value
     end subroutine IO_setScalarStr
     
     subroutine IO_setScalarLog (name, value)
       character(len=*),intent(in)                 :: name
       logical,intent(in)                          :: value
     end subroutine IO_setScalarLog
     
  end interface
  
  interface IO_getScalar
     
     subroutine IO_getScalarReal (name, value)
       character(len=*),intent(in)                :: name
       real,intent(out)                           :: value
     end subroutine IO_getScalarReal
     
     subroutine IO_getScalarInt(name, value)
       character(len=*),intent(in)                :: name
       integer,intent(out)                        :: value
     end subroutine IO_getScalarInt
     
     subroutine IO_getScalarStr (name, value)
       character(len=*),intent(in)                :: name
       character(len=*),intent(out)               :: value
     end subroutine IO_getScalarStr
     
     subroutine IO_getScalarLog (name, value)
       character(len=*),intent(in)                 :: name
       logical,intent(out)                         :: value
     end subroutine IO_getScalarLog
     
  end interface
  
  interface IO_getPrevScalar
     
     subroutine IO_getPrevScalarReal (name, value, error)
       character(len=*),intent(in)                :: name
       real,intent(out)                           :: value
       integer, intent(out),optional                       :: error
     end subroutine IO_getPrevScalarReal
     
     subroutine IO_getPrevScalarInt(name, value, error)
       character(len=*),intent(in)                :: name
       integer,intent(out)                        :: value
       integer, intent(out),optional                       :: error
     end subroutine IO_getPrevScalarInt
     
     subroutine IO_getPrevScalarStr (name, value, error)
       character(len=*),intent(in)                :: name
       character(len=*),intent(out)               :: value
       integer, intent(out), optional             :: error
     end subroutine IO_getPrevScalarStr
     
     subroutine IO_getPrevScalarLog (name, value, error)
       character(len=*),intent(in)                 :: name
       logical,intent(out)                         :: value
       integer, intent(out),optional                       :: error
     end subroutine IO_getPrevScalarLog
     
  end interface
#endif

  interface
     subroutine IO_init()
     end subroutine IO_init
  end interface

  interface
     subroutine IO_finalize()
     end subroutine IO_finalize
  end interface

  interface
     subroutine IO_output( simTime, dt, nstep, nbegin, endRun, outputType)
       real, intent(in) :: simTime, dt
       integer, intent(in) :: nstep, nbegin
       logical, intent(out) :: endRun
       integer, intent(in), optional :: outputType
     end subroutine IO_output
  end interface

  interface 
     subroutine IO_outputFinal()
     end subroutine IO_outputFinal
  end interface

  interface
     subroutine IO_outputInitial( nbegin, initialSimTime)
       integer, intent(in) ::  nbegin
       real, intent(in)    :: initialSimTime
     end subroutine IO_outputInitial
  end interface

  interface
     subroutine IO_readCheckpoint()
     end subroutine IO_readCheckpoint
  end interface
     
  interface 
     subroutine IO_readParticles()
     end subroutine IO_readParticles
  end interface

  interface
     subroutine IO_readUserArray ()
     end subroutine IO_readUserArray
  end interface

  interface
     subroutine IO_sendOutputData()
     end subroutine IO_sendOutputData
  end interface

  interface
     subroutine IO_updateScalars()
     end subroutine IO_updateScalars
  end interface

  interface
     subroutine IO_writeCheckpoint()
     end subroutine IO_writeCheckpoint
  end interface

  interface
     subroutine IO_writeIntegralQuantities ( isfirst, simTime)
       integer, intent(in) ::  isfirst
       real, intent(in) :: simTime
     end subroutine IO_writeIntegralQuantities
  end interface

  interface
     subroutine IO_writeParticles( particlesToCheckpoint)
       logical, intent(in) :: particlesToCheckpoint
     end subroutine IO_writeParticles
  end interface

  interface
     subroutine IO_writePlotfile( forced)
       logical, optional, intent(in) :: forced
     end subroutine IO_writePlotfile
  end interface
  
  interface
     subroutine IO_writeUserArray ()
     end subroutine IO_writeUserArray
  end interface

  interface
     subroutine IO_initRPsFromCheckpoint ( filename, ierr)
       character(len=*),intent(IN) :: filename
       integer,intent(OUT) :: ierr
     end subroutine IO_initRPsFromCheckpoint
  end interface

  interface
     subroutine IO_checkForPlot(wrotePlot)
       implicit none
       logical, intent(out) :: wrotePlot
     end subroutine IO_checkForPlot
  end interface

  interface 
     subroutine IO_startRayWrite()
       implicit none
     end subroutine IO_startRayWrite
  end interface
  
  interface
     subroutine IO_writeRays(numRays, rayTags, posBuffer, powerBuffer, numPos)
       implicit none
       
       integer, intent(in) :: numRays
       integer, intent(in) :: rayTags(:)
       real, intent(in) :: posBuffer(:,:,:)
       real, intent(in) :: powerBuffer(:,:)
       integer, intent(in) :: numPos(:)
       
     end subroutine IO_writeRays
  end interface

  interface 
     subroutine IO_endRayWrite()
       implicit none
     end subroutine IO_endRayWrite
  end interface

  
end Module IO_interface


     
