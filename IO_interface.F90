!!****h* source/IO/IO_interface
!!
!! This is the header file for the IO module that defines its
!! public interfaces.
!!
!!***

Module IO_interface
#ifndef IO_H
#define IO_H


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
     subroutine IO_init(myPE, numProcs)
       integer, intent(in) :: myPE
       integer, intent(in) :: numProcs
     end subroutine IO_init
  end interface

  interface
     subroutine IO_finalize()
     end subroutine IO_finalize
  end interface

  interface
     subroutine IO_output(myPE, numProcs, simTime, dt, nstep, nbegin, endRun, outputType)
       real, intent(in) :: simTime, dt
       integer, intent(in) :: nstep, nbegin
       integer, intent(in) :: myPE, numProcs
       logical, intent(out) :: endRun
       integer, intent(in), optional :: outputType
     end subroutine IO_output
  end interface

  interface 
     subroutine IO_outputFinal(myPE, numProcs)
       integer, intent(in) :: myPE, numProcs
     end subroutine IO_outputFinal
  end interface

  interface
     subroutine IO_outputInitial(myPE, numProcs, nbegin, initialSimTime)
       integer, intent(in) :: myPE, numProcs, nbegin
       real, intent(in)    :: initialSimTime
     end subroutine IO_outputInitial
  end interface

  interface
     subroutine IO_readCheckpoint(myPE, numProcs)
       integer, intent(in) :: myPE, numProcs
     end subroutine IO_readCheckpoint
  end interface
     
  interface 
     subroutine IO_readParticles(myPE, numProcs)
       integer, intent(in) :: myPE, numProcs
     end subroutine IO_readParticles
  end interface

  interface
     subroutine IO_readUserArray (myPE)
       integer, intent(in) :: myPE
     end subroutine IO_readUserArray
  end interface

  interface
     subroutine IO_sendOutputData()
     end subroutine IO_sendOutputData
  end interface

  interface
     subroutine IO_updateScalars(myPE, numProcs)
       integer, intent(in) :: myPE, numProcs
     end subroutine IO_updateScalars
  end interface

  interface
     subroutine IO_writeCheckpoint(myPE, numProcs)
       integer, intent(in) :: myPE, numProcs
     end subroutine IO_writeCheckpoint
  end interface

  interface
     subroutine IO_writeIntegralQuantities (myPE, isfirst, simTime)
       integer, intent(in) :: myPE, isfirst
       real, intent(in) :: simTime
     end subroutine IO_writeIntegralQuantities
  end interface

  interface
     subroutine IO_writeOrbitInfo (myPE, isfirst, simTime)
       integer, intent(in) :: myPE, isfirst
       real, intent(in) :: simTime
     end subroutine IO_writeOrbitInfo
  end interface

  interface
     subroutine IO_writeParticles(myPE, numProcs, particlesToCheckpoint)
       integer, intent(in) :: myPE, numProcs
       logical, intent(in) :: particlesToCheckpoint
     end subroutine IO_writeParticles
  end interface

  interface
     subroutine IO_writePlotfile(myPE, numProcs, forced)
       integer, intent(in) :: myPE, numProcs
       logical, optional, intent(in) :: forced
     end subroutine IO_writePlotfile
  end interface
  
  interface
     subroutine IO_writeUserArray (myPE)
       integer, intent(in) :: myPE
     end subroutine IO_writeUserArray
  end interface

  interface
     subroutine IO_initRPsFromCheckpoint (myPE, filename, ierr)
       integer, intent(in) :: myPE
       character(len=*),intent(IN) :: filename
       integer,intent(OUT) :: ierr
     end subroutine IO_initRPsFromCheckpoint
  end interface


end Module IO_interface


     
