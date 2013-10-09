!!****h* source/Simulation/SimulationMain/MultiTidal
!!
!! This is the header file for the Multitidal module
!! that defines its public interfaces.
!!***
Module Multitidal_interface
  implicit none
#include "constants.h"
  interface
     subroutine Multitidal_findExtrema(blockID,ivar,flag,extrema)
       integer,intent(IN):: blockID,ivar,flag
       real,intent(INOUT):: extrema
     end subroutine Multitidal_findExtrema
  end interface

end Module Multitidal_interface

