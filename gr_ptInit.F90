!!****if* source/Grid/GridParticles/gr_ptInit
!!
!! NAME
!!
!!  gr_ptInit
!!
!! SYNOPSIS
!!
!!  gr_ptInit()
!!
!! DESCRIPTION
!!
!!  Initialize values for all data in the module gr_ptData,
!!  and allocate the scratch buffers
!!
!! ARGUMENTS
!!
!!
!!***

#include "Flash.h"
#ifdef FLASH_EDEP
#include "EnergyDeposition.h"
#endif

subroutine gr_ptInit()
  use gr_ptData, ONLY :   gr_ptDestBuf,gr_ptSourceBuf,gr_ptMaxPerProc,&
                          gr_ptRemove,gr_ptRemoveAlgo,&
                          gr_ptNumToReduce,gr_ptSieveCheckFreq,&
                          gr_ptLogLevel,gr_ptKeepLostParticles,gr_ptMaxVirtualCount
#ifndef FLASH_GRID_UG
  use gr_ptData, ONLY : gr_ptMaxPerProcUpperThresh, gr_ptMaxPerProcLowerThresh, &
       gr_ptMaxPerProcBlockFactor,gr_ptMaxPerProcBlockNoFuzz, gr_ptRefineOnPtMaxPerProc
  use Grid_data, ONLY : gr_meshMe, gr_refineOnParticleCount
#endif  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none 
  integer, save :: maxPerProc, propCount
 

  call RuntimeParameters_get("gr_ptRemove",gr_ptRemove)
  call RuntimeParameters_get("gr_ptRemoveAlgo",gr_ptRemoveAlgo)
  call RuntimeParameters_get("gr_ptNumToReduce",gr_ptNumToReduce)
  call RuntimeParameters_get("gr_ptSieveCheckFreq",gr_ptSieveCheckFreq)
  call RuntimeParameters_get("keepLostParticles",gr_ptKeepLostParticles)

!  call RuntimeParameters_get("pt_logLevel",gr_ptLogLevel)


  gr_ptMaxPerProc=1
  propCount=1

#ifdef NPART_PROPS
  if (NPART_PROPS > 1) then
     call RuntimeParameters_get("pt_maxPerProc",maxPerProc)
     propCount=max(NPART_PROPS,propCount)
     gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)
  endif
#endif
  gr_ptMaxVirtualCount=gr_ptMaxPerProc

#ifdef RAY_ATTR_COUNT
  call RuntimeParameters_get("ed_maxRayCount",maxPerProc)
  propCount=max(RAY_ATTR_COUNT,propCount)
  gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)
#endif
  ! JFG
  !gr_ptKeepLostParticles=.false.

  ! These are deallocated in gr_ptFinalize
  allocate(gr_ptDestBuf(propCount,gr_ptMaxPerProc))
  allocate(gr_ptSourceBuf(propCount,gr_ptMaxPerProc))

#ifndef FLASH_GRID_UG
  ! additional RPs for maxPerProc-based refinement criteria - KW
  call RuntimeParameters_get("gr_ptMaxPerProcUpperThresh",gr_ptMaxPerProcUpperThresh)
  call RuntimeParameters_get("gr_ptMaxPerProcLowerThresh",gr_ptMaxPerProcLowerThresh)
  call RuntimeParameters_get("gr_ptMaxPerProcBlockFactor",gr_ptMaxPerProcBlockFactor)
  call RuntimeParameters_get("gr_ptMaxPerProcBlockNoFuzz",gr_ptMaxPerProcBlockNoFuzz)
  call RuntimeParameters_get("gr_ptRefineOnPtMaxPerProc", gr_ptRefineOnPtMaxPerProc)
  if (gr_ptRefineOnPtMaxPerProc .AND. .NOT. gr_refineOnParticleCount) then
     call Logfile_stampMessage( &
          "WARNING: Ignoring gr_ptRefineOnPtMaxPerProc because RP refine_on_particle_count is .FALSE.")
     gr_ptRefineOnPtMaxPerProc = .FALSE.
  end if
#endif
  return
end subroutine gr_ptInit
