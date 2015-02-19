!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxByTime,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK,&
                        gr_eosModeNow, gr_maxRefine, gr_globalComm
  use tree, ONLY : newchild, refine, derefine, stay, lnblocks
!!$  use physicaldata, ONLY : force_consistency
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Grid_interface, ONLY : Grid_fillGuardCells
  use Particles_interface, only: Particles_sinkMarkRefineDerefine
  ! Added by JFG
  use tree, ONLY: lrefine_max 
  use Grid_interface, ONLY: Grid_markRefineSpecialized, Grid_getMinCellSize
  use Driver_data, ONLY: dr_simTime
  use Simulation_data, ONLY: sim_objRadius, &
      sim_fluffRefineCutoff, sim_fluffDampCutoff, sim_cylinderRadius, &
      sim_xCenter, sim_yCenter, sim_zCenter, sim_kind, stvec, &
      sim_softenRadius, sim_fixedPartTag, sim_windNCells, &
      sim_tDelay, sim_periDist, sim_ptMass, sim_cylinderType, &
      sim_maxBlocks
  use pt_sinkInterface, ONLY : pt_sinkGatherGlobal
  use Multitidal_interface, ONLY : Multitidal_findExtrema
  use Particles_sinkData, ONLY : localnpf, particles_global
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  ! End JFG
  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref,max_blocks,ierr
  logical,save :: gcMaskArgsLogged = .FALSE.
  integer,save :: eosModeLast = 0
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask

  ! Added by JFG
  real,dimension(7) :: specs
  real,dimension(3) :: pvec
  real :: mcs, flow_dist, newton, max_dens

  call PhysicalConstants_get("Newton", newton)
  call Grid_getMinCellSize(mcs)
  ! End JFG


  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if
  
  if(gr_lrefineMaxByTime) then
     call gr_setMaxRefineByTime()
  end if

  call MPI_ALLREDUCE(lnblocks,max_blocks,1,MPI_INTEGER,MPI_SUM,gr_globalComm,ierr)
  if (max_blocks .gt. sim_maxBlocks) then
     gr_maxRefine = gr_maxRefine - 1
  end if

  if (gr_eosModeNow .NE. eosModeLast) then
     gcMaskArgsLogged = .FALSE.
     eosModeLast = gr_eosModeNow
  end if

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.
!!$  gcMask(NUNK_VARS+1:maskSize) = .TRUE.


  if (.NOT.gcMaskArgsLogged) then
     call Logfile_stampVarMask(gcMask, .true., '[Grid_markRefineDerefine]', 'gcArgs')
  end if

!!$  force_consistency = .FALSE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskArgsLogged,&
       selectBlockType=ACTIVE_BLKS)
     gcMaskArgsLogged = .TRUE.
!!$  force_consistency = .TRUE.

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! For PARAMESH2, call gr_markRefineDerefine here if it hasn't been called above.
  ! This is necessary to make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

  !! JFG
  if(gr_maxRefine .lt. lrefine_max) call gr_enforceMaxRefine(gr_maxRefine)

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)
  
  call Particles_sinkMarkRefineDerefine()

  if (dr_simTime .eq. 0.d0) then
      if (sim_kind .eq. 'polytrope') then
          specs = (/ sim_xCenter, sim_yCenter, sim_zCenter, sim_objRadius, 0., 0., 0. /)
          call Grid_markRefineSpecialized(INRADIUS, 4, specs(1:4), gr_maxRefine)
      elseif (sim_kind .eq. 'powerlaw') then
          specs = (/ sim_xCenter, sim_yCenter, sim_zCenter, sim_softenRadius, 0., 0., 0. /)
          call Grid_markRefineSpecialized(INRADIUS, 4, specs(1:4), gr_maxRefine)
          specs = (/ sim_xCenter + stvec(1), sim_yCenter + stvec(2), sim_zCenter + stvec(3), &
                     sim_objRadius, 0., 0., 0. /)
          call Grid_markRefineSpecialized(INRADIUS, 4, specs(1:4), gr_maxRefine)
      endif
  else
      if (sim_kind .eq. 'polytrope') then
          call Multitidal_findExtrema(DENS_VAR, 1, max_dens)
          specs = (/ real(DENS_VAR), 0., 1., 0., 0., 0., 0. /) ! First mark everything for derefine
          call Grid_markRefineSpecialized(THRESHOLD, 3, specs(1:3), -1)
          specs = (/ real(DENS_VAR), sim_fluffRefineCutoff*max_dens, -1., 0., 0., 0., 0. /) ! Then mark stuff that satisfies threshold for refine
          call Grid_markRefineSpecialized(THRESHOLD, 3, specs(1:3), gr_maxRefine)
      else
          specs = (/ real(DENS_VAR), 0., 1., 0., 0., 0., 0. /) ! First mark everything for derefine
          call Grid_markRefineSpecialized(THRESHOLD, 3, specs(1:3), -1)
          specs = (/ real(DENS_VAR), sim_fluffRefineCutoff, -1., 0., 0., 0., 0. /) ! Then mark stuff that satisfies threshold for refine
          call Grid_markRefineSpecialized(THRESHOLD, 3, specs(1:3), gr_maxRefine)
      endif
  endif

  if (sim_kind .eq. 'cylinder') then
      if (sim_cylinderType .eq. 1) then
          specs = (/ sim_xCenter + stvec(1) - sim_cylinderRadius, &
                     sim_xCenter + stvec(1) + sim_cylinderRadius, &
                     sim_yCenter + stvec(2) - 2.d0*mcs, &
                     sim_yCenter + stvec(2) + 2.d0*mcs, &
                     sim_zCenter + stvec(3) - sim_cylinderRadius, &
                     sim_zCenter + stvec(3) + sim_cylinderRadius, 0. /)
      else
          flow_dist = sim_periDist*((dr_simTime+sim_tDelay)/&
                      (2.d0*PI*dsqrt(sim_periDist**3.d0/(newton*sim_ptMass))))**(2.d0/3.d0)
          specs = (/ sim_xCenter + flow_dist - 2.d0*mcs, &
                     sim_xCenter + flow_dist + 2.d0*mcs, &
                     sim_yCenter - sim_cylinderRadius, &
                     sim_yCenter + sim_cylinderRadius, &
                     sim_zCenter - sim_cylinderRadius, &
                     sim_zCenter + sim_cylinderRadius, 0. /)
      endif

      call Grid_markRefineSpecialized(RECTANGLE, 7, specs(1:7), gr_maxRefine)
  elseif (sim_kind .eq. 'wind') then
      call pt_sinkGatherGlobal()
      do i = 1, localnpf
          if (idnint(particles_global(TAG_PART_PROP,i)) .ne. sim_fixedPartTag) then
              pvec(1:3) = particles_global(POSX_PART_PROP:POSZ_PART_PROP,i)
          endif
      enddo

      specs = (/ sim_xCenter, sim_yCenter, sim_zCenter, sim_softenRadius, 0., 0., 0. /)
      call Grid_markRefineSpecialized(INRADIUS, 4, specs(1:4), gr_maxRefine)
      specs = (/ pvec(1), pvec(2), pvec(3), sim_windNCells*mcs, 0., 0., 0. /)
      call Grid_markRefineSpecialized(INRADIUS, 4, specs(1:4), gr_maxRefine)
  elseif (sim_kind .eq. 'polytrope') then
      specs = (/ sim_xCenter, sim_yCenter, sim_zCenter, sim_objRadius, 0., 0., 0. /)
      call Grid_markRefineSpecialized(INRADIUS, 4, specs(1:4), gr_maxRefine)
  endif
  
  return
end subroutine Grid_markRefineDerefine

