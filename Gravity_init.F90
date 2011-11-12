!!****if* source/physics/Gravity/GravityMain/Poisson/Multipole/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!
!! 
!! SYNOPSIS
!!
!!  Gravity_init(integer(IN) :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!!  ARGUMENTS
!!
!!  myPE - local processor number
!!
!!***

subroutine Gravity_init(myPE)

  use Gravity_data
  use Simulation_data, ONLY: sim_objMass, obj_radius, obj_ipos, sim_tRelax
  use Grid_interface, ONLY: Grid_getNumProcs, Grid_getListOfBlocks
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Gravity_interface, ONLY: Gravity_potentialListOfBlocks
  use gr_mpoleData, ONLY: Mtot, Xcm, Ycm, Zcm, oXcm, oYcm, oZcm
  use Driver_data, ONLY: dr_restart
  use IO_interface, ONLY: IO_getScalar
  use tree, ONLY: lrefine_max

  implicit none
  double precision :: newton, newx, newy, newz, period, a, ecc_anom, start_dist
  integer :: i
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: myPE
  character(len=MAX_STRING_LENGTH) :: strGeometry
  character(len=100) :: logstr

  ! Everybody should know these
  grv_myPE = myPE
  call Grid_getNumProcs(grv_numProcs)


  call PhysicalConstants_get("Newton", newton)
  
  if (dr_restart) then
  ! NEED TO ADD OBJECT VECTORS HERE
      call IO_getScalar("Mtot", Mtot)
      call IO_getScalar("Xcm", Xcm)
      call IO_getScalar("Ycm", Ycm)
      call IO_getScalar("Zcm", Zcm)
      call IO_getScalar("oXcm", oXcm)
      call IO_getScalar("oYcm", oYcm)
      call IO_getScalar("oZcm", oZcm)
      call IO_getScalar("ptxpos", grv_ptvec(1))
      call IO_getScalar("ptypos", grv_ptvec(2))
      call IO_getScalar("ptzpos", grv_ptvec(3))
      call IO_getScalar("ptxvel", grv_ptvec(4))
      call IO_getScalar("ptyvel", grv_ptvec(5))
      call IO_getScalar("ptzvel", grv_ptvec(6))
      call IO_getScalar("obxpos", grv_obvec(1))
      call IO_getScalar("obypos", grv_obvec(2))
      call IO_getScalar("obzpos", grv_obvec(3))
      call IO_getScalar("obxvel", grv_obvec(4))
      call IO_getScalar("obyvel", grv_obvec(5))
      call IO_getScalar("obzvel", grv_obvec(6))
      call IO_getScalar("optxpos", grv_optvec(1))
      call IO_getScalar("optypos", grv_optvec(2))
      call IO_getScalar("optzpos", grv_optvec(3))
      call IO_getScalar("optxvel", grv_optvec(4))
      call IO_getScalar("optyvel", grv_optvec(5))
      call IO_getScalar("optzvel", grv_optvec(6))
      call IO_getScalar("oobxpos", grv_oobvec(1))
      call IO_getScalar("oobypos", grv_oobvec(2))
      call IO_getScalar("oobzpos", grv_oobvec(3))
      call IO_getScalar("oobxvel", grv_oobvec(4))
      call IO_getScalar("oobyvel", grv_oobvec(5))
      call IO_getScalar("oobzvel", grv_oobvec(6))
      call IO_getScalar("dynrefinemax", grv_dynRefineMax)
  endif

  call RuntimeParameters_get("geometry", strGeometry)
  call RuntimeParameters_mapStrToInt(strGeometry, grav_geometry)
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  call RuntimeParameters_get("ptmass", grv_ptmass)
  call RuntimeParameters_get("sim_periBeta", peri_beta)
  call RuntimeParameters_get("sim_startBeta", start_beta)
  call RuntimeParameters_get("sim_periodFac", period_fac)
  call RuntimeParameters_get("grv_comCutoff",grv_comCutoff)
  call RuntimeParameters_get("grv_cfl",grv_cfl)

  call RuntimeParameters_get("sim_periTime", peri_time)
  call RuntimeParameters_get("sim_orbEcc", orb_ecc)
  call RuntimeParameters_get("useGravity", useGravity)
  call RuntimeParameters_get("updateGravity", updateGravity)

  grav_poisfact = 4.d0 * PI * newton
  grv_factor = -newton * grv_ptmass
  grv_thresh = 1.d-10

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call Gravity_potentialListOfBlocks(blockCount,blockList)
  call Bound_mass(blockCount, blockList)
  if (peri_beta .eq. 0.d0) then
      call RuntimeParameters_get("sim_periDist", peri_dist)
  else
      peri_dist = obj_radius(obj_ipos)/peri_beta*(grv_ptmass/grv_bound)**(1.d0/3.d0)
  endif
  if (period_fac .gt. 0.d0) then
      a = peri_dist/(1.d0 - orb_ecc)
      period = sqrt(4.d0*PI**2.d0/newton/(grv_ptmass + grv_bound)*a**3.d0)
      peri_time = period_fac*period + sim_tRelax
  endif
  if (start_beta .gt. 0.d0) then
      a = peri_dist/(1.d0 - orb_ecc)
      period = sqrt(4.d0*PI**2.d0/newton/(grv_ptmass + grv_bound)*a**3.d0)
      start_dist = obj_radius(obj_ipos)/start_beta*(grv_ptmass/grv_bound)**(1.d0/3.d0)
      if (start_dist .gt. 2.d0*a - peri_dist) call Driver_abortFlash('start_dist too large!')
      ecc_anom = acos((a - start_dist)/a/orb_ecc)
      peri_time = abs(ecc_anom - orb_ecc*sin(ecc_anom))*period/2.d0/PI + sim_tRelax
  endif
  if (.not. dr_restart) then
      call calc_orbit(sim_tRelax, grv_bound, grv_ptmass, grv_obvec, grv_ptvec)
      grv_optvec = grv_ptvec
      grv_oobvec = grv_obvec
      grv_dynRefineMax = lrefine_max
  endif
  if (myPE .eq. MASTER_PE) then
      write(logstr, fmt='(A30, 2ES15.8)') 'Period, pericenter time:', period, peri_time
      call Logfile_stampMessage(myPE, logstr)
      write(logstr, fmt='(A30, 2ES15.8)') 'Obj. radius, pericenter dist:', obj_radius(obj_ipos), peri_dist
      call Logfile_stampMessage(myPE, logstr)
      !do i = 1, 30
      !    call calc_orbit(1000.d0*i, newx, newy, newz, period)
      !    write(logstr, fmt='(2ES11.4)'), newx, newy
      !    call Logfile_stampMessage(myPE, logstr)
      !enddo
  endif

  return
end subroutine Gravity_init
