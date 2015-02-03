!!****if* source/Simulation/SimulationMain/MultiTidalPoly/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId, 
!!                       integer :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up MultiTidalPoly.
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  myPE -          current processor number
!!
!!
!!***

subroutine Simulation_initBlock (blockId, myPE)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData, Grid_fillGuardCells, Grid_getDeltas
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

#ifdef LOADPROFILE
#include "Starprof.h"
#endif
  
  integer,intent(IN) ::  blockId
  integer,intent(IN) ::  myPE
  
  integer  ::  i, j, k, jLo, jHi
  integer  ::  ii, jj, kk, put
  real     ::  distInv, xDist, yDist, zDist, &
               bhxDist, bhyDist, bhzDist, bhDist
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac, softening_radius
  real     ::  x1,x2,x3,cos_ang,sin_ang,lambda,radius,rot
  real     ::  dx, dy, dz, magx, magy, magz, magp, divb
  real     ::  vx, vy, vz, p, rho, e, ek, t, mp, kb, newton
  real     ::  dist, gam, rho0, T0, rsc, rho0in, T0in
  integer  ::  istat

  real, allocatable,dimension(:) :: xCoord,xCoordL,xCoordR,&
                                    yCoord,yCoordL,yCoordR,&
                                    zCoord,zCoordL,zCoordR
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  character(4), save :: unitSystem
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
#ifdef FL_NON_PERMANENT_GUARDCELLS
  real, pointer, dimension(:,:,:,:) :: solnData
#endif
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Az,Ax,Ay
#else
  real, allocatable, dimension(:,:,:) :: Az,Ax,Ay
#endif

  logical :: gcell = .true., within_radius

#ifdef LOADPROFILE
  real, dimension(sim_tableCols) :: sumVars
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn = 0.0
#else
  real, dimension(2) :: sumVars
#endif

  call PhysicalConstants_get("proton mass", mp)
  call PhysicalConstants_get("Boltzmann", kb)
  call PhysicalConstants_get("Newton", newton)
  
  call RuntimeParameters_get("sink_softening_radius", softening_radius)
  call RuntimeParameters_get("UnitSystem", unitSystem)

  ! Anninos 2012
  if (sim_useRadialProfile) then
      rsc = 1.2d16
      rho0 = 1.3d-21*sim_xRayFraction
      T0 = 0.4d0*obj_mu*newton*sim_ptMass*mp/kb/rsc
      rho0in = rho0*(softening_radius/rsc)**(-1.5d0)
      T0in = 0.5d0*T0*(softening_radius/rsc)**(-1.d0)
  else
      rsc = 1.d0
      rho0 = 0.d0
      T0 = 0.d0
      rho0in = 0.d0
      T0in = 0.d0
  endif

  !do i = 1, np
  !   print *, obj_radius(i), obj_rhop(i), eint(i) 
  !enddo
  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX), stat=istat)
  allocate(xCoordL(sizeX),stat=istat)
  allocate(xCoordR(sizeX),stat=istat)

  allocate(yCoord(sizeY), stat=istat)
  allocate(yCoordL(sizeY),stat=istat)
  allocate(yCoordR(sizeY),stat=istat)

  allocate(zCoord(sizeZ), stat=istat)
  allocate(zCoordL(sizeZ),stat=istat)
  allocate(zCoordR(sizeZ),stat=istat)

  xCoord  = 0.0
  xCoordL = 0.0
  xCoordR = 0.0

  yCoord  = 0.0
  yCoordL = 0.0
  yCoordR = 0.0

  zCoord  = 0.0
  zCoordL = 0.0
  zCoordR = 0.0

#ifndef FIXEDBLOCKSIZE
  if (NDIM == 2) then
     allocate(Ax(sizeX+1,sizeY+1,1),stat=istat)
     allocate(Ay(sizeX+1,sizeY+1,1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,1),stat=istat)
  elseif (NDIM == 3) then
     allocate(Ax(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
     allocate(Ay(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
  endif
#endif

  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS,blockId,CENTER,    sim_gCell,zCoord, sizeZ)
     call Grid_getCellCoords(KAXIS,blockId,LEFT_EDGE, sim_gCell,zCoordL,sizeZ)
     call Grid_getCellCoords(KAXIS,blockId,RIGHT_EDGE,sim_gCell,zCoordR,sizeZ)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS,blockId,CENTER,    sim_gCell,yCoord, sizeY)
     call Grid_getCellCoords(JAXIS,blockId,LEFT_EDGE, sim_gCell,yCoordL,sizeY)
     call Grid_getCellCoords(JAXIS,blockId,RIGHT_EDGE,sim_gCell,yCoordR,sizeY)
  endif

  call Grid_getCellCoords(IAXIS,blockId,CENTER,    sim_gCell,xCoord, sizeX)
  call Grid_getCellCoords(IAXIS,blockId,LEFT_EDGE, sim_gCell,xCoordL,sizeX)
  call Grid_getCellCoords(IAXIS,blockId,RIGHT_EDGE,sim_gCell,xCoordR,sizeX)

  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     ! Find a real difference between z's if problem is >= 3D
     if (NDIM > 2) then
        if (k .eq. 1) then
           dzz = zCoord(2) - zCoord(1) 
        else
           dzz = zCoord(k) - zCoord(k-1) 
        endif
     ! Otherwise this problem is <= 2D, so dzz is meaningless
     else
       dzz = 0.0
     endif
     zz = zCoord(k)

     ! Always use first ptvecs element as particle for determining radial
     ! profile (if enabled)
     bhzDist = zCoord(k) - (sim_zCenter + ptvecs(1,3))
     
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a real difference between y's if problem is >= 2D
        if (NDIM > 1) then
           if (j .eq. 1) then
              dyy = yCoord(2) - yCoord(1) 
           else
              dyy = yCoord(j) - yCoord(j-1) 
           endif
        ! Otherwise this problem is <= 1D, so dyy is meaningless
        else
          dyy = 0.0
        endif
        yy = yCoord(j)
        
        bhyDist = yCoord(j) - (sim_yCenter + ptvecs(1,2))

        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
           if (i .eq. 1) then
              dxx = xCoord(2) - xCoord(1) 
           else
              dxx = xCoord(i) - xCoord(i-1) 
           endif
           
           sumVars = 0.

#ifdef LOADPROFILE
           xn = 0.
#endif
           
           !
           !       Break the cell into sim_nSubZones^NDIM sub-zones, and look up the
           !       appropriate quantities along the 1d profile for that subzone.  
           !
           !       Have the final values for the zone be equal to the average of
           !       the subzone values.
           ! 

           within_radius = .false.

           bhxDist = xCoord(i) - (sim_xCenter + ptvecs(1,1))

           bhDist = dsqrt(bhxDist**2 + bhyDist**2 + bhzDist**2)

           if (sim_kind .ne. 'cylinder' .and. sim_kind .ne. 'wind') then
              do kk = 1, sim_nSubZones
                 zz    = zCoord(k) + dzz*((2*kk - 1)*sim_inSubInv - 0.5)
                 zDist = zz - (sim_zCenter + stvec(3))
                 
                 do jj = 1, sim_nSubZones
                    yy    = yCoord(j) + dyy*((2*jj - 1)*sim_inSubInv - 0.5)
                    yDist = yy - (sim_yCenter + stvec(2))
                    
                    do ii = 1, sim_nSubZones
                       xx    = xCoord(i) + dxx*((2*ii - 1)*sim_inSubInv - 0.5)
                       xDist = xx - (sim_xCenter + stvec(1))
                       
                       dist = dsqrt( xDist**2 + yDist**2 + zDist**2 )

                       distInv = 1. / max( dist, 1.d-10 )
                       !
                       !  a point at `dist' is frac-way between jLo and jHi.   We do a
                       !  linear interpolation of the quantities at jLo and jHi and sum those.
                       ! 
#ifdef LOADPROFILE
                       call sim_find (sim_table(:,R_PROF), sim_tableRows, dist, jLo)
                       frac = 0.
                       if (jLo .eq. 0) then
                          jLo = 1
                          jHi = 1
                       else if (jLo .ge. sim_tableRows) then
                          jLo = sim_tableRows
                          jHi = sim_tableRows
                       else
                          jHi = jLo + 1
                          frac = (dist - sim_table(jLo,R_PROF)) / & 
                               (sim_table(jHi,R_PROF)-sim_table(jLo,R_PROF))
                       endif
                       sumVars = sumVars + sim_table(jLo,:) + frac*(sim_table(jHi,:) - sim_table(jLo,:))
#else
                       call sim_find (obj_radius, obj_ipos, dist, jLo)
                       if (jLo .le. obj_ipos) then
                           within_radius = .true.
                           if (jLo .eq. 0) then
                              jLo = 1
                              jHi = 1
                              frac = 0.
                           else if (jLo .eq. obj_ipos) then
                              jLo = obj_ipos
                              jHi = obj_ipos
                              frac = 0.
                           else
                              jHi = jLo + 1
                              frac = (dist - obj_radius(jLo)) / & 
                                   (obj_radius(jHi)-obj_radius(jLo))
                           endif
                           sumVars(1) = sumVars(1) + & 
                                obj_rhop(jLo) + frac*(obj_rhop(jHi) - obj_rhop(jLo))
                           sumVars(2) = sumVars(2) +  & 
                                obj_prss(jLo) + frac*(obj_prss(jHi) - obj_prss(jLo))
                       endif
#endif
                    enddo
                 enddo
              enddo
           endif
           
#ifdef LOADPROFILE
           sumVars = sumVars*sim_inszd
           xx = xCoord(i) - sim_xCenter
           yy = yCoord(j) - sim_yCenter
           zz = zCoord(k) - sim_zCenter
           dist = sqrt(xx**2 + yy**2 + zz**2)
           if (dist .le. sim_objRadius) then
              sumVars(RHO_PROF) = max (sumVars(RHO_PROF), sim_rhoAmbient)
#ifdef H1_SPEC
              xn(H1_SPEC) = sumVars(H1_PROF)
#endif
#ifdef HE3_SPEC
              xn(HE3_SPEC) = sumVars(HE3_PROF)
#endif
#ifdef HE4_SPEC
              xn(HE4_SPEC) = sumVars(HE4_PROF)
#endif
#ifdef C12_SPEC
              xn(C12_SPEC) = sumVars(C12_PROF)
#endif
#ifdef N14_SPEC
              xn(N14_SPEC) = sumVars(N14_PROF)
#endif
#ifdef O16_SPEC
              xn(O16_SPEC) = sumVars(O16_PROF)
#endif
#ifdef NE20_SPEC
              xn(NE20_SPEC) = sumVars(NE20_PROF)
#endif
#ifdef MG24_SPEC
              xn(MG24_SPEC) = sumVars(MG24_PROF)
#endif
#ifdef SI28_SPEC
              xn(SI28_SPEC) = sumVars(SI28_PROF)
#endif
              xn = xn / sum(xn)
           else
              sumVars(RHO_PROF) = sim_rhoAmbient
              sumVars(TEMP_PROF) = sim_tAmbient
              xn(H1_SPEC) = 1.0
           endif
#else
           vx  = 0.0d0
           vy  = 0.0d0
           vz  = 0.0d0

           if (within_radius) then
               rho = max (sumVars(1)*sim_inszd, sim_rhoAmbient)
               p   = max (sumVars(2)*sim_inszd, sim_pAmbient)
               t   = p/(rho/mp/obj_mu*kb)
               if (sim_fixedParticle .eq. 0 .or. sim_fixedParticle .eq. 1) then
                   vx = stvec(4)
                   vy = stvec(5)
                   vz = stvec(6)
               endif
           else
               if (bhDist .le. softening_radius) then
                   rho = max (rho0in*(bhDist/softening_radius)**3.d0, sim_rhoAmbient)
                   t   = max (T0in*(bhDist/softening_radius)**2.d0, sim_tAmbient)
               else
                   rho = max (rho0*(bhDist/rsc)**(-1.5d0), sim_rhoAmbient)
                   t   = max (T0*(bhDist/rsc)**(-1.d0), sim_tAmbient)
               endif
               ! From Anninos 2012
               !rho = max (1.3d-21*1.d-1*(bhDist/1.2d16)**(-1.125d0), sim_rhoAmbient)
               !t   = max (1.d8*(bhDist/1.2d16)**(-0.75d0), sim_tAmbient)
               p   = t*rho/mp/obj_mu*kb
               if (sim_fixedParticle .eq. 2) then
                   vx = ptvecs(1,4)
                   vy = ptvecs(1,5)
                   vz = ptvecs(1,6)
               endif
           endif

           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           !
           !  assume gamma-law equation of state
           !
           e   = p/(obj_gamc-1.)
           e   = e/rho + ek
           e   = max (e, sim_smallP)
#endif
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k

#ifdef FL_NON_PERMANENT_GUARDCELLS

#ifdef LOADPROFILE
           do put=1,NSPECIES
              if (xn(SPECIES_BEGIN+put-1) == 0.0) xn(SPECIES_BEGIN+put-1) = sim_smallX
              solnData(SPECIES_BEGIN+put-1,i,j,k)=xn(SPECIES_BEGIN+put-1)
           enddo

           solnData(DENS_VAR,i,j,k)=sumVars(RHO_PROF)
           solnData(TEMP_VAR,i,j,k)=sumVars(TEMP_PROF)
#else
           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
           solnData(TEMP_VAR,i,k,k)=t
           solnData(GAME_VAR,i,j,k)=obj_gamc
           solnData(GAMC_VAR,i,j,k)=obj_gamc
           ! Need to add species setting for non-permanent guard cells...
#endif
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz
#else
#ifdef LOADPROFILE
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, sumVars(RHO_PROF))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, sumVars(TEMP_PROF))    
           do put=1,NSPECIES
              if (xn(SPECIES_BEGIN+put-1) == 0.0) xn(SPECIES_BEGIN+put-1) = sim_smallX
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+put-1,&
                   EXTERIOR,axis,xn(SPECIES_BEGIN+put-1))
           enddo
#else
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, t)    
           do put=SPECIES_BEGIN,SPECIES_END
              if (obj_xn(put) == 0.0) obj_xn(put) = sim_smallX
              call Grid_putPointData(blockID,CENTER,put,&
                   EXTERIOR,axis,obj_xn(put))
           enddo
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, obj_gamc)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, obj_gamc)
#endif
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
#endif
        enddo
     enddo
  enddo

  ! JFG --- Now initialize the magnetic fields, based on magnetoHD setup

  !------------------------------------------------------------------------------
  ! Construct Az at each cell corner
  ! Bx = dAz/dy - dAy/dz
  ! By = dAx/dz - dAz/dx
  ! Bz = dAy/dx - dAx/dy
  Az = 0.
  Ax = 0.
  Ay = 0.

  rot = atan(sim_rx/sim_ry)
  cos_ang = cos(rot)
  sin_ang = sin(rot)

  if (cos_ang >= sin_ang) then
     lambda = (sim_xMax-sim_xMin)*cos_ang
  else
     lambda = (sim_zMax-sim_zMin)*sin_ang
  endif

  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)
  dz = del(3)


  if (NDIM == 2) then
     k = 1
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

#if NFACE_VARS > 1
           ! x Coord at cell corner
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              x1 = xCoordL(i)
           else
              x1 = xCoordR(i-1)
           endif

           ! y Coord at cell corner
           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              x2 = yCoordL(j)
           else
              x2 = yCoordR(j-1)
           endif
#else
           ! x Coord at cell center
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              x1 = xCoord(i)
           else
              x1 = xCoord(i-1) + dx
           endif

           ! y Coord at cell center
           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              x2 = yCoord(j)
           else
              x2 = yCoord(j-1) + dy
           endif
#endif
           ! define radius of the field loop
           radius = sqrt((x1-sim_xCenter)**2 + (x2-sim_yCenter)**2)

           if (radius <= sim_fieldLoopRadius) then
              Ax(i,j,k) = 0.
              Ay(i,j,k) = 0.
              Az(i,j,k) = sim_Az_initial*(sim_fieldLoopRadius - radius)
           else
              Ax(i,j,k) = 0.
              Ay(i,j,k) = 0.
              Az(i,j,k) = 0.
           endif
        enddo
     enddo
  elseif (NDIM == 3) then
     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

              !! Rotated coordinate system
              !! / x1 \    / cos  0  sin \
              !! | x2 |  = |  0   1   0  |
              !! \ x3 /    \ -sin 0  cos /

#if NFACE_VARS > 0
              ! x Coord at cell corner
              if (i <=blkLimitsGC(HIGH,IAXIS)) then
                 xx = xCoordL(i)
              else
                 xx = xCoordR(i-1)
              endif

              ! y Coord at cell corner
              if (j <=blkLimitsGC(HIGH,JAXIS)) then
                 yy = yCoordL(j)
              else
                 yy = yCoordR(j-1)
              endif

              ! z Coord at cell corner
              if (k <=blkLimitsGC(HIGH,KAXIS)) then
                 zz = zCoordL(k)
              else
                 zz = zCoordR(k-1)
              endif

#else
              ! x Coord at cell center
              if (i <=blkLimitsGC(HIGH,IAXIS)) then
                 xx = xCoord(i)
              else
                 xx = xCoord(i-1) + dx
              endif

              ! y Coord at cell center
              if (j <=blkLimitsGC(HIGH,JAXIS)) then
                 yy = yCoord(j)
              else
                 yy = yCoord(j-1) + dy
              endif

              ! z Coord at cell center
              if (k <=blkLimitsGC(HIGH,KAXIS)) then
                 zz = zCoord(k)
              else
                 zz = zCoord(k-1) + dz
              endif
#endif
              ! define radius of the field loop
              radius = sqrt((xx-sim_xCenter)**2 + (yy-sim_yCenter)**2)

              if (radius <= sim_fieldLoopRadius) then
                 Ax(i,j,k) = 0.
                 Ay(i,j,k) = 0.
                 Az(i,j,k) = sim_Az_initial*(sim_fieldLoopRadius - radius)
              else
                 Ax(i,j,k) = 0.
                 Ay(i,j,k) = 0.
                 Az(i,j,k) = 0.
              endif
              !! For Ax and Ay
              !x1 = (cos_ang*xCoord(i) + sin_ang*zz) ! with rotation
              !!x1 = xx !without any rotation
              !x2 = yy


              !do while (x1 > 0.5*lambda)
              !   x1 = x1 - lambda
              !enddo
              !do while (x1 < -0.5*lambda)
              !   x1 = x1 + lambda
              !enddo

              !radius = sqrt(x1**2 + x2**2)

              !if (radius < sim_fieldLoopRadius) then
              !   Ax(i,j,k) = sim_Az_initial*(sim_fieldLoopRadius - radius)*(-sin_ang)
              !   print *, 'Ax', Ax(i,j,k), radius, sim_fieldLoopRadius
              !endif
              !Ay(i,j,k) = 0.

              !! For Az
              !x1 = (cos_ang*xx + sin_ang*zCoord(k)) ! with rotation
              !!x1 = xx !without any rotation
              !x2 = yy


              !do while (x1 > 0.5*lambda)
              !   x1 = x1 - lambda
              !enddo
              !do while (x1 < -0.5*lambda)
              !   x1 = x1 + lambda
              !enddo

              !radius = sqrt(x1**2 + x2**2)

              !if (radius < sim_fieldLoopRadius) then
              !   Az(i,j,k) = sim_Az_initial*(sim_fieldLoopRadius - radius)*cos_ang
              !   print *, 'Az', Az(i,j,k), radius, sim_fieldLoopRadius
              !endif

           enddo
        enddo
     enddo
  endif

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif


  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
#if NFACE_VARS > 0
           !! In this case we initialized Az using the cell-cornered coordinates.
           if (sim_killdivb) then
              if (NDIM == 2) then
                 facexData(MAG_FACE_VAR,i,j,k)= (Az(i,j+1,k)-Az(i,j,k))/dy
                 faceyData(MAG_FACE_VAR,i,j,k)=-(Az(i+1,j,k)-Az(i,j,k))/dx
              elseif (NDIM == 3) then
                 facexData(MAG_FACE_VAR,i,j,k)= -(Ay(i,j,k+1)-Ay(i,j,k))/dz + (Az(i,j+1,k)-Az(i,j,k))/dy
                 faceyData(MAG_FACE_VAR,i,j,k)=  (Ax(i,j,k+1)-Ax(i,j,k))/dz - (Az(i+1,j,k)-Az(i,j,k))/dx
                 facezData(MAG_FACE_VAR,i,j,k)= -(Ax(i,j+1,k)-Ax(i,j,k))/dy + (Ay(i+1,j,k)-Ay(i,j,k))/dx
              endif
           endif
#else
           !! In this case we initialized Az using the cell-centered coordinates.
           if (NDIM == 2) then
              magx = 0.5*(Az(i,j+1,k)-Az(i,j-1,k))/dy
              magy =-0.5*(Az(i+1,j,k)-Az(i-1,j,k))/dx
           elseif (NDIM == 3) then
              magx = -0.5*((Ay(i,j,k+1)-Ay(i,j,k-1))/dz + (Az(i,j+1,k)-Az(i,j-1,k))/dy)
              magy =  0.5*((Ax(i,j,k+1)-Ax(i,j,k-1))/dz - (Az(i+1,j,k)-Az(i-1,j,k))/dx)
              magz = -0.5*((Ax(i,j+1,k)-Ax(i,j-1,k))/dy + (Ay(i+1,j,k)-Ay(i-1,j,k))/dx)
           endif
           if (unitSystem .eq. 'cgs') then
               magp = .5/(4.*PI)*dot_product((/ magx, magy, magz /),&
                                             (/ magx, magy, magz /))
           else
               magp = .5*dot_product((/ magx, magy, magz /),&
                                     (/ magx, magy, magz /))
           endif
           divb = 0.
#ifdef FL_NON_PERMANENT_GUARDCELLS
           solnData(MAGX_VAR,i,j,k) = magx
           solnData(MAGY_VAR,i,j,k) = magy
           solnData(MAGZ_VAR,i,j,k) = magz
           solnData(MAGP_VAR,i,j,k) = magp
           solnData(DIVB_VAR,i,j,k) = divb
#else
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k

           call Grid_putPointData(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, magx)
           call Grid_putPointData(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, magy)
           call Grid_putPointData(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, magz)
           call Grid_putPointData(blockId, CENTER, MAGP_VAR, EXTERIOR, axis, magp)
           call Grid_putPointData(blockId, CENTER, DIVB_VAR, EXTERIOR, axis, divb)
#endif
#endif
        enddo
     enddo
  enddo


#if NFACE_VARS > 0
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           magx = 0.5*(facexData(MAG_FACE_VAR,i,j,k)+facexData(MAG_FACE_VAR,i+1,j,k))
           magy = 0.5*(faceyData(MAG_FACE_VAR,i,j,k)+faceyData(MAG_FACE_VAR,i,j+1,k))
           if (NDIM == 3) then
              magz = 0.5*(facezData(MAG_FACE_VAR,i,j,k)+facezData(MAG_FACE_VAR,i,j,k+1))
           endif

#if NDIM == 1
           divb = 0.
#elif NDIM >= 2
           divb = &
                     (facexData(MAG_FACE_VAR,i+1,j,  k  ) - facexData(MAG_FACE_VAR,i,j,k))/dx &
                   + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/dy
#if NDIM == 3
           divb = divb + (facezData(MAG_FACE_VAR,i,  j,  k+1) - facezData(MAG_FACE_VAR,i,j,k))/dz
#endif
#endif

           ! Update the magnetic pressure
           if (unitSystem .eq. 'cgs') then
               magp = .5/(4.*PI)*dot_product((/ magx, magy, magz /),&
                                             (/ magx, magy, magz /))
           else
               magp = .5*dot_product((/ magx, magy, magz /),&
                                     (/ magx, magy, magz /))
           endif

#ifdef FL_NON_PERMANENT_GUARDCELLS
           solnData(MAGX_VAR,i,j,k) = magx
           solnData(MAGY_VAR,i,j,k) = magy
           solnData(MAGZ_VAR,i,j,k) = magz
           solnData(MAGP_VAR,i,j,k) = magp
           solnData(DIVB_VAR,i,j,k) = divb
#else
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k

           call Grid_putPointData(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, magx)
           call Grid_putPointData(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, magy)
           call Grid_putPointData(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, magz)
           call Grid_putPointData(blockId, CENTER, MAGP_VAR, EXTERIOR, axis, magp)
           call Grid_putPointData(blockId, CENTER, DIVB_VAR, EXTERIOR, axis, divb)
#endif
        enddo
     enddo
  enddo
#endif !NFACE_VARS

#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  deallocate(xCoord)
  deallocate(xCoordL)
  deallocate(xCoordR)

  deallocate(yCoord)
  deallocate(yCoordL)
  deallocate(yCoordR)

  deallocate(zCoord)
  deallocate(zCoordL)
  deallocate(zCoordR)

#ifndef FIXEDBLOCKSIZE
  deallocate(Az)
  deallocate(Ax)
  deallocate(Ay)
#endif

  return
end subroutine Simulation_initBlock


!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(nn) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, nn, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: nn
  integer, intent(OUT):: i
  real, intent(IN)    :: x(nn), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(nn)) then

     i = nn+1

  else

     il = 1
     ir = nn
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif

  return
end subroutine sim_find
