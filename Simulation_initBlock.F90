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
    Grid_getCellCoords, Grid_putPointData
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Multispecies_interface, ONLY:  Multispecies_getSumFrac, Multispecies_getSumInv, Multispecies_getAvg
  
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
  real     ::  distInv, xDist, yDist, zDist
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek, t, mp, kb
  real     ::  dist, gam
  integer  ::  istat

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
#ifdef FL_NON_PERMANENT_GUARDCELLS
  real, dimension(:,:,:,:),pointer :: solnData
#endif

  logical :: gcell = .true.

#ifdef LOADPROFILE
  real, dimension(sim_tableCols) :: sumVars
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: xn = 0.0
#else
  real, dimension(2) :: sumVars
#endif

  call PhysicalConstants_get("proton mass", mp)
  call PhysicalConstants_get("Boltzmann", kb)
  
  call Multispecies_getAvg(GAMMA, gam)


  !do i = 1, np
  !   print *, obj_radius(i), obj_rhop(i), eint(i) 
  !enddo
  
  ! get the coordinate information for the current block from the database

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
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

           do kk = 1, sim_nSubZones
              zz    = zCoord(k) + dzz*((2*kk - 1)*sim_inSubInv - 0.5)
              zDist = zz - sim_zCenter
              
              do jj = 1, sim_nSubZones
                 yy    = yCoord(j) + dyy*((2*jj - 1)*sim_inSubInv - 0.5)
                 yDist = yy - sim_yCenter
                 
                 do ii = 1, sim_nSubZones
                    xx    = xCoord(i) + dxx*((2*ii - 1)*sim_inSubInv - 0.5)
                    xDist = xx - sim_xCenter
                    
                    dist    = dsqrt( xDist**2 + yDist**2 + zDist**2 )
                    distInv = 1. / max( dist, 1.E-10 )
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
#endif
                 enddo
              enddo
           enddo
           
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
           rho = max (sumVars(1)*sim_inszd, sim_rhoAmbient)
           p   = max (sumVars(2)*sim_inszd, sim_pAmbient)
           vx  = 0.0d0
           vy  = 0.0d0
           vz  = 0.0d0
           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           !
           !  assume gamma-law equation of state
           !
           t   = p/(rho/mp/obj_mu*kb)
           e   = p/(gam-1.)
           e   = e/rho + ek
           e   = max (e, sim_smallP)
#endif
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
#ifdef LOADPROFILE
           solnData(DENS_VAR,i,j,k)=sumVars(RHO_PROF)
           solnData(TEMP_VAR,i,j,k)=sumVars(TEMP_PROF)
           do put=1,NSPECIES
              if (xn(SPECIES_BEGIN+put-1) == 0.0) xn(SPECIES_BEGIN+put-1) = sim_smallX
              solnData(SPECIES_BEGIN+put-1,i,j,k)=xn(SPECIES_BEGIN+put-1)
           enddo
#else
           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
           solnData(TEMP_VAR,i,k,k)=t
           solnData(GAME_VAR,i,j,k)=gam
           solnData(GAMC_VAR,i,j,k)=gam
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
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, gam)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, gam)
#endif
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
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
