!!****if* source/Simulation/SimulationMain/MultiTidalBiPoly/Simulation_initBlock
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
!!  a specified block.  This version sets up MultiTidalBiPoly.
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
  
  integer,intent(IN) ::  blockId
  integer,intent(IN) ::  myPE
  
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk, put
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumP, sumMu, sumGam
  double precision,dimension(SPECIES_BEGIN:SPECIES_END) ::	sumXn
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac, efrac
  real     ::  vx, vy, vz, p, rho, e, ek, t, mp, kb
  real     ::  dist, gam, mu
  logical  ::  validGeom
  integer  ::  istat

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData

  logical :: gcell = .true.

  call PhysicalConstants_get("proton mass", mp)
  call PhysicalConstants_get("Boltzmann", kb)
  
!  call Multispecies_getAvg(GAMMA, gam)


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
           
           sumRho = 0.
           sumP   = 0.
		   sumGam = 0.
		   sumMu  = 0.
		   sumXn(H1_SPEC)   = 0.
		   sumXn(HE4_SPEC)  = 0.
           sumXn(CORE_SPEC) = 0.
           !
           !       Break the cell into sim_nSubZones^NDIM sub-zones, and look up the
           !       appropriate quantities along the 1d profile for that subzone.  
           !
           !       Have the final values for the zone be equal to the average of
           !       the subzone values.
           ! 

           do kk = 0, (sim_nSubZones-1)*K3D
              zz    = zCoord(k) + (kk*sim_inSubzm1-.5)*dzz 
              zDist = (zz - (sim_zMax - sim_zMin) / 2.d0) * K3D
              
              do jj = 0, (sim_nSubZones-1)*K2D
                 yy    = yCoord(j) + (jj*sim_inSubzm1-.5)*dyy
                 yDist = (yy - (sim_yMax - sim_yMin) / 2.d0) * K2D
                 
                 do ii = 0, (sim_nSubZones-1)
                    xx    = xCoord(i) + (ii*sim_inSubzm1-.5)*dxx
                    xDist = xx - (sim_xMax - sim_xMin) / 2.d0
                    
                    dist    = dsqrt( xDist**2 + yDist**2 + zDist**2 )
                    distInv = 1. / max( dist, 1.E-10 )
                    call sim_find (obj_radius, obj_ipos, dist, jLo)
                    !
                    !  a point at `dist' is frac-way between jLo and jHi.   We do a
                    !  linear interpolation of the quantities at jLo and jHi and sum those.
                    ! 
                    if (jLo .eq. 0) then
                       jLo = 1
                       jHi = 1
                       frac = 0.
                    else if (jLo .eq. obj_ipos) then
                       jLo = obj_ipos
                       jHi = obj_ipos
                       frac = 0.
                       efrac = 1.
                    else
                       jHi = jLo + 1
                       frac = (dist - obj_radius(jLo)) / & 
                            (obj_radius(jHi)-obj_radius(jLo))
						if (jHi .le. obj_ipoi) then
							efrac = 0.
						else if (jLo .gt. obj_ipoi) then
							efrac = 1.
						else if ((jLo .eq. obj_ipoi) .and. (dist .gt. obj_radius(obj_ipoi))) then
							efrac = (obj_radius(jHi) - dist) / &
								(obj_radius(jHi)-obj_radius(jLo))
						else
							efrac = 0.
						endif
                    endif
                    ! 
                    !   Now total these quantities.   Note that  v is a radial velocity; 
                    !   we multiply by the tangents of the appropriate angles to get
                    !   the projections in the x, y and z directions.
                    !
                    sumP = sumP +  & 
                         obj_prss(jLo) + frac*(obj_prss(jHi) - obj_prss(jLo))
                    
                    sumRho = sumRho + & 
                         obj_rhop(jLo) + frac*(obj_rhop(jHi) - obj_rhop(jLo))

                    sumGam = sumGam + &
                         gammac + efrac*(gammae - gammac)
                    
                    sumMu  = sumMu + &
                         obj_muc + efrac*(obj_mue -obj_muc)

                    sumXn(H1_SPEC) = sumXn(H1_SPEC) + &
                         0.0 + efrac*(0.7 - 0.0)
                    
                    sumXn(HE4_SPEC) = sumXn(HE4_SPEC) + &
                         0.0 + efrac*(0.3 - 0.0)
                    
                    sumXn(CORE_SPEC) = sumXn(CORE_SPEC) + &
                         1.0 + efrac*(0.0 - 1.0)
                 enddo
              enddo
           enddo
           
           rho = max (sumRho * sim_inszd, sim_rhoAmbient)
           p   = max (sumP   * sim_inszd, sim_pAmbient)
           vx  = 0.0d0
           vy  = 0.0d0
           vz  = 0.0d0
           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           !
           !  assume gamma-law equation of state
           !
		   gam = sumGam * sim_inszd
		   mu  = sumMu  * sim_inszd
		   obj_xn(H1_SPEC)   = sumXn(H1_SPEC)   * sim_inszd 
		   obj_xn(HE4_SPEC)  = sumXn(HE4_SPEC)  * sim_inszd
		   obj_xn(CORE_SPEC) = sumXn(CORE_SPEC) * sim_inszd
!           t   = p/(rho/mp/obj_mu*kb)
           t   = p/(rho/mp/mu*kb)
           e   = p/(gam-1.)
           e   = e/rho + ek
           e   = max (e, sim_smallP)
           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
           solnData(TEMP_VAR,i,k,k)=t
           solnData(GAME_VAR,i,j,k)=gam
           solnData(GAMC_VAR,i,j,k)=gam
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz
#else
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, t)    
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, gam)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, gam)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
           do put=SPECIES_BEGIN,SPECIES_END
              if (obj_xn(put) == 0.0) obj_xn(put) = sim_smallX
              call Grid_putPointData(blockID,CENTER,put,&
                   EXTERIOR,axis,obj_xn(put))
           enddo
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
