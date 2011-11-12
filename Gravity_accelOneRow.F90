!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2),intent(in):: pos, 
!!                      integer, intent(in) :: sweepDir, 
!!                      integer, intent(in) :: blockID, 
!!                      integer, intent(in) :: numCells, 
!!                      double precision(numCells),intent(out) :: grav, 
!!                      integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  test against sweep_x,
!!                                 sweep_y, and sweep_z
!!  blockID   -     The local identifier of the block to work on
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!
!!
!!***


subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, ptgrav, potentialIndex)

    use Driver_data, ONLY: dr_simTime
    use Gravity_data, ONLY: grv_factor, grv_thresh, grv_ptvec, grv_obvec, grv_optvec, grv_oobvec
    use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
        Grid_releaseBlkPtr, Grid_getCellCoords, Grid_getBlkIndexLimits, Grid_getMinCellSize
    use Simulation_data, ONLY: sim_tRelax, sim_sinkRadius, sim_sinkRef
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use gr_mpoleData, ONLY: Xcm, Ycm, Zcm, oXcm, oYcm, oZcm
    use tree, ONLY : lrefine_max
    implicit none

#include "Flash.h"
#include "constants.h"

    integer, dimension(2), intent(in) :: pos
    integer, intent(in)               :: sweepDir, blockID,  numCells
    double precision, intent(inout)               :: grav(numCells)
    double precision, intent(inout)               :: ptgrav(numCells)
    integer, intent(IN),optional      :: potentialIndex
    double precision            :: blockSize(MDIM)
    double precision, POINTER, DIMENSION(:,:,:,:) :: solnVec

    integer         :: ii, iimin, iimax, lb
    double precision            :: gpot(numCells), delxinv
    double precision            :: dens(numCells), xgrid(numCells)
    double precision, parameter :: onesixth = 1.e0/6.e0
    integer         :: potVar

#ifdef FIXEDBLOCKSIZE
    double precision,dimension(GRID_KHI_GC) :: zCenter
    double precision,dimension(GRID_JHI_GC) :: yCenter
    double precision,dimension(GRID_IHI_GC) :: xCenter
#else
    double precision,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
    double precision :: dr32, tmpdr32, drs, min_cell, inner_rad

    integer :: sizeX,sizeY,sizez

    integer :: j,k,mode
    logical :: gcell = .true.
    double precision :: period
    double precision :: accelx, accely, accelz, tinitial, ptx, pty, ptz, deld
    !==================================================



    call Grid_getBlkPhysicalSize(blockID, blockSize)


    call Grid_getBlkPtr(blockID, solnVec)

    !! IF a variable index is explicitly specified, assume that as the potential
    !! otherwise use the default current potential GPOT_VAR  
    if(present(potentialIndex)) then
        potVar=potentialIndex
        if (potVar .eq. GPOT_VAR) then
            ptx = Xcm + (grv_ptvec(1) - grv_obvec(1))
            pty = Ycm + (grv_ptvec(2) - grv_obvec(2))
            ptz = Zcm + (grv_ptvec(3) - grv_obvec(3))
            mode = 1
        elseif (potVar .eq. GPOL_VAR) then
            ptx = oXcm + (grv_optvec(1) - grv_oobvec(1))
            pty = oYcm + (grv_optvec(2) - grv_oobvec(2))
            ptz = oZcm + (grv_optvec(3) - grv_oobvec(3))
            mode = 2
        endif
    else
        potVar=GPOT_VAR
        ptx = Xcm + (grv_ptvec(1) - grv_obvec(1))
        pty = Ycm + (grv_ptvec(2) - grv_obvec(2))
        ptz = Zcm + (grv_ptvec(3) - grv_obvec(3))
        mode = 1
    endif

    iimin   = 1
    iimax   = numCells


    !Get row of potential values and compute inverse of zone spacing  
    if (sweepDir == SWEEP_X) then                     ! x-direction

     delxinv = dble(NXB) / blockSize(IAXIS)
     
     gpot(:) = solnVec(potVar,:,pos(1),pos(2))
     dens(:) = solnVec(DENS_VAR,:,pos(1),pos(2))
     
    elseif (sweepDir == SWEEP_Y) then                 ! y-direction
     
     delxinv = dble(NYB) / blockSize(JAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),:,pos(2))
     dens(:) = solnVec(DENS_VAR,pos(1),:,pos(2))
     
    else                                            ! z-direction
     
     delxinv = dble(NZB) / blockSize(KAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),pos(2),:)
     dens(:) = solnVec(DENS_VAR,pos(1),pos(2),:)
     
    endif

    !-------------------------------------------------------------------------------

    !               Compute gravitational acceleration


    !**************** first-order differences
    !                 preserves conservation

    !delxinv = delxinv / 12.d0

    !do ii = iimin+2, iimax-2
    !    grav(ii) = delxinv * (-gpot(ii-2) + 8.d0*gpot(ii-1) - 8.d0*gpot(ii+1) + gpot(ii+2))
    !enddo

    !grav(iimin+1) = delxinv * (gpot(iimin) - gpot(iimin+2))
    !grav(iimax-1) = delxinv * (gpot(iimax-2) - gpot(iimax))

    !grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
    !grav(iimax) = grav(iimax-1)

    !Bryan, Norman, Stone 1995
    delxinv = 0.5e0 * delxinv
    do ii = iimin+1, iimax-1
        if ((dens(ii+1) - dens(ii))*(dens(ii) - dens(ii-1)) .lt. 0.d0) then
            deld = sign(min(abs(dens(ii-1) - dens(ii+1)), 2.d0*abs(dens(ii) - dens(ii-1)), 2.d0*abs(dens(ii) - dens(ii+1))), &
                dens(ii+1) - dens(ii-1))
        else
            deld = 0.d0
        endif
        grav(ii) = delxinv * (gpot(ii-1) - gpot(ii+1) + deld/dens(ii)/12.d0*(gpot(ii-1) - 2.d0*gpot(ii) + gpot(ii+1)))
    enddo

    !delxinv = 0.5e0 * delxinv
    !do ii = iimin+1, iimax-1
    !    grav(ii) = delxinv * (gpot(ii-1) - gpot(ii+1))
    !enddo
    !grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
    !grav(iimax) = grav(iimax-1)

    !xgrid(iimin) = 0
    !do ii = iimin+1, iimax
    !    xgrid(ii) = xgrid(ii-1) + 1.d0/delxinv
    !enddo
    !call fss004(xgrid, numCells, gpot, grav)
    !grav = -grav

    !Now include the point mass
    call RuntimeParameters_get('tinitial',tinitial)
    if (dr_simTime .ge. tinitial + sim_tRelax) then
      j=pos(1)
      k=pos(2)
#ifndef FIXEDBLOCKSIZE
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      sizeX=blkLimitsGC(HIGH,IAXIS)
      sizeY=blkLimitsGC(HIGH,JAXIS)
      sizeZ=blkLimitsGC(HIGH,KAXIS)
      allocate(xCenter(sizeX))
      allocate(yCenter(sizeY))
      allocate(zCenter(sizeZ))
#else
      sizeX=GRID_IHI_GC
      sizeY=GRID_JHI_GC
      sizeZ=GRID_KHI_GC
#endif
      if (mode .le. 2) call com_accel(ptx, pty, ptz, mode, accelx, accely, accelz)
      
      zCenter = 0.
      yCenter = 0.
      if (NDIM == 3) then 
         call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
         zCenter = zCenter - ptz
      endif
      if (NDIM >= 2) then
         call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
         yCenter = yCenter - pty
      endif
      call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
      xCenter = xCenter - ptx
      
      call Grid_getMinCellSize(min_cell)
      inner_rad = sim_sinkRadius*2.**(lrefine_max-sim_sinkRef)*min_cell
      if (sweepDir .eq. SWEEP_X) then                       ! x-component

         tmpdr32 = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 

         do ii = 1, numCells

            dr32 = sqrt(xCenter(ii)*xCenter(ii) + tmpdr32)
            if (dr32 .le. inner_rad) then
                drs = inner_rad
            else
                drs = dr32
            endif
            dr32 = drs*drs*dr32

            if (mode .le. 2) ptgrav(ii) = grv_factor*xCenter(ii)/dr32 - accelx
         enddo


      else if (sweepDir .eq. SWEEP_Y) then          ! y-component

         Tmpdr32 = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 

         do ii = 1, numCells
            
            dr32 = sqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
            if (dr32 .le. inner_rad) then
                drs = inner_rad
            else
                drs = dr32
            endif
            dr32 = drs*drs*dr32

            if (mode .le. 2) ptgrav(ii) = grv_factor*yCenter(ii)/dr32 - accely
         enddo

      else if (sweepDir .eq. SWEEP_Z) then          ! z-component

         tmpdr32 = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 

         do ii = 1, numCells
            
            dr32 = sqrt(zCenter(ii)*zCenter(ii) + tmpdr32)           
            if (dr32 .le. inner_rad) then
                drs = inner_rad
            else
                drs = dr32
            endif
            dr32 = drs*drs*dr32
            
            if (mode .le. 2) ptgrav(ii) = grv_factor*zCenter(ii)/dr32 - accelz
         enddo

      endif
    endif

    !==============================================================================
#ifndef FIXEDBLOCKSIZE
    deallocate(xCenter)
    deallocate(yCenter)
    deallocate(zCenter)
#endif

    call Grid_releaseBlkPtr(blockID, solnVec)
    return
   
end subroutine Gravity_accelOneRow


