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
    use Gravity_data, ONLY: grv_factor, grv_thresh, grv_ptvec, grv_obvec, grv_optvec, grv_oobvec, grv_mode, &
        grv_hptvec, grv_hobvec, grv_exactvec, grv_oexactvec, grv_obaccel, &
        grv_ptaccel, grv_hptaccel, grv_optaccel, grv_oobaccel
    use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
        Grid_releaseBlkPtr, Grid_getCellCoords, Grid_getBlkIndexLimits, Grid_getMinCellSize
    use Simulation_data, ONLY: sim_tRelax, sim_softenRadius, sim_fluffDampCutoff
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use tree, ONLY : lrefine_max
    use gr_mpoleData, ONLY : twelfth
    use Grid_data, ONLY : gr_meshMe
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
    integer         :: potVar, denVar

#ifdef FIXEDBLOCKSIZE
    double precision,dimension(GRID_KHI_GC) :: zCenter
    double precision,dimension(GRID_JHI_GC) :: yCenter
    double precision,dimension(GRID_IHI_GC) :: xCenter
#else
    double precision,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
    double precision :: dr32, tmpdr32, min_cell

    integer :: sizeX,sizeY,sizez

    integer :: j,k
    logical :: gcell = .true.
    double precision :: period
    double precision :: tinitial, ptx, pty, ptz, deld
    double precision, dimension(3) :: obaccel, ptaccel
    !==================================================

    call Grid_getBlkPhysicalSize(blockID, blockSize)

    call Grid_getBlkPtr(blockID, solnVec)

    !! IF a variable index is explicitly specified, assume that as the potential
    !! otherwise use the default current potential GPOT_VAR  
    if(present(potentialIndex)) then
        potVar=potentialIndex
    else
        potVar=GPOT_VAR
    endif

    if (potVar .eq. GPOT_VAR) then
        denVar = DENS_VAR
    else
        denVar = ODEN_VAR
    endif

    if (grv_mode .eq. 1) then
        ptx = grv_oexactvec(1) + (grv_optvec(1) - grv_oobvec(1))
        pty = grv_oexactvec(2) + (grv_optvec(2) - grv_oobvec(2))
        ptz = grv_oexactvec(3) + (grv_optvec(3) - grv_oobvec(3))
        obaccel = grv_oobaccel
        ptaccel = grv_optaccel
    elseif (grv_mode .eq. 2) then
        ptx = (grv_oexactvec(1) + grv_exactvec(1))/2.d0 + (grv_hptvec(1) - grv_hobvec(1))
        pty = (grv_oexactvec(2) + grv_exactvec(2))/2.d0 + (grv_hptvec(2) - grv_hobvec(2))
        ptz = (grv_oexactvec(3) + grv_exactvec(3))/2.d0 + (grv_hptvec(3) - grv_hobvec(3))
        ptaccel = grv_hptaccel
    elseif (grv_mode .eq. 3) then
        ptx = grv_exactvec(1) + (grv_ptvec(1) - grv_obvec(1))
        pty = grv_exactvec(2) + (grv_ptvec(2) - grv_obvec(2))
        ptz = grv_exactvec(3) + (grv_ptvec(3) - grv_obvec(3))
        obaccel = grv_obaccel
        ptaccel = grv_ptaccel
    endif

    if (grv_mode .ne. 2) then
        iimin   = 1
        iimax   = numCells

        !Get row of potential values and compute inverse of zone spacing  
        if (sweepDir == SWEEP_X) then                     ! x-direction
            delxinv = dble(NXB) / blockSize(IAXIS)
            
            gpot(:) = solnVec(potVar,:,pos(1),pos(2))
            dens(:) = solnVec(denVar,:,pos(1),pos(2))

            grav = -obaccel(1)
        elseif (sweepDir == SWEEP_Y) then                 ! y-direction
            delxinv = dble(NYB) / blockSize(JAXIS)
            
            gpot(:) = solnVec(potVar,pos(1),:,pos(2))
            dens(:) = solnVec(denVar,pos(1),:,pos(2))

            grav = -obaccel(2)
        else                                            ! z-direction
            delxinv = dble(NZB) / blockSize(KAXIS)
            
            gpot(:) = solnVec(potVar,pos(1),pos(2),:)
            dens(:) = solnVec(denVar,pos(1),pos(2),:)

            grav = -obaccel(3)
        endif

        !-------------------------------------------------------------------------------

        !               Compute gravitational acceleration


        !**************** first-order differences
        !                 preserves conservation

        !delxinv = 0.5e0 * delxinv
        !
        !do ii = iimin+1, iimax-1
        !   grav(ii) = delxinv * (gpot(ii-1) - gpot(ii+1))
        !enddo
        !
        !grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
        !grav(iimax) = grav(iimax-1)

        !Bryan, Norman, Stone 1995
        delxinv = 0.5e0 * delxinv
        do ii = iimin+1, iimax-1
            if (dens(ii) .lt. sim_fluffDampCutoff) cycle
            if ((dens(ii+1) - dens(ii))*(dens(ii) - dens(ii-1)) .lt. 0.d0) then
                deld = sign(min(dabs(dens(ii-1) - dens(ii+1)), 2.d0*dabs(dens(ii) - dens(ii-1)), 2.d0*dabs(dens(ii) - dens(ii+1))), &
                    dens(ii+1) - dens(ii-1))
            else
                deld = 0.d0
            endif
            grav(ii) = grav(ii) + delxinv * (gpot(ii-1) - gpot(ii+1) + deld/dens(ii)*twelfth*(gpot(ii-1) - 2.d0*gpot(ii) + gpot(ii+1)))
        enddo
    endif

    !Now include the point mass
    call RuntimeParameters_get('tinitial',tinitial)
    ptgrav = 0.d0
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

      if (sweepDir .eq. SWEEP_X) then                       ! x-component

         tmpdr32 = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 
         ptgrav = -ptaccel(1)

         do ii = 1, numCells
            if (dens(ii) .lt. sim_fluffDampCutoff) cycle
            dr32 = dsqrt(xCenter(ii)*xCenter(ii) + tmpdr32)
            if (dr32 .lt. sim_softenRadius) then
                dr32 = sim_softenRadius*sim_softenRadius*dr32
            else
                dr32 = dr32*dr32*dr32
            endif
            ptgrav(ii) = ptgrav(ii) + grv_factor*xCenter(ii)/dr32
         enddo


      else if (sweepDir .eq. SWEEP_Y) then          ! y-component

         tmpdr32 = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 
         ptgrav = -ptaccel(2)

         do ii = 1, numCells
            if (dens(ii) .lt. sim_fluffDampCutoff) cycle
            dr32 = dsqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
            if (dr32 .lt. sim_softenRadius) then
                dr32 = sim_softenRadius*sim_softenRadius*dr32
            else
                dr32 = dr32*dr32*dr32
            endif
            ptgrav(ii) = ptgrav(ii) + grv_factor*yCenter(ii)/dr32
         enddo

      else if (sweepDir .eq. SWEEP_Z) then          ! z-component

         tmpdr32 = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 
         ptgrav = -ptaccel(3)

         do ii = 1, numCells
            if (dens(ii) .lt. sim_fluffDampCutoff) cycle
            dr32 = dsqrt(zCenter(ii)*zCenter(ii) + tmpdr32)           
            if (dr32 .lt. sim_softenRadius) then
                dr32 = sim_softenRadius*sim_softenRadius*dr32
            else
                dr32 = dr32*dr32*dr32
            endif
            ptgrav(ii) = ptgrav(ii) + grv_factor*zCenter(ii)/dr32
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


