!!****f* source/physics/Gravity/Gravity_computeDt
!!
!! NAME
!!
!!  Gravity_computeDt
!!  
!! SYNOPSIS
!!
!!  Gravity_computeDt(integer(IN)        :: blockID,
!!                    
!!                    real (OUT)         :: dt_grav,
!!                    integer(:)(INOUT)  :: dt_minloc(5))
!!
!! DESCRIPTION
!!
!!  Compute the timestep limiter due to the gravitational solver.
!!
!! ARGUMENTS
!!
!!  MyPE:          local processor number
!!  dt_grav:       Will Return the limiting timestep. Should be
!!                 set to a large value (1.D99) on input.
!!  dt_minloc(5):  An array to receive information about which
!!                 processor, block, and zone was responsible
!!                 for setting the limiting timestep.  The order
!!                 is i, j, k, b, p, where (i,j,k) = zone
!!                 indices, b = local block ID, and p = PE #.
!!                 This routine should only modify these values
!!                 if it changes dt_grav.
!!  blockID:       The local ID of the block to compute the
!!                 limiter on.
!!
!!***

subroutine Gravity_computeDt (blockID, dt_grav, dt_minloc)

!==============================================================================
    use Grid_interface, ONLY : Grid_getMinCellSize, Grid_getBlkPtr, &
        Grid_getBlkIndexLimits, Grid_releaseBlkPtr
    use Gravity_data, ONLY : grv_obvec, grv_ptvec, grv_exactvec, grv_ptaccel, &
        grv_ptmass, grv_cfl
    use Hydro_data, ONLY: hy_cfl
    use Driver_data, ONLY: dr_simTime
    use Simulation_data, ONLY: sim_tRelax, sim_fluffDampCutoff, sim_softenRadius
    use gr_mpoleData, ONLY: twelfth
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  
    implicit none

#include "Flash.h"
#include "constants.h"
  
    integer, intent(IN)    ::  blockID
    integer, intent(INOUT) ::  dt_minloc(5)
    real,intent(OUT)       ::  dt_grav

    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    real, DIMENSION(GRID_ILO_GC:GRID_IHI_GC,                   &
                    GRID_JLO_GC:GRID_JHI_GC,                   &
                    GRID_KLO_GC:GRID_KHI_GC,MDIM) :: grav, ptgrav
    real, dimension(:,:,:,:),pointer :: solnData
    double precision,dimension(MDIM) :: blockSize
    real :: mcs, delxinv, deld, newton, gm, dr32
    integer :: i, j, k
    double precision,dimension(GRID_KHI_GC) :: zCenter
    double precision,dimension(GRID_JHI_GC) :: yCenter
    double precision,dimension(GRID_IHI_GC) :: xCenter
    
    call Grid_getBlkPhysicalSize(blockID, blockSize)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    delxinv = dble(NXB) / blockSize(IAXIS)
    delxinv = 0.5e0 * delxinv

    grav = 0.d0
    ptgrav = 0.d0

    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)        
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)                         
            do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                if (solnData(DENS_VAR,i,j,k) .lt. sim_fluffDampCutoff) cycle

                if ((solnData(DENS_VAR,i+1,j,k) - solnData(DENS_VAR,i,j,k))*&
                    (solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i-1,j,k)) .lt. 0.d0) then
                    deld = sign(min(dabs(solnData(DENS_VAR,i-1,j,k) - solnData(DENS_VAR,i+1,j,k)), &
                        2.d0*dabs(solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i-1,j,k)), &
                        2.d0*dabs(solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i+1,j,k))), &
                        solnData(DENS_VAR,i+1,j,k) - solnData(DENS_VAR,i-1,j,k))
                else
                    deld = 0.d0
                endif
                grav(i,j,k,1) = dabs(delxinv * &
                    (solnData(GPOT_VAR,i-1,j,k) - solnData(GPOT_VAR,i+1,j,k) + &
                    deld/solnData(DENS_VAR,i,j,k)*twelfth*(solnData(GPOT_VAR,i-1,j,k) - &
                    2.d0*solnData(GPOT_VAR,i,j,k) + solnData(GPOT_VAR,i+1,j,k))))

                if ((solnData(DENS_VAR,i,j+1,k) - solnData(DENS_VAR,i,j,k))*&
                    (solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i,j-1,k)) .lt. 0.d0) then
                    deld = sign(min(dabs(solnData(DENS_VAR,i,j-1,k) - solnData(DENS_VAR,i,j+1,k)), &
                        2.d0*dabs(solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i,j-1,k)), &
                        2.d0*dabs(solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i,j+1,k))), &
                        solnData(DENS_VAR,i,j+1,k) - solnData(DENS_VAR,i,j-1,k))
                else
                    deld = 0.d0
                endif
                grav(i,j,k,2) = dabs(delxinv * &
                    (solnData(GPOT_VAR,i,j-1,k) - solnData(GPOT_VAR,i,j+1,k) + &
                    deld/solnData(DENS_VAR,i,j,k)*twelfth*(solnData(GPOT_VAR,i,j-1,k) - &
                    2.d0*solnData(GPOT_VAR,i,j,k) + solnData(GPOT_VAR,i,j+1,k))))

                if ((solnData(DENS_VAR,i,j,k+1) - solnData(DENS_VAR,i,j,k))*&
                    (solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i,j,k-1)) .lt. 0.d0) then
                    deld = sign(min(dabs(solnData(DENS_VAR,i,j,k-1) - solnData(DENS_VAR,i,j,k+1)), &
                        2.d0*dabs(solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i,j,k-1)), &
                        2.d0*dabs(solnData(DENS_VAR,i,j,k) - solnData(DENS_VAR,i,j,k+1))), &
                        solnData(DENS_VAR,i,j,k+1) - solnData(DENS_VAR,i,j,k-1))
                else
                    deld = 0.d0
                endif
                grav(i,j,k,3) = dabs(delxinv * &
                    (solnData(GPOT_VAR,i,j,k-1) - solnData(GPOT_VAR,i,j,k+1) + &
                    deld/solnData(DENS_VAR,i,j,k)*twelfth*(solnData(GPOT_VAR,i,j,k-1) - &
                    2.d0*solnData(GPOT_VAR,i,j,k) + solnData(GPOT_VAR,i,j,k+1))))
            enddo
        enddo
    enddo

    if (dr_simTime .gt. sim_tRelax) then
        call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xCenter, GRID_IHI_GC)
        xCenter = xCenter - (grv_exactvec(1) + (grv_ptvec(1) - grv_obvec(1)))
        call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yCenter, GRID_JHI_GC)
        yCenter = yCenter - (grv_exactvec(2) + (grv_ptvec(2) - grv_obvec(2)))
        call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zCenter, GRID_KHI_GC)
        zCenter = zCenter - (grv_exactvec(3) + (grv_ptvec(3) - grv_obvec(3)))

        call PhysicalConstants_get("Newton", newton)
        gm = -newton*grv_ptmass

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)        
            do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)                         
                do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                    if (solnData(DENS_VAR,i,j,k) .lt. sim_fluffDampCutoff) cycle
                    dr32 = dsqrt(xCenter(i)*xCenter(i) + yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k))
                    if (dr32 .lt. sim_softenRadius) then
                        dr32 = sim_softenRadius*sim_softenRadius*dr32
                    else
                        dr32 = dr32*dr32*dr32
                    endif
                    ptgrav(i,j,k,:) = dabs(gm*(/ xCenter(i), yCenter(j), zCenter(k) /)/dr32 - grv_ptaccel)
                enddo
            enddo
        enddo
    endif

    call Grid_releaseBlkPtr(blockID, solnData)
  
    call Grid_getMinCellSize(mcs)
    if (dr_simTime .gt. sim_tRelax) then
        dt_grav = min(huge(dt_grav), minval(grv_cfl*dsqrt(mcs/grav)), minval(grv_cfl*dsqrt(mcs/ptgrav)))!, minval(hy_cfl*mcs/abs(grv_obvec(4:6) - grv_ptvec(4:6))))
    else
        dt_grav = min(huge(dt_grav), minval(grv_cfl*dsqrt(mcs/grav)))
    endif
    !dt_grav = huge(dt_grav)
    
    return

end subroutine Gravity_computeDt
