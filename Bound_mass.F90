subroutine Bound_mass(blockCount, blockList)
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize
    use Gravity_data, ONLY: grv_thresh, grv_bound, grv_boundvec, grv_ptmass, grv_obvec, grv_ptvec, grv_peakvec
    use gr_mpoleData, ONLY: X_centerofmass, Y_centerofmass, Z_centerofmass
    use gr_isoMpoleData, ONLY: Xcm, Ycm, Zcm
    use Driver_data, ONLY: dr_simTime, dr_globalMe
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Simulation_data, ONLY: sim_tRelax, sim_fluffDampCutoff, sim_softenRadius
    implicit none 
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ, lb, istat, i, j, k, ierr, it
    double precision,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    logical :: gcell = .true.
    double precision, dimension(:,:,:,:),pointer :: solnData
    double precision, dimension(7) :: lsum, gsum
    double precision, dimension(3) :: bndvel, ptvel, bndpos, relvel
    double precision :: dvol, delm, xx, yy, zz, bndv2, ptv2, ptrad
    double precision :: G, tinitial, dx, delxinv

    call PhysicalConstants_get("Newton", G)
    call RuntimeParameters_get('tinitial',tinitial)

    do it = 1, 1 !Converges quite quickly
        lsum = 0.d0
        gsum = 0.d0

        do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
            sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
            allocate(xCoord(sizeX),stat=istat)
            sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
            allocate(yCoord(sizeY),stat=istat)
            sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
            allocate(zCoord(sizeZ),stat=istat)

            if (NDIM == 3) call Grid_getCellCoords&
                                (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
            if (NDIM >= 2) call Grid_getCellCoords&
                                (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
            call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
            call Grid_getBlkPtr(blockList(lb),solnData)
            dx = xCoord(2) - xCoord(1)
            dvol = dx**3.d0
            delxinv = 0.5e0 / dx
            do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                zz = zCoord(k)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    yy = yCoord(j)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

                        if (solnData(DENS_VAR,i,j,k) .lt. sim_fluffDampCutoff) cycle !Don't include any mass below the damping cutoff

                        !if (it .eq. 1) then
                            bndvel = solnData(VELX_VAR:VELZ_VAR,i,j,k) - grv_peakvec(4:6)
                        !else
                        !    bndvel = solnData(VELX_VAR:VELZ_VAR,i,j,k) - relvel
                        !endif
                        bndv2 = 0.5d0*sum(bndvel**2.d0)
                        if (-solnData(GPOT_VAR,i,j,k) .ge. bndv2) then
                            xx = xCoord(i)

                            !if (dr_simTime .gt. sim_tRelax) then
                            !    bndpos = (/ xx, yy, zz /) - grv_peakvec(1:3) - grv_ptvec(1:3) + grv_obvec(1:3)
                            !    ptrad = dsqrt(sum(bndpos**2.d0))
                            !    ptvel = solnData(VELX_VAR:VELZ_VAR,i,j,k) + grv_obvec(4:6) - grv_ptvec(4:6)
                            !    ptv2 = sum(ptvel**2.d0)
                            !endif

                            !if (dr_simTime .le. sim_tRelax .or. &
                            !    ptv2 - G*grv_ptmass/ptrad .gt. bndv2 - solnData(GPOT_VAR,i,j,k)) then
                            delm = solnData(DENS_VAR,i,j,k)*dvol
                            lsum(1)  = lsum(1)  + delm
                            lsum(2)  = lsum(2)  + xx*delm
                            lsum(3)  = lsum(3)  + yy*delm
                            lsum(4)  = lsum(4)  + zz*delm
                            lsum(5)  = lsum(5)  + solnData(VELX_VAR,i,j,k)*delm
                            lsum(6)  = lsum(6)  + solnData(VELY_VAR,i,j,k)*delm
                            lsum(7)  = lsum(7)  + solnData(VELZ_VAR,i,j,k)*delm
                        endif
                    enddo
                enddo
            enddo

            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo

        call MPI_ALLREDUCE(lsum, gsum, 7, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
        grv_bound = gsum(1)
        grv_boundvec = gsum(2:7) / gsum(1)
    enddo

end subroutine Bound_mass
