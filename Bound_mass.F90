subroutine Bound_mass(blockCount, blockList)
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize, Grid_getMyPE
    use Gravity_data, ONLY: grv_thresh, grv_bound, grv_boundvec
    implicit none 
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    double precision :: dvol, delm, xx, yy, zz
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ, lb, istat, i, j, k, ierr
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    logical :: gcell = .true.
    real, dimension(:,:,:,:),pointer :: solnData
    double precision, dimension(7) :: lsum, gsum

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
        dvol = (xCoord(2) - xCoord(1))**3.d0
        do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
            zz = zCoord(k)
            do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                yy = yCoord(j)
                do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                    if (-solnData(GPOT_VAR,i,j,k) .ge. sum(solnData(VELX_VAR:VELZ_VAR,i,j,k)**2.d0)) then
                        xx = xCoord(i)
                        delm = solnData(DENS_VAR,i,j,k)*dvol
                        lsum(1) = lsum(1) + delm
                        lsum(2) = lsum(2) + xx*delm
                        lsum(3) = lsum(3) + yy*delm
                        lsum(4) = lsum(4) + zz*delm
                        lsum(5) = lsum(5) + solnData(VELX_VAR,i,j,k)*delm
                        lsum(6) = lsum(6) + solnData(VELY_VAR,i,j,k)*delm
                        lsum(7) = lsum(7) + solnData(VELZ_VAR,i,j,k)*delm
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
end subroutine Bound_mass
