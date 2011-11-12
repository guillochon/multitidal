subroutine Orbit_energy(blockCount, blockList)

    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize, Grid_getDeltas,&
        Grid_getBlkBoundBox
    use Gravity_data, ONLY: grv_thresh, grv_ener, grv_ptvec, grv_obvec, grv_ptmass, grv_tot_ener, &
        grv_totmass, grv_exactvec, grv_factor
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use gr_mpoleData, ONLY: X_centerofmass, Y_centerofmass, Z_centerofmass
    use gr_isoMpoleData, ONLY: Xcm, Ycm, Zcm
    use Driver_data, ONLY : dr_dtOld
    use Simulation_data, ONLY : sim_softenRadius, sim_fluffDampCutoff

    implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    double precision :: ener, dvol, newton, r, vel2, xx, yy, zz, tot_ener
    double precision :: x, y, z, min_cell
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ, lb, istat, i, j, k, ierr, imin, imax, jmin, jmax, kmin, kmax
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    logical :: gcell = .true.
    real, dimension(:,:,:,:),pointer :: solnData
    real,dimension(MDIM) :: delta
    real,dimension(LOW:HIGH,MDIM) :: bndBox
    double precision :: bvx, bvy, bvz, inner_rad

    call PhysicalConstants_get("Newton", newton)
    call Grid_getMinCellSize(min_cell)
    ener = 0.d0
    tot_ener = 0.d0
    
    inner_rad = sim_softenRadius

    do lb = 1, blockCount
        call Grid_getDeltas(blockList(lb),delta)
        call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkBoundBox(blockList(lb), bndBox)
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
        dvol = delta(IAXIS) * delta(JAXIS) * delta(KAXIS)
        kmax = blkLimits(HIGH,KAXIS)
        kmin = blkLimits(LOW,KAXIS)  
        jmax = blkLimits(HIGH,JAXIS)
        jmin = blkLimits(LOW,JAXIS)
        imax = blkLimits(HIGH,IAXIS)
        imin = blkLimits(LOW,IAXIS)

        do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
            if (NDIM == 3) zz = bndBox(LOW,KAXIS) + (k-kmin+0.5)*delta(KAXIS)
            do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                if (NDIM >= 2) yy = bndBox(LOW,JAXIS) + (j-jmin+0.5)*delta(JAXIS)
                do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                    if (solnData(DENS_VAR,i,j,k) .lt. sim_fluffDampCutoff) cycle
                    xx = bndBox(LOW,IAXIS) + (i-imin+0.5)*delta(IAXIS)
                    x = xx - grv_exactvec(1) - grv_ptvec(1) + grv_obvec(1)
                    y = yy - grv_exactvec(2) - grv_ptvec(2) + grv_obvec(2)
                    z = zz - grv_exactvec(3) - grv_ptvec(3) + grv_obvec(3)
                    r = max(dsqrt(x**2.d0 + y**2.d0 + z**2.d0), inner_rad**2.d0/dsqrt(x**2.d0 + y**2.d0 + z**2.d0))
                    tot_ener = tot_ener - dvol*solnData(DENS_VAR,i,j,k)*(newton*grv_ptmass/r)
                    if (-solnData(GPOT_VAR,i,j,k) .ge. 0.5d0*sum(solnData(VELX_VAR:VELZ_VAR,i,j,k)**2.d0)) then
                        ener = ener - dvol*solnData(DENS_VAR,i,j,k)*(newton*grv_ptmass/r)
                    endif
                enddo
            enddo
        enddo

        call Grid_releaseBlkPtr(blockList(lb), solnData)
        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
    enddo

    grv_ener = 0.d0
    call MPI_ALLREDUCE(ener, grv_ener, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    grv_tot_ener = 0.d0
    call MPI_ALLREDUCE(tot_ener, grv_tot_ener, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
end subroutine Orbit_energy
