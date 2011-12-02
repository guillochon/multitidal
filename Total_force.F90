subroutine Total_force(blockCount, blockList)
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords
    use Gravity_data, ONLY: grv_obvec, grv_ptvec, grv_obaccel, grv_ptaccel, grv_hptaccel, &
        grv_optaccel, grv_oobaccel, grv_oexactvec, grv_exactvec, grv_hobvec, grv_hptvec, &
        grv_oobvec, grv_optvec, grv_orb3D, grv_ptmass, grv_optmass
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Simulation_data, ONLY: sim_fluffDampCutoff, sim_softenRadius
    use gr_mpoleData, ONLY: twelfth
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use Grid_data, ONLY: gr_meshMe
    implicit none 
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ, lb, istat, i, j, k, ierr, it, potVar, denVar
    double precision,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    logical :: gcell = .true.
    double precision, dimension(:,:,:,:),pointer :: solnData
    double precision, dimension(4) :: lsum, gsum
    double precision, dimension(3) :: deld, offset, ptpos
    double precision :: dvol, delm, dr32, newton
    double precision :: tinitial, dx, delxinv, cell_mass, ptmass

    call RuntimeParameters_get('tinitial',tinitial)
    call PhysicalConstants_get("Newton", newton)

    do it = 1, 3
        lsum = 0.d0
        gsum = 0.d0

        if (it .eq. 1) then
            potVar = GPOL_VAR
            denVar = ODEN_VAR
            offset = grv_oexactvec(1:3) + grv_optvec(1:3) - grv_oobvec(1:3)
            ptmass = grv_optmass
        elseif (it .eq. 2) then
            offset = (grv_oexactvec(1:3) + grv_exactvec(1:3))/2.d0 + grv_hptvec(1:3) - grv_hobvec(1:3)
            ptmass = (grv_optmass + grv_ptmass) / 2.d0
        else
            potVar = GPOT_VAR
            denVar = DENS_VAR
            offset = grv_exactvec(1:3) + grv_ptvec(1:3) - grv_obvec(1:3)
            ptmass = grv_ptmass
        endif

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
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                        if (solnData(denVar,i,j,k) .lt. sim_fluffDampCutoff) cycle !Don't include any mass below the damping cutoff

                        cell_mass = solnData(denVar,i,j,k) * dvol
                        lsum(1) = lsum(1) + cell_mass

                        ptpos = (/ xCoord(i), yCoord(j), zCoord(k) /) - offset
                        dr32 = dsqrt(sum(ptpos**2.d0))
                        if (dr32 .lt. sim_softenRadius) then
                            dr32 = sim_softenRadius*sim_softenRadius*dr32
                        else
                            dr32 = dr32*dr32*dr32
                        endif
                        lsum(2:4) = lsum(2:4) - cell_mass * newton*ptmass/dr32*ptpos

                        if (it .eq. 2) cycle

                        if ((solnData(denVar,i+1,j,k) - solnData(denVar,i,j,k))*&
                            (solnData(denVar,i,j,k) - solnData(denVar,i-1,j,k)) .lt. 0.d0) then
                            deld(1) = sign(min(dabs(solnData(denVar,i-1,j,k) - solnData(denVar,i+1,j,k)), &
                                2.d0*dabs(solnData(denVar,i,j,k) - solnData(denVar,i-1,j,k)), &
                                2.d0*dabs(solnData(denVar,i,j,k) - solnData(denVar,i+1,j,k))), &
                                solnData(denVar,i+1,j,k) - solnData(denVar,i-1,j,k))
                        else
                            deld(1) = 0.d0
                        endif
                        if ((solnData(denVar,i,j+1,k) - solnData(denVar,i,j,k))*&
                            (solnData(denVar,i,j,k) - solnData(denVar,i,j-1,k)) .lt. 0.d0) then
                            deld(2) = sign(min(dabs(solnData(denVar,i,j-1,k) - solnData(denVar,i,j+1,k)), &
                                2.d0*dabs(solnData(denVar,i,j,k) - solnData(denVar,i,j-1,k)), &
                                2.d0*dabs(solnData(denVar,i,j,k) - solnData(denVar,i,j+1,k))), &
                                solnData(denVar,i,j+1,k) - solnData(denVar,i,j-1,k))
                        else
                            deld(2) = 0.d0
                        endif
                        if ((solnData(denVar,i,j,k+1) - solnData(denVar,i,j,k))*&
                            (solnData(denVar,i,j,k) - solnData(denVar,i,j,k-1)) .lt. 0.d0) then
                            deld(3) = sign(min(dabs(solnData(denVar,i,j,k-1) - solnData(denVar,i,j,k+1)), &
                                2.d0*dabs(solnData(denVar,i,j,k) - solnData(denVar,i,j,k-1)), &
                                2.d0*dabs(solnData(denVar,i,j,k) - solnData(denVar,i,j,k+1))), &
                                solnData(denVar,i,j,k+1) - solnData(denVar,i,j,k-1))
                        else
                            deld(3) = 0.d0
                        endif
                        lsum(5:7) = lsum(5:7) + cell_mass * delxinv * (/ &
                            solnData(potVar,i-1,j,k) - solnData(potVar,i+1,j,k) + deld(1)/solnData(denVar,i,j,k)*twelfth*&
                                (solnData(potVar,i-1,j,k) - 2.d0*solnData(potVar,i,j,k) + solnData(potVar,i+1,j,k)), &
                            solnData(potVar,i,j-1,k) - solnData(potVar,i,j+1,k) + deld(2)/solnData(denVar,i,j,k)*twelfth*&
                                (solnData(potVar,i,j-1,k) - 2.d0*solnData(potVar,i,j,k) + solnData(potVar,i,j+1,k)), &
                            solnData(potVar,i,j,k-1) - solnData(potVar,i,j,k+1) + deld(3)/solnData(denVar,i,j,k)*twelfth*&
                                (solnData(potVar,i,j,k-1) - 2.d0*solnData(potVar,i,j,k) + solnData(potVar,i,j,k+1)) /)
                    enddo
                enddo
            enddo

            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
        enddo

        call MPI_ALLREDUCE(lsum, gsum, 7, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

        if (it .eq. 1) then
            grv_optaccel = gsum(2:4) / gsum(1)
            if (.not. grv_orb3D) grv_optaccel(3) = 0.d0
            grv_oobaccel = gsum(5:7) / gsum(1)
            !grv_oobaccel = 0.d0
        elseif (it .eq. 2) then
            grv_hptaccel = gsum(2:4) / gsum(1)
            if (.not. grv_orb3D) grv_hptaccel(3) = 0.d0
        else
            grv_ptaccel = gsum(2:4) / gsum(1)
            if (.not. grv_orb3D) grv_ptaccel(3) = 0.d0
            grv_obaccel = gsum(5:7) / gsum(1)
            !grv_obaccel = 0.d0
        endif
    enddo

    !if (gr_meshMe .eq. MASTER_PE) then
    !    print *, 'grv_optaccel', grv_optaccel
    !    print *, 'grv_oobaccel', grv_oobaccel
    !    print *, 'grv_hptaccel', grv_hptaccel
    !    print *, 'grv_ptaccel', grv_ptaccel
    !    print *, 'grv_obaccel', grv_obaccel
    !    call Driver_abortFlash('done')
    !endif

end subroutine Total_force
