subroutine Total_force(blockCount, blockList)
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_findExtrema
    use Gravity_data, ONLY: grv_obvec, grv_ptvec, grv_obaccel, grv_ptaccel, grv_hptaccel, &
        grv_optaccel, grv_oobaccel, grv_oexactvec, grv_exactvec, grv_hobvec, grv_hptvec, &
        grv_oobvec, grv_optvec, grv_orb3D, grv_ptmass, grv_optmass, grv_densCut, &
        grv_comCutoff, grv_comPeakCut, grv_ompoleaccel, grv_mpoleaccel, grv_ompolevec, &
        grv_mpolevec, grv_totmass, grv_ototmass, grv_o2obaccel, grv_o2mpoleaccel, grv_mode
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Simulation_data, ONLY: sim_fluffDampCutoff, sim_softenRadius, sim_totForceInv, sim_totForceSub
    use gr_mpoleData, ONLY: twelfth
    use PhysicalConstants_interface, ONLY: PhysicalConstants_get
    use Grid_data, ONLY: gr_meshMe, gr_meshComm
    !use gr_mpoleInterface, ONLY: gr_mpoleGradPot, gr_mpoleGradOldPot
    use Driver_data, ONLY: dr_dt, dr_dtOld

    implicit none 
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ, lb, istat, i, j, k, ierr, it, potVar, denVar, ii, jj, kk
    double precision,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    logical :: gcell = .true.
    double precision, dimension(:,:,:,:),pointer :: solnData
    double precision, dimension(12) :: lsum, gsum
    double precision, dimension(3) :: deld, ooffset, hoffset, noffset, ptpos, tempaccel
    double precision :: dvol, dr32, newton, dtfac, dtfac2, ptmass, xx, yy, zz
    double precision :: tinitial, dx, delxinv, cell_mass, denscut, ldenscut, extrema

    call RuntimeParameters_get('tinitial',tinitial)
    call PhysicalConstants_get("Newton", newton)

    ldenscut = -huge(0.d0)
    do lb = 1,blockCount
       call Grid_findExtrema (blockList(lb), DENS_VAR, 1, extrema)
       if (extrema .gt. ldenscut) ldenscut = extrema
    enddo

    call MPI_ALLREDUCE(ldenscut, denscut, 1, FLASH_REAL, MPI_MAX, gr_meshComm, ierr)

    denscut = denscut*grv_comPeakCut

    dtfac = dr_dt/dr_dtOld
    dtfac2 = dtfac/2.d0

    ooffset = grv_exactvec(1:3) + grv_optvec(1:3) - grv_oobvec(1:3)
    hoffset = dtfac2*(grv_exactvec(1:3) - grv_oexactvec(1:3)) + grv_exactvec(1:3) + grv_hptvec(1:3) - grv_hobvec(1:3)
    noffset = dtfac *(grv_exactvec(1:3) - grv_oexactvec(1:3)) + grv_exactvec(1:3) + grv_ptvec(1:3) - grv_obvec(1:3)

    do it = 1, 3
        lsum = 0.d0
        gsum = 0.d0

        if (it .eq. 1) then
#ifdef GPO2_VAR
            potVar = GPO2_VAR
            denVar = ODE2_VAR
#else
            cycle
#endif
        elseif (it .eq. 2) then
            potVar = GPOL_VAR
            denVar = ODEN_VAR
        else
            potVar = GPOT_VAR
            denVar = DENS_VAR
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
                        if (solnData(denVar,i,j,k) .lt. sim_fluffDampCutoff) cycle

                        cell_mass = solnData(denVar,i,j,k) * dvol
                        lsum(1) = lsum(1) + cell_mass

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

                        lsum(2:4) = lsum(2:4) + cell_mass * delxinv * (/ &
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

        call MPI_ALLREDUCE(lsum(1:4), gsum(1:4), 4, FLASH_REAL, MPI_SUM, gr_meshComm, ierr)

        if (it .eq. 1) then
            grv_o2mpoleaccel = gsum(2:4) / gsum(1)
        elseif (it .eq. 2) then
            grv_ompoleaccel = gsum(2:4) / gsum(1)
        else
            grv_mpoleaccel = gsum(2:4) / gsum(1)
        endif
        !grv_o2mpoleaccel = 0.d0
        !grv_ompoleaccel = 0.d0
        !grv_mpoleaccel = 0.d0
    enddo

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
        dvol = (dx/sim_totForceSub)**3.d0
        delxinv = 0.5e0 / dx
        do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
            do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                    if (solnData(denVar,i,j,k) .lt. sim_fluffDampCutoff) cycle !Don't include any mass below the damping cutoff
                    do kk = 1, sim_totForceSub
                       zz    = zCoord(k) + dx*((2*kk - 1)*sim_totForceInv - 0.5)
                       do jj = 1, sim_totForceSub
                          yy    = yCoord(j) + dx*((2*jj - 1)*sim_totForceInv - 0.5)
                          do ii = 1, sim_totForceSub
                             xx    = xCoord(i) + dx*((2*ii - 1)*sim_totForceInv - 0.5)
                            
                             cell_mass = solnData(DENS_VAR,i,j,k) * dvol
                             ptpos = (/ xx, yy, zz /) - ooffset
                             ptmass = grv_ptmass
                             dr32 = dsqrt(sum(ptpos**2.d0))
                             if (dr32 .lt. sim_softenRadius) then
                                 dr32 = sim_softenRadius*sim_softenRadius*dr32
                             else
                                 dr32 = dr32*dr32*dr32
                             endif
                             lsum(1) = lsum(1) + cell_mass
                             lsum(2:4) = lsum(2:4) - cell_mass * newton*ptmass/dr32*ptpos

                             cell_mass = (dtfac2*(solnData(DENS_VAR,i,j,k) - solnData(ODEN_VAR,i,j,k)) + &
                                 solnData(DENS_VAR,i,j,k)) * dvol
                             ptpos = (/ xx, yy, zz /) - hoffset
                             ptmass = (dtfac2*(grv_ptmass - grv_optmass) + grv_ptmass)
                             dr32 = dsqrt(sum(ptpos**2.d0))
                             if (dr32 .lt. sim_softenRadius) then
                                 dr32 = sim_softenRadius*sim_softenRadius*dr32
                             else
                                 dr32 = dr32*dr32*dr32
                             endif
                             lsum(5) = lsum(5) + cell_mass
                             lsum(6:8) = lsum(6:8) - cell_mass * newton*ptmass/dr32*ptpos

                             cell_mass = (dtfac*(solnData(DENS_VAR,i,j,k) - solnData(ODEN_VAR,i,j,k)) + &
                                 solnData(DENS_VAR,i,j,k)) * dvol
                             ptpos = (/ xx, yy, zz /) - noffset
                             ptmass = (dtfac*(grv_ptmass - grv_optmass) + grv_ptmass)
                             dr32 = dsqrt(sum(ptpos**2.d0))
                             if (dr32 .lt. sim_softenRadius) then
                                 dr32 = sim_softenRadius*sim_softenRadius*dr32
                             else
                                 dr32 = dr32*dr32*dr32
                             endif
                             lsum(9) = lsum(9) + cell_mass
                             lsum(10:12) = lsum(10:12) - cell_mass * newton*ptmass/dr32*ptpos
                          enddo
                       enddo
                    enddo
                enddo
            enddo
        enddo

        call Grid_releaseBlkPtr(blockList(lb), solnData)
        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
    enddo

    call MPI_ALLREDUCE(lsum, gsum, 12, FLASH_REAL, MPI_SUM, gr_meshComm, ierr)

    grv_optaccel = gsum(2:4)   / gsum(1)
    grv_hptaccel = gsum(6:8)   / gsum(5)
    grv_ptaccel  = gsum(10:12) / gsum(9)

    !if (grv_mode .eq. 1) then
    !    call gr_mpoleGradPot(ooffset, grv_optaccel)
    !    grv_optaccel = grv_optaccel*grv_optmass/grv_ototmass
    !    call gr_mpoleGradOldPot(hoffset, grv_hptaccel)
    !    call gr_mpoleGradPot(hoffset, tempaccel)
    !    grv_hptaccel = (dtfac2*(tempaccel - grv_hptaccel) + tempaccel)*&
    !        (dtfac2*(grv_ptmass - grv_optmass) + grv_ptmass)/(dtfac2*(grv_totmass - grv_ototmass) + grv_totmass)
    !    call gr_mpoleGradOldPot(noffset, grv_ptaccel)
    !    call gr_mpoleGradPot(noffset, tempaccel)
    !    grv_ptaccel = (dtfac*(tempaccel - grv_ptaccel) + tempaccel)*&
    !        (dtfac*(grv_ptmass - grv_optmass) + grv_ptmass)/(dtfac*(grv_totmass - grv_ototmass) + grv_totmass)
    !else
    !    call gr_mpoleGradOldPot(ooffset, grv_optaccel)
    !    call gr_mpoleGradPot(noffset, grv_ptaccel)
    !    grv_optaccel = grv_optaccel*grv_optmass/grv_ototmass
    !    grv_ptaccel = grv_ptaccel*grv_ptmass/grv_totmass
    !    call gr_mpoleGradOldPot(hoffset, grv_hptaccel)
    !    call gr_mpoleGradPot(hoffset, tempaccel)
    !    grv_hptaccel = (grv_hptaccel + tempaccel) / 2.d0
    !    grv_hptaccel = grv_hptaccel*(grv_ptmass + grv_optmass)/(grv_totmass + grv_ototmass)
    !endif

    if (.not. grv_orb3D) then
        grv_optaccel(3) = 0.d0
        grv_hptaccel(3) = 0.d0
        grv_ptaccel(3) = 0.d0
    endif
end subroutine Total_force
