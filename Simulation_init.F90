!!****if* source/Simulation/SimulationMain/WDAccrection/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!! ARGUMENTS
!!
!!   myPE -   current processor number
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***

subroutine Simulation_init()

    use Simulation_data 
    use Particles_sinkData, ONLY : particles_local, ipm, ipvx, ipvy, ipvz
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Multispecies_interface, ONLY : Multispecies_getSumFrac, Multispecies_getSumInv, Multispecies_getAvg
    use Grid_data, ONLY : gr_globalMe
    use Logfile_interface, ONLY : Logfile_stampMessage
    use tree, ONLY : lrefine_max
    use Eos_data, ONLY : eos_gasConstant
    use pt_sinkInterface, ONLY : pt_sinkCreateParticle
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use Driver_data, ONLY : dr_simTime

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Multispecies.h"
#include "Eos.h"
#include "Starprof.h"

    integer             :: ierr, i
    double precision  polyk, &
                    x(np),y(np),yp(np), &
                    mass(np),ebind(np), &
                    rhom(np),zbeta(np),ztemp(np),exact(np), &
                    xsurf,ypsurf
    double precision :: cfl, rho0, a, period, start_dist, ecc_anom, newton
    integer mode,iend,pno
    character(len=200) :: logstr
    double precision, dimension(6) :: obvec, ptvec

    call PhysicalConstants_get("Newton", newton)

    call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
    call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
    call RuntimeParameters_get('smallx', sim_smallX)
    call RuntimeParameters_get('smlrho', sim_smallRho)
    call RuntimeParameters_get('smallp', sim_smallP)
    call RuntimeParameters_get('xmin',sim_xMin)
    call RuntimeParameters_get('ymin',sim_yMin)
    call RuntimeParameters_get('zmin',sim_zMin)
    call RuntimeParameters_get('xmax',sim_xMax)
    call RuntimeParameters_get('ymax',sim_yMax)
    call RuntimeParameters_get('zmax',sim_zMax)
    call RuntimeParameters_get('tinitial',sim_tInitial)
    call RuntimeParameters_get('sim_tRelax',sim_tRelax)
    call RuntimeParameters_get('sim_relaxRate',sim_relaxRate)
    call RuntimeParameters_get('sim_softenRadius',sim_softenRadius)
    call RuntimeParameters_get('sim_accRadius',sim_accRadius)
    call RuntimeParameters_get('sim_accCoeff',sim_accCoeff)
    call RuntimeParameters_get('sim_fluffDampCoeff',sim_fluffDampCoeff)
    call RuntimeParameters_get('sim_fluffDampCutoff',sim_fluffDampCutoff)
    call RuntimeParameters_get('sim_maxBlocks',sim_maxBlocks)
    call RuntimeParameters_get("sim_periBeta", sim_periBeta)
    call RuntimeParameters_get("sim_startBeta", sim_startBeta)
    call RuntimeParameters_get("sim_periodFac", sim_periodFac)
    call RuntimeParameters_get("sim_orbEcc", sim_orbEcc)
    call RuntimeParameters_get("sim_useInitialPeakDensity", sim_useInitialPeakDensity)
    call RuntimeParameters_get("sim_ptMassRefine", sim_ptMassRefine)
    call RuntimeParameters_get("sim_ptMassRefRad", sim_ptMassRefRad)
    call RuntimeParameters_get("sim_totForceSub", sim_totForceSub)
    call RuntimeParameters_get("sim_rotFac", sim_rotFac)
    call RuntimeParameters_get("sim_rotAngle", sim_rotAngle)
    call RuntimeParameters_get("sim_tSpinup", sim_tSpinup)

    call RuntimeParameters_get("sim_powerLawScale", sim_powerLawScale)
    call RuntimeParameters_get("sim_powerLawExponent", sim_powerLawExponent)
    call RuntimeParameters_get("sim_powerLawExtent", sim_powerLawExtent)
    call RuntimeParameters_get("sim_powerLawMass", sim_powerLawMass)
    call RuntimeParameters_get("sim_powerLawTemperature", sim_powerLawTemperature)

    call RuntimeParameters_get("sim_ptMass", sim_ptMass)

    if (sim_powerLawExponent .le. -3.d0) then
        call Driver_abortFlash('Error: sim_powerLawExponent must be greater than -3.0')
    endif

    sim_xCenter = (sim_xMax + sim_xMin) / 2.d0
    sim_yCenter = (sim_yMax + sim_yMin) / 2.d0
    sim_zCenter = (sim_zMax + sim_zMin) / 2.d0

    call RuntimeParameters_get('sim_kind',sim_kind)

#ifdef LOADPROFILE
    call RuntimeParameters_get('sim_tAmbient', sim_tAmbient)
    call RuntimeParameters_get('sim_profFile',sim_profFile)

    if (gr_globalMe .eq. MASTER_PE) then
        call read_table_dims(sim_profFile, sim_tableRows, sim_tableCols)
        allocate(sim_table(sim_tableRows, sim_tableCols))
        call read_table(sim_profFile, sim_table, sim_tableRows, sim_tableCols)
        write(*,*) 'r_prof', sim_table(:,R_PROF)
        write(*,*) 'rho_prof', sim_table(:,RHO_PROF)
        write(*,*) 'temp_prof', sim_table(:,TEMP_PROF)
    endif

    call MPI_BCAST(sim_tableRows, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(sim_tableCols, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    if (gr_globalMe .ne. MASTER_PE) allocate(sim_table(sim_tableRows, sim_tableCols))
    call MPI_BCAST(sim_table, sim_tableRows*sim_tableCols, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

    sim_objRadius = sim_table(sim_tableRows,R_PROF)
    sim_objMass = sum(PI*(sim_table(2:sim_tableRows,R_PROF) + sim_table(1:sim_tableRows-1,R_PROF))**2.d0*&
                      (sim_table(2:sim_tableRows,R_PROF) - sim_table(1:sim_tableRows-1,R_PROF))*&
                      sim_table(:,RHO_PROF))
    sim_objCentDens = sim_table(1,RHO_PROF)
#else
    call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
    call RuntimeParameters_get('sim_objPolyN',sim_objPolyN)
    call RuntimeParameters_get('sim_objMass',sim_objMass)
    call RuntimeParameters_get('sim_objCentDens',sim_objCentDens)

    obj_xn = 0.d0
    obj_xn(H1_SPEC) = 0.7d0
    obj_xn(HE4_SPEC) = 0.3d0

    call Multispecies_getSumInv(A, obj_mu, obj_xn)
    obj_mu = 1.e0 / obj_mu

    if (sim_kind .eq. 'polytrope') then
        if (gr_globalMe .eq. MASTER_PE) then
            mode = 1
            call polytr(sim_objPolyN,sim_objMass/sim_msun,sim_objCentDens,polyk,obj_mu,mode, &
                x,y,yp,obj_radius,obj_rhop,mass,obj_prss,ebind, &
                rhom,ztemp,zbeta,exact,xsurf,ypsurf,np,iend,obj_ipos)
        endif
    elseif (sim_kind .eq. 'powerlaw') then
        if (gr_globalMe .eq. MASTER_PE) then
            obj_ipos = np
            rho0 = sim_powerLawMass*sim_powerLawExtent**(-3.d0 - sim_powerLawExponent)*&
                sim_powerLawScale**sim_powerLawExponent*(3.d0 + sim_powerLawExponent)/(2.d0*PI)
            do i = 1, np
                obj_radius(i) = sim_powerLawExtent*dble(i)/np
                obj_rhop(i) = rho0*(obj_radius(i)/sim_powerLawScale)**sim_powerLawExponent
                obj_prss(i) = eos_gasConstant*obj_rhop(i) * &
                                 sim_powerLawTemperature / obj_mu
            enddo 
            sim_objMass = sim_powerLawMass
            sim_objCentDens = obj_rhop(1)
        endif
    endif

    !call MPI_BCAST(obj_xn, NSPECIES, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    !call MPI_BCAST(obj_mu, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    call MPI_BCAST(obj_radius, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    call MPI_BCAST(obj_rhop, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    call MPI_BCAST(obj_prss, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    call MPI_BCAST(obj_ipos, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)              

    sim_objRadius = obj_radius(obj_ipos)
#endif

    if (gr_globalMe .eq. MASTER_PE) then
        write(logstr, fmt='(A30, ES15.8)') 'Object mass (m_sun):', sim_objMass/sim_msun
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, ES15.8)') 'Object radius:', sim_objRadius
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, ES15.8)') 'Object central density:', sim_objCentDens
        call Logfile_stampMessage(logstr)

        call RuntimeParameters_get("cfl", cfl)

        write(logstr, fmt='(A30, ES15.8)') 'Ambient CFL timestep:', cfl*(sim_xMax - sim_xMin)/NXB/2.d0**(lrefine_max-1)/dsqrt(sim_fluidGamma*sim_pAmbient/sim_rhoAmbient)
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, ES15.8)') 'Fluff CFL timestep:', cfl*(sim_xMax - sim_xMin)/NXB/2.d0**(lrefine_max-1)/dsqrt(sim_fluidGamma*sim_smallP/sim_smallRho)
        call Logfile_stampMessage(logstr)
    endif

    !Sink radius is scaled relative to the original pericenter distance
    
    sim_softenRadius = sim_softenRadius*sim_objRadius/sim_periBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
    sim_accRadius = sim_accRadius*sim_objRadius/sim_periBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
    sim_startDistance = sim_objRadius/sim_startBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)

    write(logstr, fmt='(A30, ES15.8)') 'Sink radius:', sim_softenRadius
    call Logfile_stampMessage(logstr)

    sim_inSubInv = 1.d0/2.d0**sim_nSubZones
    sim_nSubZones = 2**(sim_nSubZones - 1)
    sim_inszd      = (1.d0/sim_nSubZones)**NDIM

    sim_totForceInv = 1.d0/2.d0**sim_totForceSub
    sim_totForceSub = 2**(sim_totForceSub - 1)

    if (sim_periBeta .eq. 0.d0) then
        call RuntimeParameters_get("sim_periDist", sim_periDist)
    else
        sim_periDist = sim_objRadius/sim_periBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
    endif
    if (sim_periodFac .gt. 0.d0) then
        a = sim_periDist/(1.d0 - sim_orbEcc)
        period = dsqrt(4.d0*PI**2.d0/newton/(sim_ptMass + sim_objMass)*a**3.d0)
        sim_periTime = sim_periodFac*period + sim_tRelax
    endif
    if (sim_startBeta .gt. 0.d0) then
        a = sim_periDist/(1.d0 - sim_orbEcc)
        period = dsqrt(4.d0*PI**2.d0/newton/(sim_ptMass + sim_objMass)*a**3.d0)
        start_dist = sim_objRadius/sim_startBeta*(sim_ptMass/sim_objMass)**(1.d0/3.d0)
        if (start_dist .gt. 2.d0*a - sim_periDist) call Driver_abortFlash('start_dist too large!')
        ecc_anom = dacos((a - start_dist)/a/sim_orbEcc)
        sim_periTime = abs(ecc_anom - sim_orbEcc*dsin(ecc_anom))*period/2.d0/PI + sim_tRelax
    endif

    if (gr_globalMe .eq. MASTER_PE) then
        call calc_orbit(0.d0, sim_objMass, sim_ptMass, obvec, ptvec)
        ptvec = ptvec - obvec
        pno = pt_sinkCreateParticle(0.5d0*sim_xmax + ptvec(1), 0.5d0*sim_ymax + ptvec(2), &
            0.5d0*sim_zmax + ptvec(3), 0., 1, gr_globalMe)
        particles_local(ipm,pno) = sim_ptMass
        particles_local(ipvx,pno) = ptvec(4)
        particles_local(ipvy,pno) = ptvec(5)
        particles_local(ipvz,pno) = ptvec(6)
        pno = pt_sinkCreateParticle(0.5d0*sim_xmax, 0.5d0*sim_ymax, 0.5d0*sim_zmax, 0., 1, gr_globalMe)
        particles_local(ipm,pno) = 2.d33
    endif

    if (gr_globalMe .eq. MASTER_PE) then
        write(logstr, fmt='(A30, 2ES15.8)') 'Start distance:', start_dist
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 6ES15.8)') 'Pt. mass start pos:', ptvec
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 2ES15.8)') 'Semimajor axis, eccentricity:', a, sim_orbEcc
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 2ES15.8)') 'Period, pericenter time:', period, sim_periTime
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 2ES15.8)') 'Obj. radius, pericenter dist:', sim_objRadius, sim_periDist
        call Logfile_stampMessage(logstr)

        if (dr_simTime .eq. sim_tInitial) then
            open(unit = 11, file = 'extras.dat', status = 'unknown')
            write(11, fmt='(ES22.15)') sim_periBeta
            write(11, fmt='(ES22.15)') period
            write(11, fmt='(ES22.15)') sim_periTime
            write(11, fmt='(ES22.15)') sim_objRadius
            write(11, fmt='(ES22.15)') sim_objMass
            write(11, fmt='(ES22.15)') sim_orbEcc
            write(11, fmt='(ES22.15)') sim_ptMass
            write(logstr, fmt='(I4)')  lrefine_max
            write(11, fmt='(A)') adjustl(logstr)
            write(11, fmt='(ES22.15)') sim_xMax - sim_xMin
            close(11) 
        endif
    endif
end subroutine Simulation_init
