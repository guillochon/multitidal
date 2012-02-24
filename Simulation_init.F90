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
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Multispecies_interface, ONLY : Multispecies_getSumFrac, Multispecies_getSumInv, Multispecies_getAvg
    use Grid_data, ONLY : gr_globalMe
    use Logfile_interface, ONLY : Logfile_stampMessage
    use tree, ONLY : lrefine_max

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Multispecies.h"
#include "Eos.h"

    integer             :: i, ierr, j
    double precision    :: start_t
    double precision  order,mtot,rhoc,polyk,mu, &
                    x(np),y(np),yp(np), &
                    mass(np),ebind(np), &
                    rhom(np),zbeta(np),ztemp(np),exact(np), &
                    xsurf,ypsurf,cfl
    integer mode,iend
    character(len=100) :: logstr
    logical :: ionized

    call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
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
!    call RuntimeParameters_get('sim_objMass',sim_objMass)
    call RuntimeParameters_get('sim_objCoreMass',sim_objCoreMass)
    call RuntimeParameters_get('sim_accRadius',sim_accRadius)
    call RuntimeParameters_get('sim_accCoeff',sim_accCoeff)
    call RuntimeParameters_get('sim_objPolyN',sim_objPolyN)
    call RuntimeParameters_get('sim_objPolyN2',sim_objPolyN2)
!    call RuntimeParameters_get('sim_objCentDen',sim_objCentDen)
    call RuntimeParameters_get('sim_objcore',sim_objCore)
    call RuntimeParameters_get('sim_objenve',sim_objEnve)
    call RuntimeParameters_get('sim_objRadius',sim_objRadius)
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
    call RuntimeParameters_get("ptmass", sim_ptMass)
    call RuntimeParameters_get("cfl", cfl)

    gammac = 1.0d0+1.0d0/sim_objPolyN
    gammae = 1.0d0+1.0d0/sim_objPolyN2
    if (gr_globalMe .eq. MASTER_PE) then
        obj_xn(H1_SPEC) = 0.7
        obj_xn(HE4_SPEC) = 0.3
        obj_xn(CORE_SPEC) = 0.0
        speciesMask(H1_SPEC) = H1_SPEC
        speciesMask(HE4_SPEC) = HE4_SPEC
        speciesMask(CORE_SPEC) = UNDEFINED_INT
        call Multispecies_getSumInv(A, obj_mue, obj_xn, speciesMask)
        obj_mue = 1.e0 / obj_mue
        print*, "obj_mue", obj_mue
        obj_xn(H1_SPEC) = 0.0
        obj_xn(HE4_SPEC) = 0.0
        obj_xn(CORE_SPEC) = 1.0
        call Multispecies_getSumInv(A, obj_muc, obj_xn) !, speciesMask)
        obj_muc = 1.e0 / obj_muc
        print*, "obj_muc", obj_muc
        mode = 1
        open(42,file="input.dat",position='append')
        write(42,*) sim_objPolyN,sim_objPolyN2,sim_objCoreMass,sim_objCentDen,polyk,sim_objRadius,obj_muc,obj_mue,mode,&
!                    x,y,yp,obj_radius,obj_rhop,mass,obj_prss,ebind, &
!                    rhom,ztemp,zbeta,exact,
                    xsurf,ypsurf,np,iend,obj_ipos,sim_objCore,sim_objEnve,obj_ipoi
!        close(42)
        call polytr(sim_objPolyN,sim_objPolyN2,sim_objMass,sim_objCoreMass,sim_objCentDen,polyk,sim_objRadius,obj_muc,obj_mue,mode, &
            x,y,yp,obj_radius,obj_rhop,mass,obj_prss,ebind, &
            rhom,ztemp,zbeta,exact,xsurf,ypsurf,np,iend,obj_ipos,sim_objCore,sim_objEnve,obj_ipoi)

        write(42,*) iend,obj_ipos,obj_ipoi,sim_objMass
		
		j = 0
		i = obj_ipoi + 1
		do while ((obj_radius(i+1)-obj_radius(i))/obj_radius(i) .lt. 2.3e-3)
			j = j + 1
			i = i + 1
		enddo

		do i = obj_ipoi + 1, iend
			if (i .le. (np-j)) then
				obj_radius(i) = obj_radius(i+j)
				obj_rhop(i)   = obj_rhop(i+j)
				obj_prss(i)   = obj_prss(i+j)
!!			else
!!				obj_radius(i) = 0.0
!!				obj_rhop(i)   = 0.0
!!				obj_prss(i)   = 0.0
			endif
		enddo
		iend	 = iend - j
		obj_ipos = obj_ipos - j

		do i=1,obj_ipos
			write(42,*) obj_radius(i)
		enddo
        write(42,*) iend,obj_ipos,obj_ipoi
        close(42)

        write(logstr, fmt='(A30, 2ES15.8)') 'Ambient CFL timestep:', cfl*(sim_xMax - sim_xMin)/NXB/2.d0**(lrefine_max-1)/dsqrt(sim_fluidGammae*sim_pAmbient/sim_rhoAmbient)
        call Logfile_stampMessage(logstr)
        write(logstr, fmt='(A30, 2ES15.8)') 'Fluff CFL timestep:', cfl*(sim_xMax - sim_xMin)/NXB/2.d0**(lrefine_max-1)/dsqrt(sim_fluidGammae*sim_smallP/sim_smallRho)
        call Logfile_stampMessage(logstr)
    endif

    call MPI_BCAST(obj_xn, NSPECIES, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    call MPI_BCAST(obj_muc, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    call MPI_BCAST(obj_mue, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    call MPI_BCAST(obj_radius, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)              
    call MPI_BCAST(obj_rhop, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    call MPI_BCAST(obj_prss, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)                
    call MPI_BCAST(obj_ipos, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)              
    call MPI_BCAST(obj_ipoi, 1, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(gammac, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(gammae, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(sim_objMass, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

    !Sink radius is scaled relative to the original pericenter distance
    
    sim_softenRadius = sim_softenRadius*obj_radius(obj_ipos)/sim_periBeta*(sim_ptMass/sim_objMass/sim_msun)**(1.d0/3.d0)
    sim_accRadius = sim_accRadius*obj_radius(obj_ipos)/sim_periBeta*(sim_ptMass/sim_objMass/sim_msun)**(1.d0/3.d0)
    sim_startDistance = obj_radius(obj_ipos)/sim_startBeta*(sim_ptMass/sim_objMass/sim_msun)**(1.d0/3.d0)

    write(logstr, fmt='(A30, ES15.8)') 'Sink radius:', sim_softenRadius
    call Logfile_stampMessage(logstr)

    sim_inSubInv = 1.d0/2.d0**sim_nSubZones
    sim_nSubZones = 2**(sim_nSubZones - 1)
    sim_inszd      = (1.d0/sim_nSubZones)**NDIM

    sim_totForceInv = 1.d0/2.d0**sim_totForceSub
    sim_totForceSub = 2**(sim_totForceSub - 1)

end subroutine Simulation_init
