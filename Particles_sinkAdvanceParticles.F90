!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkAdvanceParticles
!!
!! NAME
!!
!!  Particles_sinkAdvanceParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAdvanceParticles(real, INTENT(in)  :: dr_dt)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   dr_dt : 
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

! used to be AdvanceSinkParticles
subroutine Particles_sinkAdvanceParticles(dr_dt)

  ! Advance sink particles through timestep dr_dt

  use Particles_sinkData
  use pt_sinkInterface, only: pt_sinkAccelGasOnSinks, pt_sinkAccelSinksOnSinks, &
      pt_sinkGetSubCycleTimeStep
  use Particles_data, only: pt_indexCount, pt_indexList
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_globalMe
  use Cosmology_interface, ONLY :  Cosmology_getParams, Cosmology_getRedshift, Cosmology_getOldRedshift
  use Driver_interface, ONLY : Driver_getSimTime, Driver_abortFlash
  use Grid_interface, ONLY : Grid_moveParticles, Grid_sortParticles
  ! Added by JFG
  use pt_sinkInterface, only: pt_sinkGatherGlobal
  use Simulation_data, only: sim_fixedPartTag
  ! End JFG
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  implicit none

  real, INTENT(in) :: dr_dt

  real, dimension(maxsinks) :: ax_gas, ay_gas, az_gas

  integer             :: i, nsubcycles, j
  logical, save       :: first_call = .true.
  real                :: dt, dt_global, t_sub, wterm
  real, save          :: dt_subOld
  real                :: local_min_radius, local_max_accel
  logical             :: end_subcycling
  integer, save       :: integrator
  logical, parameter  :: debug = .false.
  real                :: max_part_vel
  character(len=80)   :: sink_integrator
  real                :: alpha, dotalpha, aterm, bterm, woldterm
  real, save          :: Hubble, OmegaM, OmegaCurv, OmegaL, OmegaB
  real                :: redshift_old, sOld, a0, simTime
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  
  ! this should be an interface
  real                :: pt_sinkBNCOEFF, pt_sinkCNCOEFF, pt_sinkDNCOEFF
  
  ! Added by JFG
  real, dimension(6)  :: fixed_vec
  integer             :: fixedi
  ! End JFG

  interface
    subroutine pt_sinkCsm_integrateFriedmann(a,dt,anew)
      real, intent(IN)  :: a
      real, intent(IN)  :: dt
      real, intent(OUT) :: anew
    end subroutine
  end interface

  call Cosmology_getOldRedshift(redshift_old)
  a0 = 1./(1. + redshift_old)

  if (first_call) then

     call RuntimeParameters_get("sink_integrator", sink_integrator)
     select case (sink_integrator)
     case ("leapfrog")
        integrator = 1
     case ("euler")
        integrator = 2
     case ("leapfrog_cosmo")
        integrator = 3
     case default
        sink_integrator = "leapfrog"
        integrator = 1
        print*, "No sink integrator specificed in flash.par..."
     end select
     if (dr_globalMe .eq. MASTER_PE) print*, "Advance Particles, sink integrator = ", trim(sink_integrator)

     call Cosmology_getParams(Hubble,OmegaM,OmegaB,OmegaL)
     OmegaCurv = OmegaM + OmegaL - 1.

     first_call = .false.
  end if

  if (localnpf .eq. 0) return

  ! clear particle accelerations
  do i=1, localnp
     particles_local(ACCX_PART_PROP,i) = 0.0
     particles_local(ACCY_PART_PROP,i) = 0.0
     particles_local(ACCZ_PART_PROP,i) = 0.0
  end do

  if(debug) then
     do i=1, localnp
        do j=1, pt_sinkParticleProps
           print*, i,j,particles_local(j,i)
        end do
     end do
  end if

  ! Compute acceleration on sinks from gas and 
  ! mapped particle (not sinks) density (e.g. dark matte)

  if(debug) then
     max_part_vel = 0.0
    do i=1, localnp
       max_part_vel = max(max_part_vel, abs(particles_local(VELX_PART_PROP, i)))
       max_part_vel = max(max_part_vel, abs(particles_local(VELY_PART_PROP, i)))
       max_part_vel = max(max_part_vel, abs(particles_local(VELZ_PART_PROP, i)))
    end do
    print*, dr_globalMe, "start of advance particles (1): max_part_vel=", max_part_vel
  end if

  ! Added by JFG
  if (sim_fixedPartTag .ne. 0) then
      call pt_sinkGatherGlobal()
      do i = 1, localnpf
          if (idnint(particles_global(TAG_PART_PROP,i)) .eq. sim_fixedPartTag) then
              fixedi = i
              exit
          else
              cycle
          endif
          call Driver_abortFlash('Error: Unable to find fixed particle tag [1]')
      enddo
      fixed_vec(1) = particles_global(POSX_PART_PROP,fixedi)
      fixed_vec(2) = particles_global(POSY_PART_PROP,fixedi)
      fixed_vec(3) = particles_global(POSZ_PART_PROP,fixedi)
      fixed_vec(4) = particles_global(VELX_PART_PROP,fixedi)
      fixed_vec(5) = particles_global(VELY_PART_PROP,fixedi)
      fixed_vec(6) = particles_global(VELZ_PART_PROP,fixedi)
  endif
  ! End JFG

  call pt_sinkAccelGasOnSinks()

  if(debug) then
     max_part_vel = 0.0
    do i=1, localnp
       max_part_vel = max(max_part_vel, abs(particles_local(VELX_PART_PROP, i)))
       max_part_vel = max(max_part_vel, abs(particles_local(VELY_PART_PROP, i)))
       max_part_vel = max(max_part_vel, abs(particles_local(VELZ_PART_PROP, i)))
    end do
    print*, dr_globalMe, "after gas on sinks: max_part_vel=", max_part_vel
  end if

  ! store acceleration from GasOnSinksAccel
  do i=1, localnp
     ax_gas(i) = particles_local(ACCX_PART_PROP, i)
     ay_gas(i) = particles_local(ACCY_PART_PROP, i)
     az_gas(i) = particles_local(ACCZ_PART_PROP, i)
  end do

  dt_global = dr_dt
  t_sub = 0.0
  nsubcycles = 0
  end_subcycling = .false.

  if (integrator .eq. 1) then   !!!!!!! Leapfrog

     do while (.not. end_subcycling)

        do i = 1, localnp
           particles_local(ACCX_PART_PROP, i) = ax_gas(i)
           particles_local(ACCY_PART_PROP, i) = ay_gas(i)
           particles_local(ACCZ_PART_PROP, i) = az_gas(i)
        end do

        ! compute sink-sink gravity
        call pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

        if(debug) then
           max_part_vel = 0.0
           do i=1, localnp
              max_part_vel = max(max_part_vel, abs(particles_local(VELX_PART_PROP, i)))
              max_part_vel = max(max_part_vel, abs(particles_local(VELY_PART_PROP, i)))
              max_part_vel = max(max_part_vel, abs(particles_local(VELZ_PART_PROP, i)))
           end do
           print*, dr_globalMe, "after sinks on sinks: max_part_vel=", max_part_vel
           if(max_part_vel .ge. 1.0e10) then
              call Driver_abortFlash("velocity too high!")
           end if
           do i=1, localnp
              print*, "accleration contribution from other sinks x", i, particles_local(ACCX_PART_PROP, i) - ax_gas(i)
              print*, "accleration contribution from other sinks y", i, particles_local(ACCY_PART_PROP, i) - ay_gas(i)
              print*, "accleration contribution from other sinks z", i, particles_local(ACCZ_PART_PROP, i) - az_gas(i)
           end do
           print*, dr_globalMe, "advance sinks: getting sub cycle timestep"
           print*, dr_globalMe, "local_min_radius=", local_min_radius
           print*, dr_globalMe, "local_max_accel=", local_max_accel
        endif

        call pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)

        ! check whether we have reached the global (hydro) timestep for exiting
        if (t_sub + dt .ge. dt_global) then
           dt = dt_global - t_sub
           end_subcycling = .true.
        end if

        if(debug) then
           print*, dr_globalMe, "advance sink particles subcycle dt = ", dt
           print*, dr_globalMe, "advance sink particles, global dt = ", dt_global
        end if

        ! First step of leapfrog:
        wterm = 0.5 * dt

        do i = 1, localnp
           particles_local(VELX_PART_PROP, i) =  & 
                particles_local(VELX_PART_PROP, i) + & 
                wterm * particles_local(ACCX_PART_PROP, i)
           particles_local(POSX_PART_PROP, i) = & 
                particles_local(POSX_PART_PROP, i) + & 
                dt*particles_local(VELX_PART_PROP, i)

           particles_local(VELY_PART_PROP, i) =  & 
                particles_local(VELY_PART_PROP, i) + & 
                wterm * particles_local(ACCY_PART_PROP, i)
           particles_local(POSY_PART_PROP, i) = & 
                particles_local(POSY_PART_PROP, i) + & 
                dt*particles_local(VELY_PART_PROP, i)

           particles_local(VELZ_PART_PROP, i) =  & 
                particles_local(VELZ_PART_PROP, i) + & 
                wterm * particles_local(ACCZ_PART_PROP, i)
           particles_local(POSZ_PART_PROP, i) = & 
                particles_local(POSZ_PART_PROP, i) + & 
                dt*particles_local(VELZ_PART_PROP, i)
        enddo

        ! restore gas contribution
        do i = 1, localnp
           particles_local(ACCX_PART_PROP, i) = ax_gas(i)
           particles_local(ACCY_PART_PROP, i) = ay_gas(i)
           particles_local(ACCZ_PART_PROP, i) = az_gas(i)
        end do

        ! Add sink-sink acceleration at new position of particle
        if(debug) print*, dr_globalMe, "advance sinks: calculating sink on sink acceleration"
        call pt_sinkAccelSinksOnSinks(local_min_radius,local_max_accel)

        ! Second step of leapfrog
        do i = 1, localnp

           particles_local(VELX_PART_PROP, i) = & 
                particles_local(VELX_PART_PROP, i) + & 
                wterm*particles_local(ACCX_PART_PROP, i)

           particles_local(VELY_PART_PROP, i) = & 
                particles_local(VELY_PART_PROP, i) + & 
                wterm*particles_local(ACCY_PART_PROP, i)

           particles_local(VELZ_PART_PROP, i) = & 
                particles_local(VELZ_PART_PROP, i) + & 
                wterm*particles_local(ACCZ_PART_PROP, i)

        end do

        ! add dt to the current time within subcycling loop
        t_sub = t_sub + dt
        nsubcycles = nsubcycles + 1

        if(debug) then 
           print*, dr_globalMe, "nsubcycles=", nsubcycles
        endif

     end do

     if(debug) then
        if(dr_globalMe .eq. MASTER_PE) then
           print*, "subcycling on master PE"
        end if
     endif

  end if    ! leapfrog integrator

  if(integrator .eq. 2) then       ! euler method integrator.

     do while(.not. end_subcycling)

        do i = 1, localnp
           particles_local(ACCX_PART_PROP,i) = ax_gas(i)
           particles_local(ACCY_PART_PROP,i) = ay_gas(i)
           particles_local(ACCZ_PART_PROP,i) = az_gas(i)
        enddo

        call pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

        call pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)

        if (t_sub + dt .ge. dt_global) then
           dt = dt_global - t_sub
           end_subcycling = .true.
        endif

        do i = 1, localnp

           particles_local(POSX_PART_PROP,i) = particles_local(POSX_PART_PROP,i) + & 
                dt * particles_local(VELX_PART_PROP,i)
           particles_local(VELX_PART_PROP,i) = particles_local(VELX_PART_PROP,i) + & 
                dt * particles_local(ACCX_PART_PROP,i)

           particles_local(POSY_PART_PROP,i) = particles_local(POSY_PART_PROP,i) + & 
                dt * particles_local(VELY_PART_PROP,i)
           particles_local(VELY_PART_PROP,i) = particles_local(VELY_PART_PROP,i) + & 
                dt * particles_local(ACCY_PART_PROP,i)

           particles_local(POSZ_PART_PROP,i) = particles_local(POSZ_PART_PROP,i) + & 
                dt * particles_local(VELZ_PART_PROP,i)
           particles_local(VELZ_PART_PROP,i) = particles_local(VELZ_PART_PROP,i) + & 
                dt * particles_local(ACCZ_PART_PROP,i)

        enddo

        t_sub = t_sub + dt
        nsubcycles = nsubcycles + 1

     enddo

  end if    ! euler integrator

  if (integrator .eq. 3) then   ! Leapfrog Cosmo

     call Driver_getSimTime(simTime)
     ! simTime should be time at start of global timestep
     simTime = simTime - dt_global
     sOld = a0

     do while (.not. end_subcycling)

        do i = 1, localnp
           particles_local(ACCX_PART_PROP, i) = ax_gas(i)
           particles_local(ACCY_PART_PROP, i) = ay_gas(i)
           particles_local(ACCZ_PART_PROP, i) = az_gas(i)
        end do

        ! compute sink-sink gravity
        if(debug) print*, dr_globalMe, "advance sinks: computing sink on sink acceleration"
        call pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

        call pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)

        if(debug) print*, dr_globalMe, "advance sinks: got subcycle timestep", dt

        ! check whether we have reached the global (hydro) timestep for exiting
        if (t_sub + dt .ge. dt_global) then
           dt = dt_global - t_sub
           end_subcycling = .true.
        end if

        alpha    = 2.0*Hubble/sOld * sqrt(OmegaM/sOld-OmegaCurv+OmegaL*sOld**2)
        dotalpha = 2.0*Hubble**2 * (-3.0*OmegaM + 2.0*OmegaCurv*sOld)/sOld**3

        call pt_sinkCsm_integrateFriedmann(a0,dt,sOld)
        ! sOld is now the scale factor at start of subcycle timestep
        a0 = sOld
        simTime = simTime + dt

        ! Advances particle positions fully over the sink particle
        ! subcycle timestep
        ! Velocities are defined in center of timestep
        do i = 1, localnp

           dt_subOld = particles_local(ipdtold,i)

           ! If dt_subOld = 0, assume this is the first time moving this
           ! particular sink particle.
           ! As with dark matter, first step is defined at t0
           ! afterwards, position and velocity are staggered, ala leapfrog
           if(dt_subOld .eq. 0.0) then
              aterm     = dt
              bterm     = 1. - alpha*dt*0.5
              wterm     = 0.5*dt
              woldterm  = 0.
           else
              aterm      = dt
              bterm      = pt_sinkBNCOEFF(dt, dt_subOld, alpha, dotalpha)
              wterm      = pt_sinkCNCOEFF(dt, dt_subOld, alpha, dotalpha)
              woldterm   = pt_sinkDNCOEFF(dt, dt_subOld, alpha, dotalpha)
           end if

           particles_local(VELX_PART_PROP,i) = bterm*particles_local(VELX_PART_PROP,i) + &
                wterm*particles_local(ACCX_PART_PROP,i) + &
                woldterm*particles_local(OACX_PART_PROP,i)
           particles_local(POSX_PART_PROP,i) = particles_local(POSX_PART_PROP,i) + &
                aterm*particles_local(VELX_PART_PROP,i)
           particles_local(OACX_PART_PROP,i) = particles_local(ACCX_PART_PROP,i)

           particles_local(VELY_PART_PROP,i)= bterm*particles_local(VELY_PART_PROP,i) + &
                wterm*particles_local(ACCY_PART_PROP,i) + &
                woldterm*particles_local(OACY_PART_PROP,i)
           particles_local(POSY_PART_PROP,i)=  particles_local(POSY_PART_PROP,i) + &
                aterm*particles_local(VELY_PART_PROP,i)
           particles_local(OACY_PART_PROP,i) = particles_local(ACCY_PART_PROP,i)

           particles_local(VELZ_PART_PROP,i)= bterm*particles_local(VELZ_PART_PROP,i) + &
                wterm*particles_local(ACCZ_PART_PROP,i) + &
                woldterm*particles_local(OACZ_PART_PROP,i)
           particles_local(POSZ_PART_PROP,i) =  particles_local(POSZ_PART_PROP,i) + &
                aterm*particles_local(VELZ_PART_PROP,i)
           particles_local(OACZ_PART_PROP,i) = particles_local(ACCZ_PART_PROP,i)

           particles_local(ipdtold,i) = dt

        end do

        ! add dt to the current time within subcycling loop
        t_sub = t_sub + dt
        nsubcycles = nsubcycles + 1

        if(debug)then
           print*, dr_globalMe, "nsubcycles=", nsubcycles
        endif

     end do

  end if

  if (dr_globalMe .eq. MASTER_PE) print*, "pt_sinkAdvanceParticles: #subcycles:", nsubcycles

  ! Added by JFG
  if (sim_fixedPartTag .ne. 0) then
      call pt_sinkGatherGlobal()
      do i = 1, localnpf
          if (idnint(particles_global(TAG_PART_PROP,i)) .eq. sim_fixedPartTag) then
              fixedi = i
              exit
          else
              cycle
          endif
          call Driver_abortFlash('Error: Unable to find fixed particle tag [2]')
      enddo
      fixed_vec(1) = particles_global(POSX_PART_PROP,fixedi) - fixed_vec(1)
      fixed_vec(2) = particles_global(POSY_PART_PROP,fixedi) - fixed_vec(2)
      fixed_vec(3) = particles_global(POSZ_PART_PROP,fixedi) - fixed_vec(3)
      fixed_vec(4) = particles_global(VELX_PART_PROP,fixedi) - fixed_vec(4)
      fixed_vec(5) = particles_global(VELY_PART_PROP,fixedi) - fixed_vec(5)
      fixed_vec(6) = particles_global(VELZ_PART_PROP,fixedi) - fixed_vec(6)
      print *, 'fixed_vec', fixed_vec
      do i = 1, localnp
          particles_local(POSX_PART_PROP,i) = particles_local(POSX_PART_PROP,i) - fixed_vec(1)
          particles_local(POSY_PART_PROP,i) = particles_local(POSY_PART_PROP,i) - fixed_vec(2)
          particles_local(POSZ_PART_PROP,i) = particles_local(POSZ_PART_PROP,i) - fixed_vec(3)
          particles_local(VELX_PART_PROP,i) = particles_local(VELX_PART_PROP,i) - fixed_vec(4)
          particles_local(VELY_PART_PROP,i) = particles_local(VELY_PART_PROP,i) - fixed_vec(5)
          particles_local(VELZ_PART_PROP,i) = particles_local(VELZ_PART_PROP,i) - fixed_vec(6)
      enddo
      call pt_sinkGatherGlobal()
  endif
  ! End JFG

  ! moved from Particles_advance
  call Grid_moveParticles(particles_local,pt_sinkParticleProps,&
    pt_maxSinksPerProc,localnp, pt_indexList, pt_indexCount, .false.)
  call Grid_sortParticles(particles_local,pt_sinkParticleProps,localnp, &
    NPART_TYPES,pt_maxSinksPerProc,particlesPerBlk,BLK_PART_PROP)
  
  return

end subroutine Particles_sinkAdvanceParticles

!===============================================================================

!               Functions for the second-order variable-timestep leapfrog
!               velocity update.

function pt_sinkBNCOEFF (DT, DT_OLD, ALPHA, DOTALPHA)
  !       This function is the "B_n" function from
  !       equation 29 in SecondOrder.tex

  implicit none

  real :: pt_sinkBNCOEFF
  real, INTENT(in)  :: DT, DT_OLD, ALPHA, DOTALPHA
  real, parameter   :: onetwelfth = 1./12.
  real, parameter   :: onesixth   = 1./6.

  pt_sinkBNCOEFF = 1. - 0.5*ALPHA*DT + DT**2*(ALPHA**2 - DOTALPHA)*onesixth
  pt_sinkBNCOEFF = pt_sinkBNCOEFF * (1. - 0.5*ALPHA*DT_OLD + &
       DT_OLD**2*(ALPHA**2+2*DOTALPHA)*onetwelfth)

  return
end function pt_sinkBNCOEFF

!===============================================================================

function pt_sinkCNCOEFF (DT, DT_OLD, ALPHA, DOTALPHA)
  !       This function is the "C_n" function from
  !       equation 31 in SecondOrder.tex

  implicit none

  real :: pt_sinkCNCOEFF
  real, INTENT(in)   :: DT, DT_OLD, ALPHA, DOTALPHA
  real, parameter    :: onesixth   = 1./6.
  real, parameter    :: onethird   = 1./3.

  pt_sinkCNCOEFF = 0.5*DT + DT_OLD*onethird + DT**2/(6*DT_OLD) - &
       ALPHA*DT*(DT+DT_OLD)*onesixth

  return

end function pt_sinkCNCOEFF

!===============================================================================

function pt_sinkDNCOEFF (DT, DT_OLD, ALPHA, DOTALPHA)
  !       This function is the "D_n" function from
  !       equation 33 in SecondOrder.tex

  implicit none

  real :: pt_sinkDNCOEFF
  real, INTENT(in)  :: DT, DT_OLD, ALPHA, DOTALPHA
  real, parameter   :: onetwelfth = 1./12.

  pt_sinkDNCOEFF = (DT_OLD**2 - DT**2)/(6*DT_OLD) - ALPHA*DT_OLD*(DT+DT_OLD)*onetwelfth

  return
end function pt_sinkDNCOEFF

!===============================================================================

subroutine pt_sinkCsm_integrateFriedmann(a,dt,anew)

  implicit none

  real, intent(IN)  :: a
  real, intent(IN)  :: dt
  real, intent(OUT) :: anew

  real :: dt2, dt6, dtn
  real :: a1
  real :: dadt,dadt1,dadt2

  interface
    subroutine pt_sinkCsm_friedmannDeriv(avar, dadtvar)
      real, intent(IN) :: avar
      real, intent(INOUT) :: dadtvar
    end subroutine
  end interface

  ! Use a Runge-Kutta method

  dt2 = 0.5*dt
  dt6 = dt/6.
  dtn = dt + dt2

  call pt_sinkCsm_friedmannDeriv(a,dadt)

  a1 = a + dt2*dadt

  call pt_sinkCsm_friedmannDeriv(a1,dadt1)

  a1 = a + dt2*dadt1

  call pt_sinkCsm_friedmannDeriv(a1,dadt2)

  a1 = a + dt*dadt2
  dadt2 = dadt1 + dadt2

  call pt_sinkCsm_friedmannDeriv(a1,dadt1)

  anew = a + dt6*(dadt + dadt1 + 2.*dadt2)

  return

end subroutine pt_sinkCsm_integrateFriedmann

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pt_sinkCsm_friedmannDeriv(avar, dadtvar)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  real, intent(IN) :: avar
  real, intent(INOUT) :: dadtvar

  real, save :: hubble, omega, lambda, curv

  logical, save :: first_call = .true.

  if(first_call) then
     call RuntimeParameters_get("HubbleConstant", hubble)
     call RuntimeParameters_get("OmegaMatter", omega)
     call RuntimeParameters_get("CosmologicalConstant", lambda)
     curv = omega + lambda - 1.0

     first_call = .false.
  end if

  dadtvar = hubble*sqrt(omega/avar - curv + lambda*avar**2)

  !===============================================================================

  return

end subroutine pt_sinkCsm_friedmannDeriv
