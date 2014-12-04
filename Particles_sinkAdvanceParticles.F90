!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkAdvanceParticles
!!
!! NAME
!!
!!  Particles_sinkAdvanceParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAdvanceParticles(real, INTENT(in) :: dr_dt)
!!
!! DESCRIPTION
!!
!!  Updates sink particle postions and velocities based on gravitational accelerations,
!!  using Leapfrog or Euler, depending on the user's choice There is also a special
!!  implementation of Leapfrog for cosmological simulations.
!!  This routine performs subcycling on sink-sink interactions, in case of very close
!!  encounters and highly eccentric orbits.
!!
!! ARGUMENTS
!!
!!   dr_dt - the current time step
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   removed call to AccelGasOnSinks (now called in Gravity_PotentialListOfBlocks)
!!   (CTSS, CF 2013)
!!
!!***

subroutine Particles_sinkAdvanceParticles(dr_dt)

  ! Advance sink particles through timestep dr_dt

  use Particles_sinkData
  use pt_sinkInterface, only: pt_sinkAccelSinksOnSinks, pt_sinkGetSubCycleTimeStep, &
                              pt_sinkGatherGlobal
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_globalMe
  use Cosmology_interface, ONLY :  Cosmology_getParams, Cosmology_getRedshift, Cosmology_getOldRedshift
  use Driver_interface, ONLY : Driver_getSimTime, Driver_abortFlash
  ! Added by JFG
  use pt_sinkInterface, only: pt_sinkGatherGlobal
  use Simulation_data, only: sim_fixedPartTag, sim_tRelax, sim_comAccel
  use Driver_data, ONLY : dr_simTime
  ! End JFG
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  implicit none

  real, INTENT(in) :: dr_dt

  real, dimension(maxsinks) :: ax_gas, ay_gas, az_gas

  integer             :: i, nsubcycles, np
  logical, save       :: first_call = .true.
  real                :: dt, dt_global, t_sub, wterm
  real, save          :: dt_subOld
  real                :: local_min_radius, local_max_accel
  logical             :: end_subcycling
  integer, save       :: integrator
  character(len=80)   :: sink_integrator
  real                :: alpha, dotalpha, aterm, bterm, woldterm
  real, save          :: Hubble, OmegaM, OmegaCurv, OmegaL, OmegaB
  real                :: redshift_old, sOld, a0, simTime

  ! this should be an interface
  real                :: pt_sinkBNCOEFF, pt_sinkCNCOEFF, pt_sinkDNCOEFF

  integer, parameter :: gather_nprops = 13
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: ACCX_PART_PROP, ACCY_PART_PROP, ACCZ_PART_PROP, &
                  POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP, &
                  VELX_PART_PROP, VELY_PART_PROP, VELZ_PART_PROP, &
                  OACX_PART_PROP, OACY_PART_PROP, OACZ_PART_PROP, &
                  DTOLD_PART_PROP /)

  ! Added by JFG
  real, dimension(6)          :: fixed_vec
  real, dimension(3,maxsinks) :: tempAccel
  integer                     :: fixedi
  ! End JFG

  interface
    subroutine pt_sinkCsm_integrateFriedmann(a,dt,anew)
      real, intent(IN)  :: a
      real, intent(IN)  :: dt
      real, intent(OUT) :: anew
    end subroutine
  end interface

  if (dr_simTime .lt. sim_tRelax) return

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

     call Cosmology_getParams(Hubble, OmegaM, OmegaB, OmegaL)
     OmegaCurv = OmegaM + OmegaL - 1.

     first_call = .false.
  end if

  if (localnpf .eq. 0) return

  ! Decide whether to update global or local list
  if (sink_AdvanceSerialComputation) then
     ! this requires that all CPUs update all particles
     ! in the global list (serial computation), but avoids communication during subcycling
     np = localnpf
     ! we only have to call pt_sinkGatherGlobal once here,
     ! instead of many times in pt_sinkAccelSinksOnSinks()
     call pt_sinkGatherGlobal(gather_propinds, gather_nprops)
     ! store acceleration from GasOnSinksAccel (now called in Gravity_PotentialListOfBlocks)
     ax_gas(1:np) = particles_global(ACCX_PART_PROP, 1:np)
     ay_gas(1:np) = particles_global(ACCY_PART_PROP, 1:np)
     az_gas(1:np) = particles_global(ACCZ_PART_PROP, 1:np)
  else
     ! this parallelizes, but requires communicating the particle
     ! positions of all particles in pt_sinkAccelSinksOnSinks()
     np = localnp
     ! store acceleration from GasOnSinksAccel (now called in Gravity_PotentialListOfBlocks)
     ax_gas(1:np) = particles_local(ACCX_PART_PROP, 1:np)
     ay_gas(1:np) = particles_local(ACCY_PART_PROP, 1:np)
     az_gas(1:np) = particles_local(ACCZ_PART_PROP, 1:np)
  endif

  ! Added by JFG
  ! NOTE: THIS ONLY CURRENTLY WORKS IF ADVANCE SERIAL IS FALSE!
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
      !fixed_vec(1) = particles_global(POSX_PART_PROP,fixedi)
      !fixed_vec(2) = particles_global(POSY_PART_PROP,fixedi)
      !fixed_vec(3) = particles_global(POSZ_PART_PROP,fixedi)
      !fixed_vec(4) = particles_global(VELX_PART_PROP,fixedi)
      !fixed_vec(5) = particles_global(VELY_PART_PROP,fixedi)
      !fixed_vec(6) = particles_global(VELZ_PART_PROP,fixedi)
  endif
  ! End JFG

  dt_global = dr_dt
  t_sub = 0.0
  nsubcycles = 0
  end_subcycling = .false.

  if (integrator .eq. 1) then   !!!!!!! Leapfrog

     do while (.not. end_subcycling)

        ! copy gas acceleration into list
        if (sink_AdvanceSerialComputation) then
          particles_global(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_global(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_global(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        else
          particles_local(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_local(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_local(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        endif

        ! compute sink-sink gravity
        call pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

        call pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)

        ! JFG
        call pt_sinkGatherGlobal()
        particles_local(ACCX_PART_PROP, 1:np) = particles_local(ACCX_PART_PROP, 1:np) - particles_global(ACCX_PART_PROP, fixedi)
        particles_local(ACCY_PART_PROP, 1:np) = particles_local(ACCY_PART_PROP, 1:np) - particles_global(ACCY_PART_PROP, fixedi)
        particles_local(ACCZ_PART_PROP, 1:np) = particles_local(ACCZ_PART_PROP, 1:np) - particles_global(ACCZ_PART_PROP, fixedi)
        ! End JFG

        ! check whether we have reached the global (hydro) timestep for exiting
        if (t_sub + dt .ge. dt_global) then
           dt = dt_global - t_sub
           end_subcycling = .true.
        end if

        ! First step of leapfrog:
        wterm = 0.5 * dt

        if (sink_AdvanceSerialComputation) then

          do i = 1, np
            particles_global(VELX_PART_PROP, i) =  & 
                particles_global(VELX_PART_PROP, i) + & 
                wterm * particles_global(ACCX_PART_PROP, i)
            particles_global(POSX_PART_PROP, i) = & 
                particles_global(POSX_PART_PROP, i) + & 
                dt*particles_global(VELX_PART_PROP, i)

            particles_global(VELY_PART_PROP, i) =  & 
                particles_global(VELY_PART_PROP, i) + & 
                wterm * particles_global(ACCY_PART_PROP, i)
            particles_global(POSY_PART_PROP, i) = & 
                particles_global(POSY_PART_PROP, i) + & 
                dt*particles_global(VELY_PART_PROP, i)

            particles_global(VELZ_PART_PROP, i) =  & 
                particles_global(VELZ_PART_PROP, i) + & 
                wterm * particles_global(ACCZ_PART_PROP, i)
            particles_global(POSZ_PART_PROP, i) = & 
                particles_global(POSZ_PART_PROP, i) + & 
                dt*particles_global(VELZ_PART_PROP, i)
          enddo

        else ! parallel

          do i = 1, np
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

        endif ! sink_AdvanceSerialComputation

        ! reset list accelerations to gas contribution only
        if (sink_AdvanceSerialComputation) then
          particles_global(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_global(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_global(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        else
          particles_local(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_local(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_local(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        endif

        ! Add sink-sink acceleration at new position of particle
        call pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

        ! JFG
        call pt_sinkGatherGlobal()
        
        if (sim_fixedPartTag .ne. 0 .and. t_sub + dt .ge. dt_global) then
            sim_comAccel(1) = particles_global(ACCX_PART_PROP,fixedi)
            sim_comAccel(2) = particles_global(ACCY_PART_PROP,fixedi)
            sim_comAccel(3) = particles_global(ACCZ_PART_PROP,fixedi)
        endif

        particles_local(ACCX_PART_PROP, 1:np) = particles_local(ACCX_PART_PROP, 1:np) - particles_global(ACCX_PART_PROP, fixedi)
        particles_local(ACCY_PART_PROP, 1:np) = particles_local(ACCY_PART_PROP, 1:np) - particles_global(ACCY_PART_PROP, fixedi)
        particles_local(ACCZ_PART_PROP, 1:np) = particles_local(ACCZ_PART_PROP, 1:np) - particles_global(ACCZ_PART_PROP, fixedi)
        ! End JFG

        ! Second step of leapfrog
        if (sink_AdvanceSerialComputation) then

          do i = 1, np

            particles_global(VELX_PART_PROP, i) = & 
                particles_global(VELX_PART_PROP, i) + & 
                wterm*particles_global(ACCX_PART_PROP, i)

            particles_global(VELY_PART_PROP, i) = & 
                particles_global(VELY_PART_PROP, i) + & 
                wterm*particles_global(ACCY_PART_PROP, i)

            particles_global(VELZ_PART_PROP, i) = & 
                particles_global(VELZ_PART_PROP, i) + & 
                wterm*particles_global(ACCZ_PART_PROP, i)

          enddo

        else ! parallel

          do i = 1, np

            particles_local(VELX_PART_PROP, i) = & 
                particles_local(VELX_PART_PROP, i) + & 
                wterm*particles_local(ACCX_PART_PROP, i)

            particles_local(VELY_PART_PROP, i) = & 
                particles_local(VELY_PART_PROP, i) + & 
                wterm*particles_local(ACCY_PART_PROP, i)

            particles_local(VELZ_PART_PROP, i) = & 
                particles_local(VELZ_PART_PROP, i) + & 
                wterm*particles_local(ACCZ_PART_PROP, i)

          enddo

        endif ! sink_AdvanceSerialComputation

        ! add dt to the current time within subcycling loop
        t_sub = t_sub + dt
        nsubcycles = nsubcycles + 1

     end do

  end if    ! leapfrog integrator

  if (integrator .eq. 2) then       ! euler method integrator.

     do while(.not. end_subcycling)

        ! copy gas acceleration into list
        if (sink_AdvanceSerialComputation) then
          particles_global(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_global(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_global(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        else
          particles_local(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_local(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_local(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        endif

        call pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

        call pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)

        if (t_sub + dt .ge. dt_global) then
           dt = dt_global - t_sub
           end_subcycling = .true.
        endif

        if (sink_AdvanceSerialComputation) then

          do i = 1, np

            particles_global(POSX_PART_PROP,i) = particles_global(POSX_PART_PROP,i) + & 
                dt * particles_global(VELX_PART_PROP,i)
            particles_global(VELX_PART_PROP,i) = particles_global(VELX_PART_PROP,i) + & 
                dt * particles_global(ACCX_PART_PROP,i)

            particles_global(POSY_PART_PROP,i) = particles_global(POSY_PART_PROP,i) + & 
                dt * particles_global(VELY_PART_PROP,i)
            particles_global(VELY_PART_PROP,i) = particles_global(VELY_PART_PROP,i) + & 
                dt * particles_global(ACCY_PART_PROP,i)

            particles_global(POSZ_PART_PROP,i) = particles_global(POSZ_PART_PROP,i) + & 
                dt * particles_global(VELZ_PART_PROP,i)
            particles_global(VELZ_PART_PROP,i) = particles_global(VELZ_PART_PROP,i) + & 
                dt * particles_global(ACCZ_PART_PROP,i)

          enddo

        else ! parallel

          do i = 1, np

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

        endif ! sink_AdvanceSerialComputation

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

        ! copy gas acceleration into list
        if (sink_AdvanceSerialComputation) then
          particles_global(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_global(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_global(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        else
          particles_local(ACCX_PART_PROP, 1:np) = ax_gas(1:np)
          particles_local(ACCY_PART_PROP, 1:np) = ay_gas(1:np)
          particles_local(ACCZ_PART_PROP, 1:np) = az_gas(1:np)
        endif

        ! compute sink-sink gravity
        call pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

        call pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)

        ! JFG
        tempAccel(1,1:np) = particles_local(ACCX_PART_PROP,1:np)
        tempAccel(2,1:np) = particles_local(ACCY_PART_PROP,1:np)
        tempAccel(3,1:np) = particles_local(ACCZ_PART_PROP,1:np)

        particles_local(ACCX_PART_PROP,1:np) = particles_local(ACCX_PART_PROP,1:np) - ax_gas(1:np)
        particles_local(ACCY_PART_PROP,1:np) = particles_local(ACCY_PART_PROP,1:np) - ay_gas(1:np)
        particles_local(ACCZ_PART_PROP,1:np) = particles_local(ACCZ_PART_PROP,1:np) - az_gas(1:np)

        call pt_sinkGatherGlobal()

        sim_comAccel(1) = particles_global(ACCX_PART_PROP,fixedi)
        sim_comAccel(2) = particles_global(ACCY_PART_PROP,fixedi)
        sim_comAccel(3) = particles_global(ACCZ_PART_PROP,fixedi)

        particles_local(ACCX_PART_PROP,1:np) = tempAccel(1,1:np)
        particles_local(ACCY_PART_PROP,1:np) = tempAccel(2,1:np)
        particles_local(ACCZ_PART_PROP,1:np) = tempAccel(3,1:np)

        if (sim_fixedPartTag .ne. 0 .and. t_sub + dt .ge. dt_global) then
            if (dr_globalMe .eq. MASTER_PE) then
                print *, 'sim_comAccel', sim_comAccel
                print *, 'a_gas', ax_gas(fixedi), ay_gas(fixedi), az_gas(fixedi)
            !    print *, 'particle', particles_global(:,fixedi)
            endif
        endif

        particles_local(ACCX_PART_PROP, 1:np) = particles_local(ACCX_PART_PROP, 1:np) - sim_comAccel(1)
        particles_local(ACCY_PART_PROP, 1:np) = particles_local(ACCY_PART_PROP, 1:np) - sim_comAccel(2)
        particles_local(ACCZ_PART_PROP, 1:np) = particles_local(ACCZ_PART_PROP, 1:np) - sim_comAccel(3)
        ! End JFG

        ! check whether we have reached the global (hydro) timestep for exiting
        if (t_sub + dt .ge. dt_global) then
           dt = dt_global - t_sub
           end_subcycling = .true.
        end if

        alpha    = 2.0*Hubble/sOld * sqrt(OmegaM/sOld-OmegaCurv+OmegaL*sOld**2)
        dotalpha = 2.0*Hubble**2 * (-3.0*OmegaM + 2.0*OmegaCurv*sOld)/sOld**3

        call pt_sinkCsm_integrateFriedmann(a0, dt, sOld)
        ! sOld is now the scale factor at start of subcycle timestep
        a0 = sOld
        simTime = simTime + dt

        ! Advances particle positions fully over the sink particle
        ! subcycle timestep
        ! Velocities are defined in center of timestep

        if (sink_AdvanceSerialComputation) then

         do i = 1, np

           dt_subOld = particles_global(DTOLD_PART_PROP,i)

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
              wterm      = pt_sinkCNCOEFF(dt, dt_subOld, alpha)
              woldterm   = pt_sinkDNCOEFF(dt, dt_subOld, alpha)
           end if

           particles_global(VELX_PART_PROP,i) = bterm*particles_global(VELX_PART_PROP,i) + &
                wterm*particles_global(ACCX_PART_PROP,i) + &
                woldterm*particles_global(OACX_PART_PROP,i)
           particles_global(POSX_PART_PROP,i) = particles_global(POSX_PART_PROP,i) + &
                aterm*particles_global(VELX_PART_PROP,i)
           particles_global(OACX_PART_PROP,i) = particles_global(ACCX_PART_PROP,i)

           particles_global(VELY_PART_PROP,i)= bterm*particles_global(VELY_PART_PROP,i) + &
                wterm*particles_global(ACCY_PART_PROP,i) + &
                woldterm*particles_global(OACY_PART_PROP,i)
           particles_global(POSY_PART_PROP,i)=  particles_global(POSY_PART_PROP,i) + &
                aterm*particles_global(VELY_PART_PROP,i)
           particles_global(OACY_PART_PROP,i) = particles_global(ACCY_PART_PROP,i)

           particles_global(VELZ_PART_PROP,i)= bterm*particles_global(VELZ_PART_PROP,i) + &
                wterm*particles_global(ACCZ_PART_PROP,i) + &
                woldterm*particles_global(OACZ_PART_PROP,i)
           particles_global(POSZ_PART_PROP,i) =  particles_global(POSZ_PART_PROP,i) + &
                aterm*particles_global(VELZ_PART_PROP,i)
           particles_global(OACZ_PART_PROP,i) = particles_global(ACCZ_PART_PROP,i)

           particles_global(DTOLD_PART_PROP,i) = dt

         end do

        else ! parallel

         do i = 1, np

           dt_subOld = particles_local(DTOLD_PART_PROP,i)

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
              wterm      = pt_sinkCNCOEFF(dt, dt_subOld, alpha)
              woldterm   = pt_sinkDNCOEFF(dt, dt_subOld, alpha)
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

           particles_local(DTOLD_PART_PROP,i) = dt

         end do

        endif ! sink_AdvanceSerialComputation

        ! add dt to the current time within subcycling loop
        t_sub = t_sub + dt
        nsubcycles = nsubcycles + 1

     end do

  end if ! leapfrog cosmo

  ! Added by JFG
  !if (sim_fixedPartTag .ne. 0) then
  !    call pt_sinkGatherGlobal()
  !    do i = 1, localnpf
  !        if (idnint(particles_global(TAG_PART_PROP,i)) .eq. sim_fixedPartTag) then
  !            fixedi = i
  !            exit
  !        else
  !            cycle
  !        endif
  !        call Driver_abortFlash('Error: Unable to find fixed particle tag [2]')
  !    enddo
  !    !fixed_vec(1) = particles_global(POSX_PART_PROP,fixedi) - fixed_vec(1)
  !    !fixed_vec(2) = particles_global(POSY_PART_PROP,fixedi) - fixed_vec(2)
  !    !fixed_vec(3) = particles_global(POSZ_PART_PROP,fixedi) - fixed_vec(3)
  !    !fixed_vec(4) = particles_global(VELX_PART_PROP,fixedi) - fixed_vec(4)
  !    !fixed_vec(5) = particles_global(VELY_PART_PROP,fixedi) - fixed_vec(5)
  !    !fixed_vec(6) = particles_global(VELZ_PART_PROP,fixedi) - fixed_vec(6)

  !    !if (sink_AdvanceSerialComputation) then
  !    !    particles_global(POSX_PART_PROP,1:localnp) = particles_global(POSX_PART_PROP,1:localnp) - fixed_vec(1)
  !    !    particles_global(POSY_PART_PROP,1:localnp) = particles_global(POSY_PART_PROP,1:localnp) - fixed_vec(2)
  !    !    particles_global(POSZ_PART_PROP,1:localnp) = particles_global(POSZ_PART_PROP,1:localnp) - fixed_vec(3)
  !    !    particles_global(VELX_PART_PROP,1:localnp) = particles_global(VELX_PART_PROP,1:localnp) - fixed_vec(4)
  !    !    particles_global(VELY_PART_PROP,1:localnp) = particles_global(VELY_PART_PROP,1:localnp) - fixed_vec(5)
  !    !    particles_global(VELZ_PART_PROP,1:localnp) = particles_global(VELZ_PART_PROP,1:localnp) - fixed_vec(6)
  !    !else
  !    !    particles_local(POSX_PART_PROP,1:localnp) = particles_local(POSX_PART_PROP,1:localnp) - fixed_vec(1)
  !    !    particles_local(POSY_PART_PROP,1:localnp) = particles_local(POSY_PART_PROP,1:localnp) - fixed_vec(2)
  !    !    particles_local(POSZ_PART_PROP,1:localnp) = particles_local(POSZ_PART_PROP,1:localnp) - fixed_vec(3)
  !    !    particles_local(VELX_PART_PROP,1:localnp) = particles_local(VELX_PART_PROP,1:localnp) - fixed_vec(4)
  !    !    particles_local(VELY_PART_PROP,1:localnp) = particles_local(VELY_PART_PROP,1:localnp) - fixed_vec(5)
  !    !    particles_local(VELZ_PART_PROP,1:localnp) = particles_local(VELZ_PART_PROP,1:localnp) - fixed_vec(6)
  !    !endif
  !endif
  ! End JFG

  ! in case of serial computation, copy global list into local list, so we leave updated
  if (sink_AdvanceSerialComputation) then
     particles_local(POSX_PART_PROP,1:localnp) = particles_global(POSX_PART_PROP,1:localnp)
     particles_local(POSY_PART_PROP,1:localnp) = particles_global(POSY_PART_PROP,1:localnp)
     particles_local(POSZ_PART_PROP,1:localnp) = particles_global(POSZ_PART_PROP,1:localnp)
     particles_local(VELX_PART_PROP,1:localnp) = particles_global(VELX_PART_PROP,1:localnp)
     particles_local(VELY_PART_PROP,1:localnp) = particles_global(VELY_PART_PROP,1:localnp)
     particles_local(VELZ_PART_PROP,1:localnp) = particles_global(VELZ_PART_PROP,1:localnp)
     particles_local(ACCX_PART_PROP,1:localnp) = particles_global(ACCX_PART_PROP,1:localnp)
     particles_local(ACCY_PART_PROP,1:localnp) = particles_global(ACCY_PART_PROP,1:localnp)
     particles_local(ACCZ_PART_PROP,1:localnp) = particles_global(ACCZ_PART_PROP,1:localnp)
     if (integrator .eq. 3) then
       particles_local(OACX_PART_PROP,1:localnp) = particles_global(OACX_PART_PROP,1:localnp)
       particles_local(OACY_PART_PROP,1:localnp) = particles_global(OACY_PART_PROP,1:localnp)
       particles_local(OACZ_PART_PROP,1:localnp) = particles_global(OACZ_PART_PROP,1:localnp)
       particles_local(DTOLD_PART_PROP,1:localnp) = particles_global(DTOLD_PART_PROP,1:localnp)
     endif
  endif

  if (dr_globalMe .eq. MASTER_PE) print*, "Particles_sinkAdvanceParticles: #subcycles:", nsubcycles

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

function pt_sinkCNCOEFF (DT, DT_OLD, ALPHA)
  !       This function is the "C_n" function from
  !       equation 31 in SecondOrder.tex

  implicit none

  real :: pt_sinkCNCOEFF
  real, INTENT(in)   :: DT, DT_OLD, ALPHA
  real, parameter    :: onesixth   = 1./6.
  real, parameter    :: onethird   = 1./3.

  pt_sinkCNCOEFF = 0.5*DT + DT_OLD*onethird + DT**2/(6*DT_OLD) - &
       ALPHA*DT*(DT+DT_OLD)*onesixth

  return

end function pt_sinkCNCOEFF

!===============================================================================

function pt_sinkDNCOEFF (DT, DT_OLD, ALPHA)
  !       This function is the "D_n" function from
  !       equation 33 in SecondOrder.tex

  implicit none

  real :: pt_sinkDNCOEFF
  real, INTENT(in)  :: DT, DT_OLD, ALPHA
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

  call pt_sinkCsm_friedmannDeriv(a, dadt)

  a1 = a + dt2*dadt

  call pt_sinkCsm_friedmannDeriv(a1, dadt1)

  a1 = a + dt2*dadt1

  call pt_sinkCsm_friedmannDeriv(a1, dadt2)

  a1 = a + dt*dadt2
  dadt2 = dadt1 + dadt2

  call pt_sinkCsm_friedmannDeriv(a1, dadt1)

  anew = a + dt6*(dadt + dadt1 + 2.*dadt2)

  return

end subroutine pt_sinkCsm_integrateFriedmann

!===============================================================================

subroutine pt_sinkCsm_friedmannDeriv(avar, dadtvar)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Cosmology_interface, ONLY : Cosmology_getParams

  implicit none

  real, intent(IN) :: avar
  real, intent(INOUT) :: dadtvar

  real, save :: Hubble, OmegaM, OmegaB, OmegaL, curv

  logical, save :: first_call = .true.

  if(first_call) then
     call Cosmology_getParams(Hubble, OmegaM, OmegaB, OmegaL)
     curv = OmegaM + OmegaL - 1.0
     first_call = .false.
  end if

  dadtvar = Hubble*sqrt(OmegaM/avar - curv + OmegaL*avar**2)

  return

end subroutine pt_sinkCsm_friedmannDeriv

!===============================================================================
