!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkAccelSinksOnSinks
!!
!! NAME
!!
!!  pt_sinkAccelSinksOnSinks
!!
!! SYNOPSIS
!!
!!  call pt_sinkAccelSinksOnSinks(real(out) :: local_min_radius,
!!                                real(out) :: local_max_accel)
!!
!! DESCRIPTION
!!
!!  Computes sinks <-> sinks gravitational accelerations.
!!  Compute the acceleration of the sink particles from other sink particles
!!  and add it to preexisting acceleration array. Send out local_min_radius
!!  and local_max_accel to constrain subcycle time step in Particles_sinkAdvanceParticles.
!!
!! ARGUMENTS
!!
!!   local_min_radius - minimum radius between sinks
!!
!!   local_max_accel - maximum acceleration between sinks
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged by Christoph Federrath, 2013
!!
!!***

subroutine pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

    use Particles_sinkData
    use pt_sinkInterface, ONLY: pt_sinkGatherGlobal, pt_sinkEwaldCorrection, &
                                pt_sinkCorrectForPeriodicBCs
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Driver_data, ONLY : dr_globalMe
    use Cosmology_interface, ONLY : Cosmology_getRedshift
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use Simulation_data, ONLY : sim_fixedPartTag, sim_gravityType

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

    real, intent(out) :: local_min_radius, local_max_accel

    character(len=80), save :: softening_type_sinks, grav_boundary_type
    logical, save :: first_call = .true.

    real, save    :: newton, softening_radius, c2
    real          :: softening_radius_comoving, slope, hinv, h2inv
    integer, save :: softeningtype
    real, save    :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz, local_min_radius_init

    integer       :: p, pf, np
    real          :: xp, yp, zp, dx, dy, dz, ax, ay, az, exc, eyc, ezc
    real          :: radius, masspf, q, kernelvalue, r3, paccx, paccy, paccz
    real          :: redshift, oneplusz3

    real          :: vxp, vyp, vzp, dvx, dvy, dvz, rsch, dvr, phi2
    integer       :: fixedi, i

    integer, parameter :: gather_nprops = 4
    integer, dimension(gather_nprops), save :: gather_propinds = &
      (/ integer :: POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP, MASS_PART_PROP /)
    ! ACCX_PART_PROP, ACCY_PART_PROP, ACCZ_PART_PROP do not have to be communicated,
    ! because they are set in Particles_sinkAdvanceParticles

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
    endif
    ! End JFG

    if (first_call) then

       call RuntimeParameters_get("sink_softening_radius", softening_radius)
       call RuntimeParameters_get("sink_softening_type_sinks", softening_type_sinks)
       select case (softening_type_sinks)
       case ("spline")
          softeningtype=1
       case ("linear")
          softeningtype=2
       case default
          softening_type_sinks = "spline"
          softeningtype = 2
          if(dr_globalMe .eq. MASTER_PE) print*, "pt_sinkAccelSinksOnSinks: invalid grav softening type specified"
       end select
       if(dr_globalMe .eq. MASTER_PE) print*, "pt_sinkAccelSinksOnSinks: grav softening type = ", &
            & trim(softening_type_sinks)

       call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

       if ((grav_boundary_type.ne."isolated").and.(grav_boundary_type.ne."periodic")) then
          call Driver_abortFlash("Sink particles can only be used with perioidic of isolated gravity type!")
       end if

       call RuntimeParameters_get("xmin", xmin)
       call RuntimeParameters_get("xmax", xmax)
       call RuntimeParameters_get("ymin", ymin)
       call RuntimeParameters_get("ymax", ymax)
       call RuntimeParameters_get("zmin", zmin)
       call RuntimeParameters_get("zmax", zmax)
       Lx = xmax-xmin
       Ly = ymax-ymin
       Lz = zmax-zmin
       local_min_radius_init = sqrt(Lx**2+Ly**2+Lz**2)

       call PhysicalConstants_get("Newton", newton)
       call PhysicalConstants_get("speed of light", c2)
       c2 = c2*c2

       first_call = .false.

    end if ! first call


    local_min_radius = local_min_radius_init
    local_max_accel = 0.0

    if (localnpf .eq. 0) return

    call Cosmology_getRedshift(redshift)
    softening_radius_comoving = softening_radius * (1.0 + redshift)

    hinv = 2.0 / softening_radius_comoving
    h2inv = hinv**2
    slope = 1.0 / softening_radius_comoving**3

    oneplusz3 = (1.0 + redshift) ** 3.0

    ! Decide whether to update global or local list
    if (sink_AdvanceSerialComputation) then
       ! this requires that all CPUs update all particles
       ! in the global list (serial computation)
       np = localnpf 
    else
       ! this parallelizes, but requires communicating the particle
       ! positions of all particles (in fact, we here call pt_sinkGatherGlobal,
       ! which additionally communicates unnecessary sink variables)
       np = localnp
       call pt_sinkGatherGlobal(gather_propinds, gather_nprops) ! this is time-consuming communication
    endif

    do p = 1, np

          paccx = 0.0
          paccy = 0.0
          paccz = 0.0

          xp = particles_global(POSX_PART_PROP,p)
          yp = particles_global(POSY_PART_PROP,p)
          zp = particles_global(POSZ_PART_PROP,p)

          vxp = particles_global(VELX_PART_PROP,p)
          vyp = particles_global(VELY_PART_PROP,p)
          vzp = particles_global(VELZ_PART_PROP,p)

          do pf = 1, localnpf ! this always goes over all sink particles

             if (pf .ne. p) then ! exclude self-interaction

                masspf = particles_global(MASS_PART_PROP,pf)
                dx = xp - particles_global(POSX_PART_PROP,pf)
                dy = yp - particles_global(POSY_PART_PROP,pf)
                dz = zp - particles_global(POSZ_PART_PROP,pf)
                dvx = vxp - particles_global(VELX_PART_PROP,pf)
                dvy = vyp - particles_global(VELY_PART_PROP,pf)
                dvz = vzp - particles_global(VELZ_PART_PROP,pf)

                if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(dx, dy, dz)

                radius = sqrt(dx**2 + dy**2 + dz**2)
                local_min_radius = min(local_min_radius, radius)

                if (radius .lt. softening_radius_comoving) then
                   if (softeningtype .eq. 1) then   ! spline softening
                      q = radius*hinv
                      if ((q.gt.1.e-5).and.(q.lt.1.0)) &
                         & kernelvalue = h2inv*(4.0/3.0*q-1.2*q**3+0.5*q**4)/radius
                      if ((q.ge.1.0)  .and.(q.lt.2.0)) &
                         & kernelvalue = h2inv * &
                         & (8.0/3.0*q-3.0*q**2+1.2*q**3-1.0/6.0*q**4-1.0/(15.0*q**2))/radius
                      ax = kernelvalue*dx
                      ay = kernelvalue*dy
                      az = kernelvalue*dz
                   end if
                   if (softeningtype.eq.2) then    ! linear kernel
                      ax = dx*slope
                      ay = dy*slope
                      az = dz*slope
                   end if
                else    
                   r3 = 1.0 / radius**3
                   if (sim_gravityType .eq. "newton") then ! Newtonian gravity of point mass
                       ax = dx*r3
                       ay = dy*r3
                       az = dz*r3
                   else
                       ! Schwarzschild metric (Gafton 2015)
                       rsch = 2.d0*newton*masspf/c2
                       dvr = (dx*dvx + dy*dvy + dz*dvz)/radius
                       phi2 = ((dx*dvy-dy*dvx)**2+(dx*dvz-dz*dvx)**2+(dz*dvy-dy*dvz)**2)/radius**4
                       ax = -(-newton*masspf*dx*r3*(1.d0-rsch/radius) + rsch*dvx*dvr/(radius*(radius - rsch)) + &
                            rsch*dx*dvr**2/(2.d0*(radius - rsch)*radius**2) - rsch*dx*phi2/radius)/(newton*masspf)
                       ay = -(-newton*masspf*dy*r3*(1.d0-rsch/radius) + rsch*dvy*dvr/(radius*(radius - rsch)) + &
                            rsch*dy*dvr**2/(2.d0*(radius - rsch)*radius**2) - rsch*dy*phi2/radius)/(newton*masspf)
                       az = -(-newton*masspf*dz*r3*(1.d0-rsch/radius) + rsch*dvz*dvr/(radius*(radius - rsch)) + &
                            rsch*dz*dvr**2/(2.d0*(radius - rsch)*radius**2) - rsch*dz*phi2/radius)/(newton*masspf)
                   endif
                end if

                if (grav_boundary_type .eq. "periodic") then
                   call pt_sinkEwaldCorrection(abs(dx), abs(dy), abs(dz), exc, eyc, ezc)
                   paccx = paccx + (ax - sign(exc,dx))*masspf
                   paccy = paccy + (ay - sign(eyc,dy))*masspf
                   paccz = paccz + (az - sign(ezc,dz))*masspf
                else
                   paccx = paccx + ax*masspf
                   paccy = paccy + ay*masspf
                   paccz = paccz + az*masspf
                endif

             end if ! pf .ne. p

          end do   ! pf

          if (sink_AdvanceSerialComputation) then

            particles_global(ACCX_PART_PROP,p) = particles_global(ACCX_PART_PROP,p) - newton*paccx*oneplusz3
            particles_global(ACCY_PART_PROP,p) = particles_global(ACCY_PART_PROP,p) - newton*paccy*oneplusz3
            particles_global(ACCZ_PART_PROP,p) = particles_global(ACCZ_PART_PROP,p) - newton*paccz*oneplusz3

            local_max_accel = max(local_max_accel, abs(particles_global(ACCX_PART_PROP,p)))
            local_max_accel = max(local_max_accel, abs(particles_global(ACCY_PART_PROP,p)))
            local_max_accel = max(local_max_accel, abs(particles_global(ACCZ_PART_PROP,p)))

          else

            particles_local(ACCX_PART_PROP,p) = particles_local(ACCX_PART_PROP,p) - newton*paccx*oneplusz3
            particles_local(ACCY_PART_PROP,p) = particles_local(ACCY_PART_PROP,p) - newton*paccy*oneplusz3
            particles_local(ACCZ_PART_PROP,p) = particles_local(ACCZ_PART_PROP,p) - newton*paccz*oneplusz3

            local_max_accel = max(local_max_accel, abs(particles_local(ACCX_PART_PROP,p)))
            local_max_accel = max(local_max_accel, abs(particles_local(ACCY_PART_PROP,p)))
            local_max_accel = max(local_max_accel, abs(particles_local(ACCZ_PART_PROP,p)))

          endif

    end do   ! local particles p

    return

end subroutine pt_sinkAccelSinksOnSinks
