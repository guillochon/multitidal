!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkMoveParticles
!!
!! NAME
!!
!!  Particles_sinkMoveParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkMoveParticles(logical(in) :: regrid)
!!
!! DESCRIPTION
!!
!!  Moves sink particles to the right block and processor depending on their
!!  physical positions in the computational domain (calls Grid_moveParticles).
!!
!! ARGUMENTS
!!
!!   regrid - logical flag indicating whether the grid changed or not
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

subroutine Particles_sinkMoveParticles(regrid)

  use Particles_sinkData, ONLY : particles_local, localnp,  & 
       sink_maxSinks, pt_sinkParticleProps
  use Particles_data, only: pt_indexList, pt_indexCount
  use Grid_interface, ONLY : Grid_moveParticles

  implicit none 
  
  logical, intent(in) :: regrid
  
  !call Grid_moveParticles(particles_local, pt_sinkParticleProps, sink_maxSinks, localnp, &
  !     pt_indexList, pt_indexCount, regrid)

end subroutine
