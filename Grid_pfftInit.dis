!!****if* source/Grid/GridSolvers/Pfft/Grid_pfftInit
!!
!! NAME 
!!
!!   Grid_pfftInit
!!
!! SYNOPSIS
!!
!!   Grid_pfftInit(integer(IN)  :: ndim,
!!                 logical(IN)  :: needMap
!!                 integer(IN)  :: globalLen(MDIM),
!!                 integer(OUT) :: localArraySize(MDIM),
!!                 integer(IN),optional :: transformType(MDIM),
!!                 integer(IN),optional :: baseDatType(0:MDIM),
!!                 integer(IN),optional :: jProcs,
!!                 integer(IN),optional :: kProcs,
!!                 integer(IN),optional :: refinementLevel,
!!                 real(IN), optional   :: region_bndBox)
!!
!! DESCRIPTION 
!!  
!!  This is the initialization routine for using Pfft. If needMap
!!  is true,  routine creates a schedule of data transfers that would
!!  map AMR or UG grid to PFFT. This schedule is remembered internally
!!  by Pfft until the Grid_pfftFinalize routine is called. If needMap
!!  is false, it is assumed the data distribution is already compatible
!!  with Pfft requirements: That is along IAXIS, the entire row is 
!!  within one processor.
!!
!!  The routine also calls gr_pfftInitMetaData, which is responsible for
!!  calculating the trignometric tables, creating communicators necessary
!!  for distributed transposes in Pfft and allocates workspace
!!
!! ARGUMENTS
!!
!!  ndim          - dimensionality of the problem
!!  needMap       - should be true if the default shape is not compatible with
!!                  requirements of the input Pfft array. Only if
!!                  this argument is true, is the map determined.
!!  globalLen     -  the globalsize of the domain
!!  localArraySize       - after mapping to pfft_grid, the local size for the domain
!!  transformType -  type of transform along each dimension
!!                   if none is specified we assume realtocomplex
!!                   in first dimension and complextocomplex in the rest
!!  baseDatType -    basic data type (rela or complex) along each dimension,
!!                   after the transform for that direction (if any),
!!                   the 0 component specifies data type before IAXIS transform.
!!  jProcs,kProcs - if they are present, they decide the shape of the 
!!                  processor grid for pfft. If they are not present,
!!                  a routine that can automatically determine the shape is 
!!                  called.
!!  refinementLevel - The block refinement level at which we will create the map.
!!  region_bndBox  - If the map is from a subset of the physical domain in AMR to
!!                   Pfft grid, this argument contains the bounding box of the region
!!
!!***
subroutine Grid_pfftInit( ndim, needMap, globalLen, localArraySize, &
     transformType, baseDatType, jProcs, kProcs, refinementLevel, region_bndBox)

#include "Pfft.h"
#include "constants.h"
#include "Flash.h"

  use gr_pfftinterface, ONLY : gr_pfftGroupUsableProcesses, gr_pfftGridPointTable, &
       gr_pfftGenMap, gr_pfftInitMetadata, gr_pfftValidateSelectedLevel
  use gr_pfftData, ONLY : pfft_inLen, pfft_globalLen, pfft_transformType, &
       pfft_pclBaseDatType, &
       pfft_procGrid, pfft_ndim, pfft_needMap, pfft_comm, pfft_me, &
       pfft_commWithTopology, pfft_usableProc, pfft_myPE, pfft_numProcs, &
       pfft_mode, pfft_isSimpleNonMappedAMR1D, pfft_regionBndBox, pfft_inRegion
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_oneRefLev, gr_meshComm, gr_meshMe, gr_meshNumProcs
  use Grid_interface, ONLY : Grid_getMaxCommonRefinement
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_stampMessage
  use Timers_interface, ONLY : Timers_start, Timers_stop

#if defined(FLASH_GRID_UG)
  use Grid_data, ONLY : gr_axisMe, gr_axisNumProcs, gr_axisComm
#elif defined(FLASH_GRID_PARAMESH)
  use tree, ONLY : lrefine_max
#endif

  implicit none
  integer, intent(IN) :: ndim
  logical, intent(IN) :: needMap
  integer, dimension(MDIM), intent(IN) :: globalLen
  integer,dimension(MDIM),intent(OUT) ::  localArraySize
  integer,dimension(MDIM),optional,intent(IN) :: transformType
  integer,dimension(0:MDIM),optional,intent(IN) :: baseDatType
  integer,optional,intent(IN) :: jProcs, kProcs
  integer,optional,intent(IN) :: refinementLevel
  real, dimension(LOW:HIGH, MDIM), optional, intent(IN) :: region_bndBox

  logical, dimension(MDIM) :: periods
  real :: reductionFactor
  integer :: i, ierr, comm, maxFullyRefinedLevel
  logical :: reorder

  include "Flash_mpi.h"
  ! Sole change by JFG: Disable timers as they were causing issues
  !call Timers_start("PFFT_Init")

  localArraySize = 0
  pfft_myPE = gr_meshMe
  pfft_numProcs = gr_meshNumProcs
  pfft_ndim = ndim
  pfft_needMap = needMap
  pfft_globalLen(1:MDIM) = globalLen(1:MDIM)
  pfft_procGrid(IAXIS) = 1   !We need a pencil grid.
  pfft_usableProc = .true.   !Specifies whether we will use this processor.

  comm = gr_meshComm
  pfft_comm(1:MDIM) = (/gr_meshComm, MPI_COMM_NULL, MPI_COMM_NULL/)
  pfft_commWithTopology = MPI_COMM_NULL
  pfft_mode = PFFT_SINGLE_LEVEL !Typically we use only a single FLASH level.
  pfft_isSimpleNonMappedAMR1D = .false.  !This can be flipped to .true. in 1d AMR case.


  !Check that input options are consistent with this simulation.
  !Create a reduced PFFT communicator if too little work (PARAMESH case). 
  !---------------------------------------------------------------------

#if defined(FLASH_GRID_UG)

  !We have FFT problems when UG axis not divisible by no. procs in that axis.
  if (any (mod(globalLen(1:MDIM), gr_axisNumProcs(1:MDIM)) /= 0) ) then
     call Driver_abortFlash("[Grid_pfftInit]: Unusable UG grid")
  end if

  !Set needMap = .false. if the current grid is already a pencil grid.
  if (.not.needMap) then
     if (gr_axisNumProcs(IAXIS) == 1) then
        pfft_me = gr_axisMe
        pfft_procGrid = gr_axisNumProcs
        pfft_comm = gr_axisComm
     else
        call Driver_abortFlash("[Grid_pfftInit]: UG is not a pencil grid")
     end if
  end if


#elif defined(FLASH_GRID_PARAMESH)

  !We can use PFFT with PARAMESH in 1d and not need to create a map
  if (.not.needMap .and. NDIM==1) then
     pfft_isSimpleNonMappedAMR1D = .true.
  end if

  !We can select a value for gr_oneRefLev, or generate a value
  !based on the maximum level at which the grid is fully refined.
  call Grid_getMaxCommonRefinement(comm, maxFullyRefinedLevel)

  if (present(refinementLevel)) then
     gr_oneRefLev = refinementLevel

     if (refinementLevel > maxFullyRefinedLevel) then
        pfft_mode = PFFT_MIXED_LEVEL  !Only set when we can't use one level.
        call Logfile_stampMessage( & 
             "[Grid_pfftInit]: Using mode involving prolonging/restricting")

        !There is a very narrow regime where using mixed levels actually works. 
        !If the chosen level is too coarse we will increase the refinement 
        !level (stored in gr_oneRefLev) to satisfy these requirements.
        call gr_pfftValidateSelectedLevel(gr_oneRefLev)

        call Driver_abortFlash("MIXED_LEVEL not yet coded")
     end if
  else
     gr_oneRefLev = maxFullyRefinedLevel
  end if
  
  call Logfile_stamp( gr_oneRefLev, &
       "[Grid_pfftInit]: Generating pencil grid at level")


  !We need to adjust the number of grid points representing the 
  !global domain.  This is because we only need enough grid 
  !points to represent block data at refinement level gr_oneRefLev.
  reductionFactor = 2.0 ** (gr_oneRefLev - lrefine_max)
  pfft_globalLen(1:MDIM) = max(1, int(globalLen(1:MDIM) * reductionFactor))


  if (needMap .or. pfft_isSimpleNonMappedAMR1D) then
     !We use pfft_comm(IAXIS) to hold PFFT's global MPI communicator. Under 
     !normal circumstances it is "comm" communicator.  However, if there are 
     !more processors than blocks it will be a reduced processor communicator.
     call gr_pfftGroupUsableProcesses(pfft_myPE, pfft_numProcs, &
          comm, pfft_comm(IAXIS))

     !In the PARAMESH case we may have a situation in which there are more
     !processors than blocks.  This may happen if the value in gr_oneRefLev is
     !very small.  Therefore, in order to prevent distributing a block too 
     !finely (and avoid a potential crash) we eliminate superfluous processors 
     !(i.e. those processors that do not own any FLASH blocks).
     if (pfft_comm(IAXIS) == MPI_COMM_NULL) then
        pfft_myPE = NONEXISTENT; pfft_numProcs = NONEXISTENT
        pfft_usableProc = .false.
        return   ! ***** END *****  (superfluous processors).
     else
        !Update pfft_myPE and pfft_numProcs.  !Ignore numProcs variable now.
        call MPI_Comm_rank(pfft_comm(IAXIS), pfft_myPE, ierr)
        call MPI_Comm_size(pfft_comm(IAXIS), pfft_numProcs, ierr)
        pfft_usableProc = .true.
     end if

  else
     call Driver_abortFlash &
          ("[Grid_pfftInit]: Set neepMap to true when not using UG!")
  end if

#endif
  !---------------------------------------------------------------------


  !At this point we know the global length of the grid at a particular 
  !refinement level, and the number of processors we will use on this grid.
  !We now need to know the processor layout.  !Initialise pfft_procGrid.
  !---------------------------------------------------------------------
  if (needMap .or. pfft_isSimpleNonMappedAMR1D) then
     !We give the user an opportunity to select a grid.  If the requested
     !processor shape is not consistent with the number of available
     !processors we will abort.  Keep in mind that we will have reduced the
     !number of available processors (for PFFT) if blocks are spread too
     !finely over the processor grid!   To help avoid an abort due to PFFT
     !removing processors, we provide two definitions (PFFT_ALL_PROCS and
     !PFFT_ONE_PROC) which allow us to distribute work along one axis only.
     if ((ndim == 3) .and. (present(jProcs) .and. present(kProcs))) then
        if ((jProcs == PFFT_ALL_PROCS) .and. (kProcs == PFFT_ONE_PROC)) then
           pfft_procGrid(JAXIS) = pfft_numProcs; pfft_procGrid(KAXIS) = 1
        else if ((jProcs == PFFT_ONE_PROC) .and. (kProcs == PFFT_ALL_PROCS)) then
           pfft_procGrid(JAXIS) = 1; pfft_procGrid(KAXIS) = pfft_numProcs
        else
           if ((jProcs * kProcs) /= pfft_numProcs) then
              call Driver_abortFlash &
                   ("Requested (jProcs * kProcs) /= no. usable procs!")
           end if
           pfft_procGrid(JAXIS) = jProcs; pfft_procGrid(KAXIS) = kProcs
        end if
        if (pfft_myPE == 0) then
           print *, "[Grid_pfftInit]: User-defined PFFT processor grid:", &
                pfft_procGrid
        end if
     else
        call gr_pfftGetProcGrid(pfft_ndim, pfft_myPE, pfft_numProcs, &
             pfft_globalLen, pfft_procGrid)
        if (pfft_myPE == 0) then
           print *, "[Grid_pfftInit]: Generating PFFT processor grid:", &
                pfft_procGrid
        end if
     end if


     !We now have a PFFT grid (pfft_procGrid) and a communicator containing 
     !all processors in the grid (pfft_comm(IAXIS)).  Now create a new 
     !communicator containing a cartesian MPI topology which links a PFFT 
     !processor ID with its coordinates in PFFT grid.
     periods(1:MDIM) = .false.; reorder = .false.
     call MPI_Cart_create(pfft_comm(IAXIS), MDIM, pfft_procGrid(1), &
          periods(1), reorder, pfft_commWithTopology, ierr)

     !Determine my coordinates in PFFT grid and store in pfft_me. (Recall 
     !that pfft_myPE may have been updated in gr_pfftGroupUsableProcesses.)
     call MPI_Cart_coords(pfft_commWithTopology, pfft_myPE, MDIM, &
          pfft_me(1), ierr)

     !Create split communicators for the data transpose.
     call MPI_Comm_split(pfft_comm(IAXIS),pfft_me(KAXIS),pfft_me(JAXIS),&
          pfft_comm(JAXIS),ierr)
     call MPI_Comm_split(pfft_comm(IAXIS),pfft_me(JAXIS),pfft_me(KAXIS),&
          pfft_comm(KAXIS),ierr)
  end if
  !---------------------------------------------------------------------


  !At this point we have a pencil grid layout (either by default or generated).
  !Now find the local array size from global size.  The formula being used is 
  !(x+y-1)/y, which will round-off the quotient to the next higher integer if 
  !x is not an exact multiple of y. For example if x=17,18 or 19 and y = 4 then 
  !the answer is 5, but if x=16 then the answer is 4.
  !---------------------------------------------------------------------
  pfft_inLen(1:MDIM) = (pfft_globalLen(1:MDIM) + pfft_procGrid(1:MDIM) - 1) / &
       pfft_procGrid(1:MDIM)
  localArraySize(1:MDIM) = pfft_inLen(1:MDIM)
  !---------------------------------------------------------------------


  if(present(transformType)) then
     pfft_transformType(1:MDIM) = transformType(1:MDIM)
  else
     pfft_transformType(IAXIS) = PFFT_REAL2C
     pfft_transformType(JAXIS:KAXIS) = PFFT_COMPLEX
  end if

  if(present(baseDatType)) then
     pfft_pclBaseDatType(:) = baseDatType(:)
  else
     select case (pfft_transformType(1))
     case(PFFT_REAL2C,PFFT_REAL2C_STUFF,PFFT_REAL2C_EXTEND)
        pfft_pclBaseDatType(0) = PFFT_PCLDATA_REAL
        pfft_pclBaseDatType(1:MDIM) = PFFT_PCLDATA_COMPLEX
     case(PFFT_COMPLEX)
        pfft_pclBaseDatType(0:MDIM) = PFFT_PCLDATA_COMPLEX
     case default                  !includes PFFT_REAL
        pfft_pclBaseDatType(0:MDIM) = PFFT_PCLDATA_REAL
     end select
     do i=2,pfft_ndim
        select case (pfft_transformType(i))
        case(PFFT_REAL2C,PFFT_REAL2C_STUFF,PFFT_REAL2C_EXTEND)
           pfft_pclBaseDatType(i:MDIM) = PFFT_PCLDATA_COMPLEX
        case(PFFT_COMPLEX)
           pfft_pclBaseDatType(i:MDIM) = PFFT_PCLDATA_COMPLEX
        case default                  !includes PFFT_REAL
           pfft_pclBaseDatType(i:MDIM) = PFFT_PCLDATA_REAL
        end select
     end do
  end if
  
  pfft_regionBndBox=0.0
  pfft_inRegion=.false.
  if(present(region_bndBox)) then
     pfft_regionBndBox=region_bndBox
     pfft_inRegion=.true.
     gr_oneRefLev = NONEXISTENT
     if (present(refinementLevel))gr_oneRefLev = refinementLevel
  end if

  call gr_pfftInitMetaData(ndim)

  !Now figure out the mapping from the original grid to Pfft grid, if needed.
  !---------------------------------------------------------------------
  if (pfft_needMap) call gr_pfftGenMap()
  !---------------------------------------------------------------------

  !call Timers_stop("PFFT_Init")
end subroutine Grid_pfftInit
