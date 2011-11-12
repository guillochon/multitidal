#if 0
MASTER_PE is the designated master processor in a run
MAX_PROCS value is needed for some allocation, keeping it in a 
  header file makes it easy to change the value if needed for 
  specific machines. The other are some general definitions for convenience
#endif
  
#define MASTER_PE 0
#define MAX_STRING_LENGTH 80
#define OUTPUT_PROP_LENGTH 24
#define REAL_FORMAT "(ES20.13)"
#define MAX_PROCS 132000
#define NO -10
#define NONEXISTENT -1
#define UNKNOWN -2
#define DUMP_IOFILE_NUM 9999

#define TAB_CHAR ACHAR(9)
#define REALSIZE 8

#define FINEST 78


#if 0
  This section defines real numeric constants that are so ubiquitous
  that it is worth having a shared definition of the values.
#endif

#define PI 3.1415926535897932384


#if 0
  This section defines the Grid geometries, which will be
  supported in future release. The first four are definitions
  for the whole grid, the last six define the geometries of 
  individual axes.
#endif

#define CARTESIAN 1
#define POLAR 2
#define CYLINDRICAL 3
#define SPHERICAL 4
#define XYZ 0
#define RAD_CYL 1
#define RAD_SPH 2
#define PHI_CYL 3
#define THETA 4
#define PHI_SPH 5 


#if 0
  This section defines the boundary conditions. Not all
  have implementations in the current release.
  The integer values must lie in the range -50..-20.
  PARAMESH_PHYSICAL_BOUNDARY should only be used when testing
  for the presence of a boundary, as in
     if (block_neighbor <= PARAMESH_PHYSICAL_BOUNDARY) then ...
  The last constant in the group is used in some places to indicate
  a surface that is not on a physical boundary. 
#endif

#define REFLECTING -31
#define OUTFLOW -32
#define PERIODIC -35
#define USER_DEFINED  -38
#define ISOLATED -33
#define HYDROSTATIC -34
#define DIRICHLET -36
#define PNEUMAN  -37
#define DIODE -39
#define GRIDBC_MG_EXTRAPOLATE -40
#define HYDROSTATIC_NVDIODE -41
#define HYDROSTATIC_NVREFL -42
#define HYDROSTATIC_NVOUT -43
#define HYDROSTATIC_NVZERO -44
#define HYDROSTATIC_F2 -45
#define HYDROSTATIC_F2_NVDIODE -46
#define HYDROSTATIC_F2_NVREFL -47
#define HYDROSTATIC_F2_NVOUT -48
#define MARSHAK -49
#define VACUUM -50

#define PARAMESH_PHYSICAL_BOUNDARY -20
#define GRIDBC_GIMME_WORK -19
#define NOT_BOUNDARY -10

#if 0
  The first three constants in this group are used to specify the edges of 
  cells or blocks, where for cell they have the obvious meaning. When used
  in connection with blocks, LEFT_EDGE 
  indicates guard cells to the left of the block, RIGHT_EDGE indicates
  guard cells to the right of the block and CENTER indicates the interior 
  of the block. The fourth and fifth constants in this group are used only 
  in connection with blocks, "WHOLE_VECTOR" is used 
  when the interior as well as the guard cells on
  both sides are referenced together, and ALLVARS is used when referencing 
  all the physical variables in the Grid data structures such as "unk"
#endif

#define LEFT_EDGE 1
#define CENTER 2
#define RIGHT_EDGE 3
#define WHOLE_VECTOR 4
#define NO_VEC 5
#define ALLVARS -1


#if 0
  This group has definition related to dimensions. The first three 
  are names of the axes. When referring to all the dimensions at once,
  ALLDIR is used, and MDIM defines the maximum number of dimensions
  supported in the code.
#endif

#define IAXIS 1
#define JAXIS 2
#define KAXIS 3
#define ALLDIR -1
#define MDIM 3


#if 0
  The next two constants are used in connection of integer block boundaries
  LOW refers to the lowest index and HIGH refers to the highest index of the 
  block. These can be used interchangeably for either the whole block 
  including guardcells or for the interior only.
#endif

#define LOW 1
#define HIGH 2


#if 0
  The next two constants are used as an argument to some subroutines that are
  called like brackets around some code, to distinguish between "preparation"
  and "cleanup" calls.
#endif

#define BEFORE 1
#define AFTER 2


#if 0
  The next two constants are used to indicate whether variables need to
  change from one form to another around guard cell exchange. If using
  conservative variables with some Grid implementations, it may be necessary to
  convert them to per-mass form before filling guard cells at fine-coarse boundaries
  to get the interpolation right. If NO_MAPPER is specified no conversions
  are done.
  These symbols are obsolete and not used in any code that is currently supported
  by the FLASH Center. They should not be used by new code.
#endif

#define NO_MAPPER 0
#define CONS_TO_PRIM 1

#if 0
  This group of constants defines options for getting a list of blocks. The 
  four refer to blocks that are on the physical boundary along the respective 
  axis. ACTIVE_BLKS refers to all blocks on which the solution is being 
  advanced. The last three are specific to Paramesh and indicate the position
  of the block in the tree.
#endif

#define IBDRY_BLKS 200
#define JBDRY_BLKS 201
#define KBDRY_BLKS 202
#define ANY_BDRY_BLKS 203
#define ACTIVE_BLKS 204
#define ALL_BLKS    205
#define LEAF 1
#define PARENT_BLK 2
#define ANCESTOR 3
#define REFINEMENT 321
#define TRAVERSED 254
#define TRAVERSED_AND_ACTIVE 278

#if 0
  These five constants are used in get/put data functions, The first 
  two are used to indicate whether to count the offset from the edge 
  that includes guardcells, or from the first interior cell. The last 
  indicate the plane for the Grid_get/putPlaneData functions.
#endif

#define INTERIOR 10
#define EXTERIOR 11
#define XYPLANE 55
#define XZPLANE 66
#define YZPLANE 77


#if 0
  This group refers to Grid data strucures, to store cell centered, face
  centered or scratch data for the physical domain. To indicate cell centered
  data we use CENTER which is defined in one of the earlier groups. WORK
  is specific to Paramesh and is used to manage a subset of cell centered
  variables. This group also has constants that can identify specific
  neighbor blocks in paramesh. MAX_GRID_DATA_STRUCT is the count of all
  all currently supported data structures
#endif

#if 0
  eventually, when we have all scratch structures in place MAX_GRID_DATA_STRUCT
  will be 9, for now it is 5
#endif

#define MAX_GRID_DATA_STRUCT 5
#define MAX_GRID_DATA_STRUCT_TMP 9
#define FACEX 3
#define FACEY 4
#define FACEZ 5
#define SCRATCH 1
#define SCRATCH_CTR 6
#define SCRATCH_FACEX 7
#define SCRATCH_FACEY 8
#define SCRATCH_FACEZ 9
#define SCRATCH_FACES 301
#define WORK 350
#define FACES 375
#define CENTER_FACES 380
#define CELL_VOLUME 382
#define CELL_FACEAREA 383

#if 0
  These constants define the grid variables that a given particle 
  property maps to through Simulation_mapParticlesVar()
#endif

#define PARTICLEMAP_UNK 1
#define PARTICLEMAP_SCRATCH 2
#define PARTICLEMAP_FACEX 3
#define PARTICLEMAP_FACEY 4
#define PARTICLEMAP_FACEZ 5
#define PARTICLEMAP_SCRATCH_CTR 6
#define PARTICLEMAP_SCRATCH_FACEX 7
#define PARTICLEMAP_SCRATCH_FACEY 8
#define PARTICLEMAP_SCRATCH_FACEZ 9

#if 0
  This constant defines the current maximum number of variables that 
  Paramesh will refine on.
#endif
#define MAXREFVARS 9


#if 0
  This group defines the supported EOS modes. MODE_RT, MODE_RP,
  and MODE_RE are for future use.
#endif
#define MODE_DENS_TEMP 101
#define MODE_DENS_EI 102
#define MODE_DENS_PRES 103
#define MODE_RT 105
#define MODE_RP 106
#define MODE_RE 107
#define MODE_EOS_NOP 55
#define MODE_EOS_WRAPPERONLY 67
#define MODE_DENS_ENTR 1207

#define MODE_DENS_TEMP_ION   30101
#define MODE_DENS_TEMP_ELE   30201
#define MODE_DENS_TEMP_RAD   30301
#define MODE_DENS_TEMP_COMP  31101
#define MODE_DENS_TEMP_ALL   31201
#define MODE_DENS_TEMP_EQUI  31301
#define MODE_DENS_TEMP_GATHER 31501

#define MODE_DENS_EI_ION     30102
#define MODE_DENS_EI_ELE     30202
#define MODE_DENS_EI_RAD     30302
#define MODE_DENS_EI_COMP    31102
#define MODE_DENS_EI_ALL     31202
#define MODE_DENS_EI_EQUI    31302
#define MODE_DENS_EI_SCATTER 31402
#define MODE_DENS_EI_GATHER  31502

#define MODE_DENS_EI_SELE_GATHER  32522
#define MODE_DENS_EI_SHOCKSELE_GATHER  33522

#define MODE_DENS_PRES_ION   30103
#define MODE_DENS_PRES_ELE   30203
#define MODE_DENS_PRES_RAD   30303
#define MODE_DENS_PRES_COMP  31103
#define MODE_DENS_PRES_ALL   31203

#if 0
These three constants define the sweep directions in the PPM algorithm.
They are also used in other directionally split solvers as well.
#endif

#define SWEEP_X 1
#define SWEEP_Y 2
#define SWEEP_Z 3
#define SWEEP_ALL 0
#define SWEEP_XYZ 1
#define SWEEP_ZYX 2
#define SWEEP_XZY 3
#define SWEEP_YZX 4
#define SWEEP_YXZ 5
#define SWEEP_ZXY 6


#if 0
  This group of constants is meant to be used with the specialized 
  refinement routines provided as reference with this release. The type
  of refinement is an argument in the routine, and these constants are
  the only valid values for that argument.
#endif

#define RECTANGLE 334
#define ELLIPSOID 335
#define THRESHOLD 336
#define INRADIUS 337
#define WITHRADIUS 338

#if 0
  These constants are used by the utilities that convert strings to 
  integer and vice-versa in the Simulation unit to map the components
  of the Grid data strucutures
#endif
#define MAPBLOCKSIZE 5000
#define MAPBLOCK_UNK 0
#define MAPBLOCK_FLUX 1
#define MAPBLOCK_PART 2
#define MAPBLOCK_SCRATCH 3
#define MAPBLOCK_FACES 4
#define MAPBLOCK_SCRATCH_CENTER 5
#define MAPBLOCK_SCRATCH_FACEX 6
#define MAPBLOCK_SCRATCH_FACEY 7
#define MAPBLOCK_SCRATCH_FACEZ 8

#if 0
  Symbols for Variable Types  
#endif
#define VARTYPE_ERROR 0
#define VARTYPE_GENERIC 1
#define VARTYPE_PER_VOLUME 2
#define VARTYPE_PER_MASS 3

#if 0
   The following constants clarify the specification of boundaries
   in domains that are not clean boxes, or in any way need physical
   boundaries somewhere inside the domain. The faces can also be combined
   with the last two constants to specify neighbors in paramesh.
#endif
#define ILO_FACE 1
#define JLO_FACE 3
#define KLO_FACE 5
#define IHI_FACE 2
#define JHI_FACE 4
#define KHI_FACE 6

#if 0
   These three number represent the fields in the surrblks datastructure of
   Paramesh. PROCNO and BLKNO is also the common way of uniquely identifying
   a block globally in Paramesh
#endif
#define BLKNO 1
#define PROCNO 2
#define TYPENO 3
#define ABSMAXNEGH 4


#if 0
   The next few constants are to facilitate the support for face centered
   variables in Uniform Grid. They are known only to Uniform Grid functions
#endif

#define NDATATYPES 4
#define CENTER_DATATYPE 1
#define FACEX_DATATYPE 2
#define FACEY_DATATYPE 3
#define FACEZ_DATATYPE 4


#if 0
  This group is constants is for use in applying boundary conditions.
  They define the indices for the array that stores the region of the
  block that has been extracted to apply boundary conditions.
#endif

#define BC_DIR 1
#define SECOND_DIR 2
#define THIRD_DIR 3
#define STRUCTSIZE 4
#define REGION_DIM 4


#if 0
  These constants are error codes for linked list get subroutines.
#endif
  
#define NORMAL 0
#define NOTFOUND -1
#define BADVALUE -2

#if 0
  IO_output has an argument that takes an integer quantity that denotes
  the type(s) of output files requested at that particular call.
#endif

#define CHECKPOINT_FILE_ONLY 1
#define PLOTFILE_ONLY 2
#define PARTICLE_FILE_ONLY 4
#define CHECKPOINT_AND_PLOTFILE (CHECKPOINT_FILE_ONLY + PLOTFILE_ONLY)
#define CHECKPOINT_AND_PARTICLEFILE (CHECKPOINT_FILE_ONLY + PARTICLE_FILE_ONLY)
#define PLOTFILE_AND_PARTICLEFILE (PLOTFILE_ONLY + PARTICLE_FILE_ONLY)
#define ALL_FILES (CHECKPOINT_FILE_ONLY + PLOTFILE_ONLY + PARTICLE_FILE_ONLY)

#ifndef PT_MAX_ATTRIBUTES
#define PT_MAX_ATTRIBUTES 10
#endif
#define PT_VAR 1
#define PT_MAP 2

#if 0
  These constants represent different modes for the flux limiter
#endif

#define FL_NONE     0
#define FL_HARMONIC 1
#define FL_MINMAX   2

#define GR_NORMAL 0
#define GR_QUIET_START 1
#define GR_PISTON 2
#define GR_HOLE 3

#if 0
These constants identify the type of communicators in use
They can either be the communicator that allows duplication of mesh
  in a simulation, or they can be directional as needed by the UG
The ones that allow duplication of the mesh have two dimensions
  one for the all processors that together have the copy of the mesh,
  and another that includes all processors that have identical rank in the
  first set of communicators.
#endif

#define GLOBAL_COMM 546
#define MESH_COMM 987
#define MESH_ACROSS_COMM 768
#define AXIS_COMM 854
#define NO_COMM 0

#if 0
  These constants represent different solvers and preconditioners for MG FLD
#endif

#define HYPRE_AMG 0
#define HYPRE_ILU 1
#define HYPRE_PCG 2
#define HYPRE_BICGSTAB 3
#define HYPRE_GMRES 4
#define HYPRE_SPLIT 5
#define HYPRE_NONE 6


#if 0
  These constants indicate the method for entering the radiation
  energy group boundaries for MGD. There is currently only one
  method of input supported - manual entry of the energy group 
  boundaries.
#endif
#define GRBD_MANUAL 0
