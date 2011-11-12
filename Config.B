#       Configuration file for the Sedov explosion problem
#       (Sedov, L. I., 1959, Similarity and Dimensional Methods
#       in Mechanics [New York:  Academic])

DATAFILES *.dat
DATAFILES *.h

REQUIRES Driver
REQUIRES physics/Hydro
REQUIRES physics/Eos/EosMain/Helmholtz
REQUIRES Simulation/SimulationComposition
REQUIRES Multispecies
#REQUIRES physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19
#REQUIRES physics/sourceTerms/Cool

#REQUIRES physics/Gravity/GravityMain/PointMass
#REQUIRES Grid/GridSolvers/Multipole
REQUIRES physics/Gravity/GravityMain/Poisson/Multipole
#REQUIRES physics/Gravity/GravityMain/Poisson/Multigrid
IF withParticles
   REQUIRES Particles/ParticlesMain
   REQUIRES Particles/ParticlesMapping
   REQUIRES Particles/ParticlesInitialization/Lattice
   REQUESTS IO/IOMain/
   REQUESTS IO/IOParticles
ENDIF

# number of species (ionmax)
SPECIES H1
#SPECIES HE3
SPECIES HE4
SPECIES C12
SPECIES N14
SPECIES O16
#SPECIES NE20
#SPECIES MG24
#SPECIES SI28
#SPECIES S32
#SPECIES AR36
#SPECIES CA40
#SPECIES TI44
#SPECIES CR48
#SPECIES FE52
#SPECIES FE54   
#SPECIES NI56
#SPECIES NEUT
#SPECIES PROT

PARAMETER refine_val_cutoff_1 REAL    0.E0
PARAMETER refine_val_cutoff_2 REAL    0.E0
PARAMETER refine_val_cutoff_3 REAL    0.E0
PARAMETER refine_val_cutoff_4 REAL    0.E0
PARAMETER sim_periDist        REAL    1.E12
PARAMETER sim_periTime        REAL    200.
PARAMETER ptxpos              REAL    0.E0
PARAMETER ptypos              REAL    0.E0
PARAMETER ptzpos              REAL    0.E0
PARAMETER ptmass              REAL    0.E0
PARAMETER sim_tAmbient        REAL    1.E4
PARAMETER sim_tRelax          REAL    0.E0
PARAMETER sim_relaxRate       REAL    0.99E0
PARAMETER larget              REAL    1.E10
PARAMETER sim_BHRadius        REAL    1.E0
PARAMETER sim_profFile        STRING  "sm.dat"

D sim_rhoAmbient     Initial ambient density
PARAMETER sim_rhoAmbient   REAL    1.

D sim_nsubzones       Number of `sub-zones' in cells for applying 1d profile
PARAMETER sim_nsubzones     INTEGER 7

PARAMETER sim_xctr         REAL 0.5E0
PARAMETER sim_yctr         REAL 0.5E0
PARAMETER sim_zctr         REAL 0.5E0

PARAMETER plot_var_13 STRING "none"
PARAMETER plot_var_14 STRING "none"
PARAMETER plot_var_15 STRING "none"
PARAMETER plot_var_16 STRING "none"
PARAMETER plot_var_17 STRING "none"
PARAMETER plot_var_18 STRING "none"
PARAMETER plot_var_19 STRING "none"
PARAMETER plot_var_20 STRING "none"
PARAMETER plot_var_21 STRING "none"
PARAMETER plot_var_22 STRING "none"
PARAMETER plot_var_23 STRING "none"

PPDEFINE HIDE_SANITIZE

