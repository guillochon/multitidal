Simulation += Simulation_data.o com_accel.o orbit.o poly2.o
Driver += Orbit_update.o bsstep.o odeint.o mmid.o pzextr.o nrutil.o nr.o nrtype.o Bound_mass.o Orbit_energy.o
Grid += gr_zoneOldPotential.o Grid_findExtrema.o
Gravity += Gravity_sendOutputData.o
IO += IO_writeOrbitInfo.o
Hydro += fss.o

bsstep.o : nrutil.o
nrutil.o : nrtype.o
odeint.o : nrutil.o
mmid.o : nrutil.o
pzextr.o : nrutil.o
Orbit_update.o : odeint.o
Simulation_init.o : Simulation_data.o 
Simulation_initBlock.o : Simulation_data.o
