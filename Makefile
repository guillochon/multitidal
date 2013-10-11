Simulation += Simulation_data.o poly2.o read_table.o Multitidal_interface.o Multitidal_findExtrema.o
Driver += nr.o nrutil.o 

nrutil.o : nrtype.o
Simulation_init.o : Simulation_data.o 
Simulation_initBlock.o : Simulation_data.o
