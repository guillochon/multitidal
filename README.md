Code originally written by James Guillochon. README last updated 10/10/2015.

This simulation setup allows one to simulate the interaction between an extended, polytropic object and a point mass. The code now utilizes the Sinks module to integrate particle trajectories, which brings the benefit of being able to include an arbitrary number of point-like perturbers. The code has been tested to work with FLASH 4.2.2.

A recommended setup for the default parameters are (for pure hydro):

	./setup MultiTidal -noclobber -maxblocks=1000 -auto -3d -objdir=object_multitidal +newMpole +uhd

One file must be copied over to the run directory for the default setup: `SpeciesList.txt`, which can be found in the object directory after running the setup command.

For MHD, substitute `+uhd` with `+usm`.

Use `+parallelIO` if it is supported for faster read/write of HDF5 files.

To load a stellar profile, add `loadProfile=True` to the setup line. This will switch to the Helmholtz equation of state, which is customized by MultiTidal to use an extended Helmholtz table. If you use this flag, you must copy the extended Helmholtz table `helm_extended_table.dat` to your run directory from the object directory.