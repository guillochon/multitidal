Code originally written by James Guillochon. README last updated 10/10/2015.

This simulation setup allows one to simulate the interaction between an extended, polytropic object and a point mass. The code now utilizes the Sinks module to integrate particle trajectories, which brings the benefit of being able to include an arbitrary number of point-like perturbers. The code has been tested to work with FLASH 4.2.2.

A recommended setup for the default parameters are (for pure hydro):

	./setup MultiTidal -noclobber -maxblocks=1000 -auto -3d -objdir=object_multitidal +newMpole +uhd

For MHD, substitute `+uhd` with `+usm`.

Use `+parallelIO` if it is supported for faster read/write of HDF5 files.