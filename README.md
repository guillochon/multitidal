Code originally written by James Guillochon. README last updated 11/24/2015.

This simulation setup allows one to simulate the interaction between an extended, polytropic object and a point mass. The code now utilizes the Sinks module to integrate particle trajectories, which brings the benefit of being able to include an arbitrary number of point-like perturbers. The code has been tested to work with FLASH 4.2.2, and is currently being tested with FLASH 4.3, the latest version.

To install, clone this repository to the Simulation path within the FLASH hierarchy:

	cd /path/to/flash/directory/FLASH4.3/source/Simulation/SimulationMain
	hg clone ssh://hg@bitbucket.org/Guillochon/multitidal

Then, change back to the main FLASH directory and run the setup script. A recommended setup for the default parameters are (for pure hydro):

	cd /path/to/flash/directory/FLASH4.3
	./setup MultiTidal -noclobber -maxblocks=1000 -auto -3d -objdir=object_multitidal +newMpole +uhd

Within this repository are a number of "example" run directories that contain files that need to be included when running this setup under various conditions. The example folder that corresponds to the setup line above is `example_setup`, make a full copy of this folder to the location you intend to produce the FLASH outputs.

	cp -R /path/to/flash/directory/FLASH4.3/source/Simulation/SimulationMain/multitidal/example_setup /desired/run/location/.

Then, to run FLASH, use `mpirun` (or similar) from that run directory, e.g.

	cd /desired/run/location
	mpirun -np 4 /path/to/flash/directory/FLASH4.3/object_multitidal/flash4

To add additional physics/features, the user needs to modify the setup line above and use a run directory appropriate to that setup:

* For MHD, substitute `+uhd` with `+usm`.

* To load a stellar profile, add `loadProfile=True` to the setup line. This will switch to the Helmholtz equation of state, which is customized by MultiTidal to use an extended Helmholtz table. If you use this flag, you must copy the extended Helmholtz table `helm_extended_table.dat` to your run directory from the object directory. An example directory for use when `loadProfile=True` is included in the repository, `example_profile_setup`.

* Use `+parallelIO` in conjunction with the other flags (if it is supported by your cluster) for faster read/write of HDF5 files.