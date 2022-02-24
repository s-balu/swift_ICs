2LPTic - base code to compute 2LPT displacement field for cosmo ICs

This creates a uniform resolution distribution of particles with positions,
velocities, IDs, and masses for a cosmological simulation with 2nd order Lagrangian perturbation theory. It writes to a GADGET SnapFormat=1 FORTRAN readable binary. It can take a CAMB input file for the transfer function.

Changes to come: 

* Write outputs in HDF5 format

I need to load the following to get the code to compile properly
* fftw2/2.1.5 
* openmpi/1.6.3
* gsl/1.15

To run an example, type

mpirun ./2LPTic ./L100_N512.param

after editing the param file parameters

* Nmesh - for 512 particles on a side, set to 512
* Nsample - for 512 particles on a side, set to 512
* FileBase - root of the output file name
* OutputDir - path to directory containing output files
* OutputFormat - 3 produces HDF5, <3 produces GADGET2 binary
* GlassTileFac - path to glass file - located in source code directory
* FileWithInputSpectrum  - path to transfer file, located in source code directory
