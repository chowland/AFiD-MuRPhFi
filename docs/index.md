# AFiD-MuRPhFi Quick Start Guide

AFiD is a highly parallel application for simulating canonical flows in a channel domain.
More technical details can be found in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007).
The code is developed jointly by the University of Twente, SURFsara, and the University of Rome "Tor Vergata".

In addition to the method described in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007), this _multi-resolution_ version of AFiD evolves a second scalar field.
This scalar field is simulated on a refined grid to allow for low diffusivity values or high Schmidt numbers.
The multi-resolution method applied is detailed in [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031) and was implemented into AFiD by Chong Shen Ng.
Interpolation between the two grids is performed using a two-point Hermite interpolation scheme, further details of which can be found [here](interpolation.md).

The high resolution grid can also be used to simulate a melting (or dissolving) solid by the phase-field method of [Hester et al (2020)](https://doi.org/10.1098/rspa.2020.0508).

This code is not fully described by either of the existing publications ([van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007) or [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031)).
For clarity, the key differences between this code and those papers is described below.

Differences with [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007):

- The temperature field (along with the salinity field on the refined grid) is no longer co-located with the wall-normal velocity. The scalar fields both take their values at the mid-points of the computational cells, as in [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031)

Differences with [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031):

- Since the current code is modified from AFiD, it is pencil-parallelized in the periodic directions. This contrasts with the previous multi-resolution code, which was slab-parallelized in the wall-normal direction.
- For reasons linked to this change in parallelization, the only terms calculated implicitly are the diffusive terms with derivatives in the wall-normal direction (e.g. $`\nu \partial_{xx}u`$). All other terms are computed explicitly.
- The multiple resolution strategy in time from [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031) is not yet implemented in the code. For now, we only rely on a CFL condition.

## Prerequisites
Before the software can be compiled, the following libraries and packages must be installed:

- A Fortran 90 compiler
- LAPACK
- MPI (OpenMPI/MPICH)
- FFTW (with `mpi` and `openmp` enabled)
- HDF5 (compiled with the `mpich` wrappers for the compilers)

A step-by-step guide to installing these on Ubuntu can be found [here](prerequisites.md).

## Building AFiD

A `Makefile` is provided to easily build the program.
Once all the pre-requisites are installed (along with GNU Make), you only need to run `make` in the command line to build AFiD.
The `Makefile` contains a variable `MACHINE` that automatically enables a range of appropriate compiler options based on the computer you are using.
This variable has a range of preset options for HPC facilities across Europe, on which the code compilation has been successful.

## Running a simulation
Once AFiD has been successfully compiled, an executable `afid` will be produced in the root directory of the repository.
You can then add this directory to your `PATH`, say on Ubuntu by adding the following line to your `.profile` file (assuming AFiD-MuRPhFi has been stored in your home directory):
```
PATH="$PATH:$HOME/AFiD-MuRPhFi"
```
Running a simulation can then be achieved by executing `afid` using `mpiexec` with the following command:
```
mpiexec -n N afid Ny Nz
```
This command tells `mpiexec` to use `N` cores, and tells `afid` that we want to decompose the domain `Ny` times in the y direction and `Nz` times in the z direction.
Note that `N`=`Ny`×`Nz` must be satisfied, and `mpiexec` must use the same MPI implementation that you used to compile `afid`.

A range of SLURM submission files are also be provided in the directory `submit_scripts` as examples for use on HPC systems.

## Post-processing

The `tools` subdirectory contains a range of helper functions in both Python and (experimentally) Julia to enable easy reading and manipulation of the data stored in the output statistics files.
Both modules are called `AFiDTools`.
More details on the helper functions are provided in the documentation, including example Jupyter notebooks and a guide to using ParaView for 3-D visualisation.