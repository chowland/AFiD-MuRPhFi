# AFiD (MultiRes) Quick Start Guide

AFiD is a highly parallel application for simulating canonical flows in a channel domain.
More technical details can be found in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007).
The code is developed jointly by the University of Twente, SURFsara, and the University of Rome "Tor Vergata".

In addition to the method described in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007), this _multi-resolution_ version of AFiD evolves a second scalar field.
This scalar field is simulated on a refined grid to allow for low diffusivity values or high Schmidt numbers.
The multi-resolution method applied is detailed in [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031) and was implemented into AFiD by Chong Shen Ng.
Interpolation between the two grids is performed using a two-point Hermite interpolation scheme, further details of which can be found [here](interpolation.md).

This code is not fully described by either of the existing publications ([van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007) or [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031)).
For clarity, the key differences between this code and those papers is described below.

Differences with [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007):

- The temperature field (along with the salinity field on the refined grid) is no longer co-located with the wall-normal velocity. The scalar fields both take their values at the mid-points of the computational cells, as in [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031)

Differences with [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031):

- Since the current code is modified from AFiD, it is pencil-parallelized in the periodic directions. This contrasts with the previous multi-resolution code, which was slab-parallelized in the wall-normal direction.
- For reasons linked to this change in parallelization, the only terms calculated implicitly are the diffusive terms with derivatives in the wall-normal direction (e.g. $`\nu \pd_{xx}u`$). All other terms are computed explicitly.
- The multiple resolution strategy in time from [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031) is not yet implemented in the code. For now, we only rely on a CFL condition.

## Prerequisites
Before the software can be compiled, the following libraries and packages must be installed:

- A Fortran 90 compiler
- BLAS
- LAPACK
- Szip
- zlib
- MPICH
- FFTW (with `mpi` and `openmp` enabled)
- HDF5 (compiled with the `mpich` wrappers for the compilers)

A step-by-step guide to installing these on Ubuntu can be found [here](prerequisites.md).

## Building AFiD

## Running a simulation

## Post-processing