# AFiD-MuRPhFi
**A**dvanced **Fi**nite **D**ifference flow solver with **Mu**ltiple **R**esolution and **Ph**ase-**Fi**eld implementations

![Build status](https://github.com/chowland/AFiD-MuRPhFi/actions/workflows/CI.yml/badge.svg)
[![Documentation](https://img.shields.io/badge/documentation-in%20progress-blue)](https://chowland.github.io/AFiD-MuRPhFi/)

## Quickstart

AFiD has the following prerequisites that must be installed to compile and run the program:
- A parallel HDF5 library wrapping an MPI implementation and a Fortran 90 compiler
- The numerical LAPACK library (or Intel MKL)
- FFTW3 (*not* the MKL implementation)

On Ubuntu, these prerequisites can be installed using the single line command
```
sudo apt install build-essential gfortran liblapack-dev libhdf5-openmpi-dev libfftw3-dev
```
AFiD can then be compiled simply by running `make`.
The [Makefile](./Makefile) contains a `MACHINE` variable and a `FLAVOUR` variable that the user can modify to reflect the libraries installed (e.g. setting `FLAVOUR=Intel` ensures the MKL libraries and Intel compiler options are used).
Presets are also available for a collection of European HPC systems using the `MACHINE` variable.

## Code description

### Key features
- Third-order Runge-Kutta time stepping
- Finite-differences calculate spatial derivatives on staggered velocity grid
- Pencil MPI decomposition based on the 2DECOMP&FFT library
- Cubic Hermite interpolation between two decoupled grids for multiple resolution
- Phase-field model to simulate the flow around melting objects following [Hester et al. (2020)](https://doi.org/10.1098/rspa.2020.0508)
- Immersed boundary method for fixed objects following [Fadlun et al. (2000)](https://doi.org/10.1006/jcph.2000.6484)
- A simple moisture model for a rapidly condensing vapour field following [Vallis et al. (2019)](https://doi.org/10.1017/jfm.2018.954)

AFiD is a highly parallel application for simulating canonical flows in a channel domain.
More technical details can be found in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007).
The code is developed jointly by the University of Twente, SURFsara, and the University of Rome "Tor Vergata".

In addition to the method described in [van der Poel et al (2015)](https://doi.org/10.1016/j.compfluid.2015.04.007), this *multi-resolution* version of AFiD evolves a second scalar field.
This scalar field is simulated on a refined grid to allow for low diffusivity values or high Schmidt numbers.
The multi-resolution method applied is detailed in [Ostilla-Monico et al (2015)](https://doi.org/10.1016/j.jcp.2015.08.031) and was implemented into AFiD by [@chongshenng](https://github.com/chongshenng).
Interpolation between the two grids is performed using a four-point Hermite interpolation scheme.

One key difference from the original version of AFiD is that the scalar fields are evaluated at the mid-points of the computational cells.
This staggered grid provides more flexibility for applications with different gravitational axes.

The refined grid can also be used to evolve a phase-field variable to model the melting and dissolving of immersed solid objects.

## Contributing
If you would like to contribute to bug fixing/feature development/documentation, please create a new branch, commit your changes and then submit a pull request.
