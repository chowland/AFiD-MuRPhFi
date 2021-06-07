# Installing prerequisites

AFiD requires the following libraries and programs to be installed before a simulation can be compiled and run:

- A Fortran 90 compiler
- BLAS
- LAPACK
- Szip
- zlib
- An MPI implementation (such as OpenMPI or MPICH)
- FFTW
- HDF5 (compiled with the MPI wrappers for the compilers)

For production runs, we recommend using the Intel Fortran compiler `ifort` as well as the BLAS and LAPACK libraries included in the Intel MKL Library.
However for small runs, development, and debugging, AFiD can be run using the GNU compiler and standard BLAS and LAPACK libraries.
These dependencies can be installed very easily on Ubuntu, with a single line command.
These instructions also work for machines running Ubuntu in the Windows Subsystem for Linux.

## Single-line prerequisite installation
```
sudo apt install build-essential gfortran libblas-dev liblapack-dev libhdf5-mpich-dev libfftw3-dev
```