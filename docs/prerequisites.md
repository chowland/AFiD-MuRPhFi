# Installing prerequisites

AFiD requires the following libraries and programs to be installed before a simulation can be compiled and run:

- A Fortran 90 compiler
- LAPACK
- An MPI implementation (such as OpenMPI or MPICH)
- FFTW
- HDF5 (compiled with the MPI wrappers for the compilers)

For production runs, we recommend consulting your HPC centre's documentation for identifying the best software stack to run with.
This may often lead to highly optimized performance using the LAPACK libraries of the Intel MKL library and the Intel Fortran compiler.

However for small runs, development, and debugging, AFiD can be run using the GNU compiler and standard BLAS and LAPACK libraries.
These dependencies can be installed very easily on Ubuntu, with a single line command.
These instructions also work for machines running Ubuntu in the Windows Subsystem for Linux.

## Single-line prerequisite installation
```
sudo apt install build-essential gfortran libblas-dev liblapack-dev libhdf5-mpich-dev libfftw3-dev
```

## Using a new HPC?
If you are trying to get AFiD up and running on a new HPC, the module system on the cluster can help you get access to all the libraries needed.
To find the modules you need, try running the command `module av hdf5` to search for e.g. HDF5.
This will highlight the available modules providing the library.

Say we found the module `hdf5/1.14-foss-2023a` through this search.
Then running `module show hdf5/1.14-foss-2023a` may also provide shell variables that should be passed to the compiler (such as `HDF5_FFLAGS` or similar).
See some of the existing examples in the `Makefile` for reference.