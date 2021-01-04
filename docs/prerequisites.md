# Installing prerequisites

AFiD requires the following libraries and programs to be installed before a simulation can be compiled and run:

- A Fortran 90 compiler
- BLAS
- LAPACK
- Szip
- zlib
- MPICH
- FFTW (with `mpi` and `openmp` enabled)
- HDF5 (compiled with the `mpich` wrappers for the compilers)

For production runs, we recommend using the Intel Fortran compiler `ifort` as well as the BLAS and LAPACK libraries included in the Intel MKL Library.
Below is a step-by-step guide to installing these prerequisites on a local machine running Ubuntu.
These instructions also work for machines running Ubuntu in the Windows Subsystem for Linux.

## Installing the Intel compilers and libraries

The Intel compilers can be downloaded from the [Intel Developer Zone](https://software.intel.com/content/www/us/en/develop/tools/parallel-studio-xe/choose-download.html) as part of Parallel Studio XE.
This software is free to download for students, classroom educators, and open source contributors.
Students simply need to register for an account on the Intel website with their academic email address.

Installation of Parallel Studio XE is pretty straightforward, and simply requires you to run the install script `install.sh`.
The whole suite of software is quite large, and the installation can be customized to save space.
If you choose only to install particular components, make sure that the Intel C/C++/Fortran compilers are all installed, along with the MKL libraries.
Once these are installed (by default in `/home/username/intel`), add the following line to your `.bashrc`:

    source /home/username/intel/bin/compilervars.sh intel64

This will ensure that the compilers and libraries are all properly linked to on your `PATH`.

*If you are working on the Windows Subsystem for Linux, the above line to link the libraries may produce an error. This may be fixed by replacing the final line of the function `remove_duplicate_paths` in `compilervars.sh` with*

    export "$arg"="${fixarg}"

Now all the Intel compilers should be installed and accessible to your system.
Test you can reach them by running `which icc` or `ifort --version`.
[fortran-lang.org](https://fortran-lang.org/learn/quickstart) has some nice, simple examples to remind you of the syntax and to test your compiler with.

## Compression libraries (Szip and zlib)

### Szip
As part of its HDF5 input/output, the DDC code takes advantage of the szip and zlib compression libraries.
We need to build these first so they can be linked to when we build the HDF5 library.
For simplicity, we will build all of the libraries we need in a directory at `/home/username/lib_papero`.
The Szip source code can be downloaded from the [HDF5 website](https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz).
We then build the library using the following commands in the terminal (making sure that we use the Intel compilers):

    export CC=icc CXX=icpc FC=ifort
    export CFLAGS='-O3 -xHost -ip' CXXFLAGS='-O3 -xHost -ip' FCFLAGS='-O3 -xHost -ip'
    tar -xvzf szip-2.1.1.tar.gz
    cd szip-2.1.1
    ./configure --prefix=/home/username/lib_papero/szip
    make
    make check
    make install

Ensure the Szip library can be found by adding its location to the `LD_LIBRARY_PATH` variable. This can be done by adding the following line to `.bashrc`:

    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/username/lib_papero/szip/lib"

### zlib
We similarly build the zlib compression by downloading the [source code](http://zlib.net/zlib-1.2.11.tar.gz), unpacking the archive, and compiling it with the Intel compiler:

    export CC=icc CFLAGS='-O3 -xHost -ip'
    tar -xvzf zlib-1.2.11.tar.gz
    cd zlib-1.2.11
    ./configure --prefix=/home/username/lib_papero/zlib
    make
    make check
    make install

## MPICH

MPICH is an implementation of the MPI standard for parallel computing.
The source code can be downloaded from the [MPICH website](https://www.mpich.org/downloads/), and then built using the Intel compilers.
This provides us with an MPI wrapper for the Fortran compiler: `mpifort`.

    export CC=icc CXX=icpc
    export F77=ifort FC=ifort
    tar -xvzf mpich-3.3.2.tar.gz
    cd mpich-3.3.2
    ./configure --enable-shared --prefix=/home/username/lib_papero/mpich
    make
    make install

The installation can be tested by also running `make installcheck`.
You should now have `mpifort` installed in the folder `/home/username/lib_papero/mpich/bin`.
We do not need to add it to the path, but it's good to know its location for when we build the HDF5 library later.

## FFTW

Intel provides its own FFT subroutines in the Math Kernel Library, but these do not support the calls used in our DDC code.
We therefore need to install our own FFTW library: downloading the [source code](http://fftw.org/download.html), and building it with the relevant options as follows.

    ./configure --enable-mpi --enable-openmp --prefix=/home/username/lib_papero/fftw
    make
    make install

## HDF5

Now we have all the prerequisites ready to finally build the HDF5 library.
This library provides a further wrapper around the Fortran compiler `h5pfc` that will link to both the HDF5 library and the MPICH installation above (if we do this next step correctly).

    export CC=/home/username/lib_papero/mpich/bin/mpicc
    export F9X=/home/username/lib_papero/mpich/bin/mpifort
    export FC=/home/username/lib_papero/mpich/bin/mpifort
    ./configure --prefix=/home/username/lib_papero/hdf5 --enable-hl --enable-fortran --enable-parallel --enable-shared --enable-build-mode=production --with-szlib=/home/username/lib_papero/szip --with-zlib=/home/username/lib_papero/zlib
    make
    make -i check   # Note this step may take a long time!
    make install

As with the MPICH installation, you can check the installation by running `make installcheck`.
Finally, we want to add the HDF5 wrappers for the compilers to the path, and we do this by adding the following line to the `.bashrc` file:

    export PATH="$PATH:/home/username/lib_papero/hdf5/bin

Now everything should be installed and linked appropriately, ready for running the DDC code.
A final check of the HDF5 installation can be seen by running `h5pfc -showconfig`, which should produce something like the following:

                SUMMARY OF THE HDF5 CONFIGURATION
                =================================

    General Information:
    -------------------
                       HDF5 Version: 1.12.0
                      Configured on: Tue Sep 15 16:17:24 CEST 2020
                      Configured by: cjh225@Megalodon
                        Host system: x86_64-unknown-linux-gnu
                  Uname information: Linux Megalodon 4.19.104-microsoft-standard #1 SMP Wed Feb 19 06:37:35 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux
                           Byte sex: little-endian
                 Installation point: /home/cjh225/lib_afid/hdf5

    Compiling Options:
    ------------------
                         Build Mode: production
                  Debugging Symbols: no
                            Asserts: no
                          Profiling: no
                 Optimization Level: high

    Linking Options:
    ----------------
                          Libraries: static, shared
      Statically Linked Executables:
                            LDFLAGS:
                         H5_LDFLAGS:
                         AM_LDFLAGS:  -L/home/cjh225/lib_afid/zlib/lib -L/home/cjh225/lib_afid/szip/lib
                    Extra libraries: -lsz -lz -ldl -lm
                           Archiver: ar
                           AR_FLAGS: cr
                             Ranlib: ranlib

    Languages:
    ----------
                                  C: yes
                         C Compiler: /home/cjh225/lib_afid/mpich/bin/mpicc ( MPICH version 3.3.2 built with icc version 19.1.2.254 (gcc version 9.3.0 compatibility))
                           CPPFLAGS:
                        H5_CPPFLAGS: -D_GNU_SOURCE -D_POSIX_C_SOURCE=200809L   -DNDEBUG -UH5_DEBUG_API
                        AM_CPPFLAGS:  -I/home/cjh225/lib_afid/zlib/include -I/home/cjh225/lib_afid/szip/include
                            C Flags:
                         H5 C Flags:   -std=c99 -Wall -Wcheck  -Wl,-s  -O3
                         AM C Flags:
                   Shared C Library: yes
                   Static C Library: yes


                            Fortran: yes
                   Fortran Compiler: /home/cjh225/lib_afid/mpich/bin/mpifort ( Intel(R) Fortran Intel(R) 64 Compiler Version 19.1.2.254 Build 20200623)
                      Fortran Flags:
                   H5 Fortran Flags:    -O3
                   AM Fortran Flags:
             Shared Fortran Library: yes
             Static Fortran Library: yes

                                C++: no

                               Java: no


    Features:
    ---------
                       Parallel HDF5: yes
    Parallel Filtered Dataset Writes: yes
                  Large Parallel I/O: yes
                  High-level library: yes
                    Build HDF5 Tests: yes
                    Build HDF5 Tools: yes
                        Threadsafety: no
                 Default API mapping: v112
      With deprecated public symbols: yes
              I/O filters (external): deflate(zlib),szip(encoder)
                                 MPE:
                       Map (H5M) API: no
                          Direct VFD: no
                  (Read-Only) S3 VFD: no
                (Read-Only) HDFS VFD: no
                             dmalloc: no
      Packages w/ extra debug output: none
                         API tracing: no
                Using memory checker: no
     Memory allocation sanity checks: no
              Function stack tracing: no
           Strict file format checks: no
        Optimization instrumentation: no