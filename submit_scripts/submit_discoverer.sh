#!/bin/bash
#SBATCH --job-name="RBC_Ra1e8"
#SBATCH --partition=cn
# Ask for 8 nodes of 128 cores each
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=8
# Set time limit to 12 hours (in seconds)
#SBATCH --time=43200
# Ask for exclusive use of the nodes you are allocated
#SBATCH --exclusive

module purge
module load hdf5/1/1.14/1.14.0-gcc-openmpi
module load lapack/latest-gcc
module load fftw/3/latest-gcc-openmpi

srun ~/AFiD-MuRPhFi/afid 32 32
