#!/bin/bash
# Set the simulation to run on 2 nodes from the rome partition
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --partition=rome
# Set a time limit of one hour for the job
# MAKE SURE YOUR TIME LIMIT IN bou.in IS SMALLER THAN THIS!
#SBATCH --time=01:00:00
# Give the job a useful name
#SBATCH --job-name=Ra1e8_RBC
# Email yourself with updates
#SBATCH --mail-user=user@myemail.nl
#SBATCH --mail-type=ALL

# Load modules for MPI and other parallel libraries
module purge
module load 2022
module load foss/2022a
module load HDF5/1.12.2-gompi-2022a

# Finally, call afid, specifying a processor mesh of 16x8
# (this assumes afid can be found in your PATH)
srun afid 16 8