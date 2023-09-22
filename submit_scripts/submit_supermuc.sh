#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J RBC_Ra1e8
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=user@myemail.nl
# Wall clock limit:
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=myaccount
#SBATCH --partition=micro
#Number of nodes and MPI tasks per node:
#SBATCH --nodes=3
#SBATCH --ntasks=144
#SBATCH --ear=off
 
#Important
module load slurm_setup
module load hdf5 fftw

mpiexec -n $SLURM_NTASKS afid 12 12 
