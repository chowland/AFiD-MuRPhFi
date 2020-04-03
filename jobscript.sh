#!/bin/bash

#SBATCH --job-name=pCF

#SBATCH --workdir=.
#SBATCH --nodes=1
#SBATCH --qos=debug
#SBATCH --tasks-per-node=48
#SBATCH -t 2:00:00

#SBATCH --output=out-%j.txt
#SBATCH --error=err-%j.txt

# cd $SLURM_SUBMIT_DIR

# This avoids inadvertent OpenMP threading
export OMP_NUM_THREADS=1
# export SLURM_STATISTICS=yes

mkdir flowmov
srun ./afid
