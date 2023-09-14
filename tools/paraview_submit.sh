#!/bin/bash
#SBATCH --exclusive
#SBATCH -p rome
#SBATCH -t 1:00:00
#SBATCH --tasks-per-node=128
#SBATCH -N 1
module load 2022
module load ParaView-server-osmesa/5.10.1-foss-2022a-mpi

srun pvserver --force-offscreen-rendering