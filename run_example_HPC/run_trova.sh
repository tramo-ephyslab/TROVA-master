#!/bin/bash -l

#SBATCH --mem=120GB
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:10:00

module --purge
module load cesga/2020
module load miniconda3/4.9.2
conda activate envname

srun -n $SLURM_NTASKS --mpi=pmi2 python  run_TROVA.py  $1 >> py_${SLURM_JOB_ID}.log

