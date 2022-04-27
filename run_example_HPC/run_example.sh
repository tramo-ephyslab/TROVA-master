#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 40
#SBATCH -p 'shared'
#SBATCH --qos=shared
#SBATCH --time=4-00:00:00

module --purge
module load cesga/2018
module load miniconda3/4.8.2
conda activate envname

srun -n $SLURM_NTASKS --mpi=pmi2  python TROVA.py input_forw_NATL.cfg >& py_${SLURM_JOB_ID}.log
 
