#!/bin/bash 
#SBATCH -J my_mpi_job 
#SBATCH -o my_mpi_job.o%j 
#SBATCH -t 50:00:00 
#SBATCH -N 1 -n 48 ###--ntasks-per-node=16 #changed --tasks-per-node to --ntasks-per-node
#SBATCH --mem=128GB

module purge
module load python/3.9

source /home/gwabeshu/miniconda3/etc/profile.d/conda.sh
conda activate /home/gwabeshu/.conda/envs/py3Env 


export I_MPI_PIN_DOMAIN=omp

PDIR="/workflow/xanthos/"
mpirun -n 20 python ${PDIR}/Basin${DirNum}.py > ${PDIR}/Logger_Basin2.txt 2>&1
