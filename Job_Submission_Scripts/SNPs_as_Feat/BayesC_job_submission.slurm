#!/bin/bash --login 
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name bayesc_write_FS_scripts
#SBATCH --output=%x_%j

cd /mnt/scratch/seguraab/yeast_project/yeast_BayesC_results
module purge
module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2



scontrol show job $SLURM_JOB_ID