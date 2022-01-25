#!/bin/bash --login 
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name geno_data_plots
#SBATCH --output=%x_%j

cd /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis
module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
python geno_data_plots.py
scontrol show job $SLURM_JOB_ID