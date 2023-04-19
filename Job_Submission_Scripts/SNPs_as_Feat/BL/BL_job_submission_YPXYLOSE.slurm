#!/bin/sh --login                 
#SBATCH --time=15:00:00                 
#SBATCH --ntasks=1                 
#SBATCH --cpus-per-task=2                 
#SBATCH --mem=36G                 
#SBATCH --job-name BL_job_YPXYLOSE                 
#SBATCH --output=%x_%j                 
cd /mnt/scratch/seguraab/yeast_project/yeast_BL_results                 
module purge                 
module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2
Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BL/BL_Model.R /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt YPXYLOSE /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv 5 10 BL_geno /mnt/scratch/seguraab/yeast_project/yeast_BL_results
scontrol show job $SLURM_JOB_ID