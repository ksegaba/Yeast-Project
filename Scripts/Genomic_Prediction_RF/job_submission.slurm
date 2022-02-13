#!/bin/bash --login 
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name map_RF_feat_to_GO
#SBATCH --output=%x_%j

cd /mnt/home/seguraab/Shiu_Lab/Project/yeast_RF_results/SNPs_as_Features
module purge
module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2

Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/gene_set_enrichment.R

scontrol show job $SLURM_JOB_ID