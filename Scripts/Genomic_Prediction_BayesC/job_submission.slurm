#!/bin/sh --login
#SBATCH --time=10:00:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --job-name wv_1kb_sliding
#SBATCH --output=%x_%j.out

module purge
module load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3
cd /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC

#Rscript -e "rmarkdown::render('BayesC_mapping.Rmd', clean=TRUE, output_format='html_document')"
Rscript window_variances.r

scontrol show job $SLURM_JOB_ID
