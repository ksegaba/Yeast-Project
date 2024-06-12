#!/usr/bin/bash
#SBATCH --array=1-35
#SBATCH --time=150:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=MxE_models
#SBATCH --output=../logs/%x_%A_%a.out

module load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3
cd /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GxE/code

X_file=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv
test_file=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt

TRAITS=($(head -n 1 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv | sed 's/,/ /g'))

## CV1 Marker-by-Environment models on target trait and top 4 corr envs (submitted May 19, 2024)
prefix=GxE_cv1_pav_5_envs
save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/markerXenv/cv1
Rscript 01_GXE_models_multialg_Kenia_noslurm.R $X_file all $test_file $prefix $save_dir mxe 10 cv1 n ${TRAITS[${SLURM_ARRAY_TASK_ID}]}

## CV2 Marker-by-Environment models on target trait and top 4 corr envs (submitted --- --, 2024)
#prefix=GxE_cv2_pav_5_envs
#save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/markerXenv/cv2
#Rscript 01_GXE_models_multialg_Kenia_noslurm.R $X_file all $test_file $prefix $save_dir mxe 10 cv2 n ${TRAITS[${SLURM_ARRAY_TASK_ID}]}

############### THESE MODELS WERE RUN ON THE COMMAND LINE ###############
## CV1 Multi-Trait models on target trait and top 4 corr envs (submitted April 23, 2024)
#save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/multi-trait/cv1
#sbatch 01_GXE_models_multialg_Kenia.R $X_file all $test_file $prefix $save_dir multi 10 cv1 n

## CV2 Multi-Trait models on target trait and top 4 corr envs (submitted April 23, 2024)
#save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/multi-trait/cv2
#sbatch 01_GXE_models_multialg_Kenia.R $X_file all $test_file $prefix $save_dir multi 10 cv2 n

################################### OLD RUNS ###################################
### ORF presence/absence (submitted May 4, 2023)
## CV1 Marker-by-Environment models on trait pairs
#X_file=Data/Peter_2018/ORFs_pres_abs.csv
#test_file=Data/Peter_2018/Test.txt
#prefix=GxE_cv1_orf_baseline
#save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/markerXenv/cv1
#sbatch 01_GXE_models_Kenia.R $X_file all $test_file $prefix $save_dir 20 cv1 n

### ORF presence/absence
## CV1 Across-Environment models on target trait and top 4 correlated environments (submitted July 6, 2023)
#X_file=Data/Peter_2018/ORFs_pres_abs.csv
#test_file=Data/Peter_2018/Test.txt
#prefix=GxE_cv1_5envs_orf_baseline
#save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/across-env/cv1
#sbatch 01_GXE_models_multialg_Kenia.R $X_file all $test_file $prefix $save_dir across 20 cv1 n

## CV1 Marker-by-Environment models on target trait and top 4 corr envs (submitted July 9, 2023)
#save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/markerXenv/cv1
#sbatch 01_GXE_models_multialg_Kenia.R $X_file all $test_file $prefix $save_dir mxe 20 cv1 n

## CV1 Multi-Trait models on target trait and top 4 corr envs (submitted July 12, 2023)
#save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/multi-trait/cv1
#sbatch 01_GXE_models_multialg_Kenia.R $X_file all $test_file $prefix $save_dir multi 20 cv1 n