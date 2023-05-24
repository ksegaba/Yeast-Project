#!/usr/bin/bash

module load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3

### ORF presence absence (submitted May 4, 2023)
## CV1 Marker-by-Environment models
X_file=Data/Peter_2018/ORFs_pres_abs.csv
test_file=Data/Peter_2018/Test.txt
prefix=GxE_cv1_orf_baseline
save_dir=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/markerXenv/cv1
sbatch GXE_models_Kenia.R $X_file all $test_file $prefix $save_dir 20 cv1 n

