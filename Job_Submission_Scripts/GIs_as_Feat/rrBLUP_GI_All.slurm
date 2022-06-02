#!/bin/sh --login
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name rrBLUP_GI_All_create_GI_feat
#SBATCH --output=%x_%j

cd /mnt/scratch/seguraab/yeast_project/GI_yeast_rrBLUP_results/
module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
#module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2

# Create GI composite features
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genetic_Interactions/composite_GI_feat.py

# One-hot encoding
#cut -d "," -f1-1000000 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All.csv > /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All_1.csv
#cut -d "," -f1000001-2000000 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All.csv > /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All_2.csv
#cut -d "," -f2000001-3002178 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All.csv > /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All_3.csv
#python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genetic_Interactions/binary_encode_GI.py
# Step 9. get the population structure
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/08_getPCs.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv

# Note that if you have too large geno matrix, try the following two steps to get the top five PCs
# Step 9_1. remove constant columns in the geno matrix
#python /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/08_1_remove_constant_columns.py /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All.csv

# Step 9_2. get the top 5 PCs
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/08_2_PCA_after_removing.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_GI_All.csv

# Step 10. genomic prediction using the genetic markers or population structure within a cross-validation scheme
# Note: 09_rrBLUP_fread_PCs.r builds a simple linear regression model using the PCs, which has similar results of models built using mixed.solve from rrBLUP
#Rscript 09_rrBLUP_fread_predict_values.r geno.csv pheno.csv all all 5 10 CVFs.csv exome_geno
#Rscript 09_rrBLUP_fread_predict_values.r PCA5_geno.csv pheno.csv all all 5 10 CVFs.csv exome_pca	
# or
#Rscript 09_rrBLUP_fread_PCs.r PCA5_geno.csv pheno.csv all all 5 10 CVFs.csv exome_pca

# Steps 11-14 are how I did the feature selection. Markers were selected using the training set, and in the end, models built using selected markers were applied on the test set. Since I did 5-fold CV, thus 1/6 of individuals will be used as test set, and the remaining 5/6 individuals will be used in the 5-fold cv.
# Step 11. hold out individuals for test set, do the stratified sampling
#python 10_holdout_test_stratified.py pheno.csv target_trait 6

# Step 12. if you want to do feature selection, then you should build models without the test set. So first, get the matrix for the training set, make the CVs file using the training individuals, then build models using the training matrices and output the coef of markers
#Rscript 11_split_geno_pheno.r geno.csv pheno.csv Test.txt
#python 07_make_CVs.py -file pheno_training.csv -cv 5 -number 10
#Rscript 09_rrBLUP.r geno_training.csv pheno_training.csv all target_trait 5 10 CVFs.csv exome_geno

# Step 13. select the number of markers based on the abs coef
#python 12_select_markers_according_to_abs_coef.py  -coef coef_file -start 250 -stop 5250 -step 250

# Step 14. genomic prediction using the genetic markers or population structure within a cross-validation scheme
# The logical for the following script is that: first X% of all the individuals will be held out as test set, which will never be used in the model training process, the remaining 1-X% will be used to train the model, using exactly the same approach as step 10. For each cv fold, the model was also applied to the test set, and after run for all the cv folds, the average r2 across n cv folds will be used for the test set.
#Rscript 13_rrBLUP_training_test_split_fread_predict_values.r geno.csv pheno.csv selected_markers target_trait Test.txt 5 10 CVFs.csv selected_markers_geno

# Step 15. get the prediction using the top 5 PCs for ramdomly selected markers. Make sure the title of the first column in your geno matrix is "ID".
#Rscript 14_random_select_subset.r geno_file start stop step total_number
#Rscript 15_rrBLUP_pca_for_subset_markers.r geno_250.csv pheno.csv selected_markers target_trait Test.txt 5 10 CVFs.csv Random_250_markers_pca

