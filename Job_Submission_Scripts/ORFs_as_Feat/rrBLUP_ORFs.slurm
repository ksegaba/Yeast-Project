#!/bin/sh --login
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name rrBLUP_ORFs_step14_baseline
#SBATCH --output=%x_%j

cd /mnt/scratch/seguraab/yeast_project/ORF_yeast_rrBLUP_results/
module purge
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2

##################### ORF Copy Number #####################
# Step 9. get the population structure
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/08_getPCs.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv
#mv PCA5_geno.csv PCA5_ORFs.csv

# Step 10. genomic prediction using the genetic markers or population structure within a cross-validation scheme
# Note: 09_rrBLUP_fread_PCs.r builds a simple linear regression model using the PCs, which has similar results of models built using mixed.solve from rrBLUP
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all all 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv ORFs
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r PCA5_ORFs_no_NA.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all all 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv pca	

# Steps 11-14 are how I did the feature selection. Markers were selected using the training set, and in the end, models built using selected markers were applied on the test set. Since I did 5-fold CV, thus 1/6 of individuals will be used as test set, and the remaining 5/6 individuals will be used in the 5-fold cv.
# Step 11. hold out individuals for test set, do the stratified sampling
#python 10_holdout_test_stratified.py pheno.csv target_trait 6 # this creates Test.txt but I didn't run this since I'm using the same file for all models

# Step 12. if you want to do feature selection, then you should build models without the test set. So first, get the matrix for the training set, make the CVs file using the training individuals, then build models using the training matrices and output the coef of markers
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/11_split_geno_pheno_fread.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt # need to make pheno_training doesn't change and that ORFs_training matches isolates from geno_training and pheno_training
#python 07_make_CVs.py -file pheno_training.csv -cv 5 -number 10 # I don't need to run this. I'm going to use the same CFVs.csv I've been using for all models
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/geno_training.csv /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/pheno_training.csv all all 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv ORFs
#mv /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/geno_training.csv /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/ORFs_training.csv
#mv /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/ORFs_training.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018

# Step 13. select the number of markers based on the abs coef
# 2/8/22: haven't run this yet
#python 12_select_markers_according_to_abs_coef.py  -coef coef_file -start 250 -stop 5250 -step 250

#TRAITS=( YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDDMSO YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDNYSTATIN YPDANISO50 YPDFLUCONAZOLE YPDSDS YPDBENOMYL200 YPDFORMAMIDE4 YPDSODIUMMETAARSENITE YPDBENOMYL500 YPDFORMAMIDE5 )
#for trait in ${TRAITS[*]}
#do
#echo "${trait} baseline ORF copy number rrBLUP model"
# Step 14. genomic prediction using the genetic markers or population structure within a cross-validation scheme
# The logical for the following script is that: first X% of all the individuals will be held out as test set, which will never be used in the model training process, the remaining 1-X% will be used to train the model, using exactly the same approach as step 10. For each cv fold, the model was also applied to the test set, and after run for all the cv folds, the average r2 across n cv folds will be used for the test set.
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/13_rrBLUP_training_test_split_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all $trait /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv ORFs_baseline
#done
# 2/8/22: will have to create a Write_FS script for feature selction for each trait individually

# Step 15. get the prediction using the top 5 PCs for ramdomly selected markers. Make sure the title of the first column in your ORFs matrix is "ID".
# 2/8/22: haven't run this yet
#Rscript 14_random_select_subset.r ORFs_file start stop step total_number
#Rscript 15_rrBLUP_pca_for_subset_markers.r ORFs_250.csv pheno.csv selected_markers $trait Test.txt 5 10 CVFs.csv Random_250_markers_pca


##################### ORF Presence Absence ##################### # 2/23/22 Will have to rerun the above because I forgot to add pres_abs to save names in step 10 to ORFs and pca
# Step 9. get the population structure
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/08_getPCs.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv
#mv PCA5_geno.csv PCA5_ORFs_pres_abs.csv

# Step 10. genomic prediction using the genetic markers or population structure within a cross-validation scheme
# Note: 09_rrBLUP_fread_PCs.r builds a simple linear regression model using the PCs, which has similar results of models built using mixed.solve from rrBLUP
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all all 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv ORFs_pres_abs # previous just ORFs as save name; possibly messed up the Coef_ORFs_trait.csv files for the previous section (ORF copy number)
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r PCA5_ORFs_pres_abs.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all all 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv pca_pres_abs	

# Steps 11-14 are how I did the feature selection. Markers were selected using the training set, and in the end, models built using selected markers were applied on the test set. Since I did 5-fold CV, thus 1/6 of individuals will be used as test set, and the remaining 5/6 individuals will be used in the 5-fold cv.
# Step 11. hold out individuals for test set, do the stratified sampling
#python 10_holdout_test_stratified.py pheno.csv target_trait 6 # this creates Test.txt but I didn't run this since I'm using the same file for all models

# Step 12. if you want to do feature selection, then you should build models without the test set. So first, get the matrix for the training set, make the CVs file using the training individuals, then build models using the training matrices and output the coef of markers
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/11_split_geno_pheno_fread.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt # need to make pheno_training doesn't change and that ORFs_training matches isolates from geno_training and pheno_training
#python 07_make_CVs.py -file pheno_training.csv -cv 5 -number 10 # I don't need to run this. I'm going to use the same CFVs.csv I've been using for all models
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/geno_training.csv /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/pheno_training.csv all all 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv ORFs_pres_abs
#mv geno_training.csv ORFs_pres_abs_training.csv
#cp ORFs_pres_abs_training.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018

# Step 13. select the number of markers based on the abs coef
# 2/23/22: haven't run this yet
#python 12_select_markers_according_to_abs_coef.py  -coef coef_file -start 250 -stop 5250 -step 250

TRAITS=( YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDDMSO YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDNYSTATIN YPDANISO50 YPDFLUCONAZOLE YPDSDS YPDBENOMYL200 YPDFORMAMIDE4 YPDSODIUMMETAARSENITE YPDBENOMYL500 YPDFORMAMIDE5 )
for trait in ${TRAITS[*]}
do
echo "${trait} baseline ORF presence/absence rrBLUP model"
# Step 14. genomic prediction using the genetic markers or population structure within a cross-validation scheme
# The logical for the following script is that: first X% of all the individuals will be held out as test set, which will never be used in the model training process, the remaining 1-X% will be used to train the model, using exactly the same approach as step 10. For each cv fold, the model was also applied to the test set, and after run for all the cv folds, the average r2 across n cv folds will be used for the test set.
Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/13_rrBLUP_training_test_split_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all $trait /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv ORFs_pres_abs_baseline
done
# 2/23/22: will have to create a Write_FS script for feature selction for each trait individually

# Step 15. get the prediction using the top 5 PCs for ramdomly selected markers. Make sure the title of the first column in your ORFs matrix is "ID".
# 2/23/22: haven't run this yet
#Rscript 14_random_select_subset.r ORFs_file start stop step total_number
#Rscript 15_rrBLUP_pca_for_subset_markers.r ORFs_250.csv pheno.csv selected_markers $trait Test.txt 5 10 CVFs.csv Random_250_markers_pca

scontrol show job $SLURM_JOB_ID