#!/bin/bash --login 
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name rf_write_FS_scripts
#SBATCH --output=%x_%j

#### NOTE: This slurm submission script contains steps to run random forest (baseline and feature selection) on yeast genotype data
cd /mnt/scratch/seguraab/yeast_project/

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4

# step 0. reformat data
#python /mnt/home/seguraab/Shiu_Lab/Project/00_Reformat_data_ML-Pipeline.py

# step 1. define a test set
#echo "Define test set"
TRAITS=( YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDSDS YPDSODIUMMETAARSENITE YPDNYSTATIN YPDFLUCONAZOLE YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDDMSO YPDANISO50  YPDBENOMYL200 YPDFORMAMIDE4 YPDBENOMYL500 YPDFORMAMIDE5 )

for trait in ${TRAITS[*]} # for each trait (condition) 
do
#echo "${trait}"
# I want to use the same test set as in rrBLUP so I need to make sure that test_set.py does this. If not, then I'll just reuse my rrBLUP test set
###python /mnt/home/seguraab/Shiu_Lab/Project/ML-Pipeline/test_set.py -df geno_rf_${trait}.csv -type r -p 0.16 -y_name Y -sep , -save test_rf_${trait}.txt
# The code will result in a random test set each time therefore I will not use this code
# So, test_rf_${trait}.txt is way different from Test.txt, hence I will use Test.txt in Step 3


# step 2. feature selection (best predictors; not for baseline)
echo "Feature selection for ${trait}"
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
path2=/mnt/scratch/seguraab/yeast_project/yeast_rf_results
path3=/mnt/scratch/seguraab/yeast_project/yeast_rf_results
save=/mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -trait ${trait} -start 500 -stop 50000 -step 500 -runFS n -runRF y
done


# step 3. train and apply classification machine learning model (Random Forest)
#echo "Train and apply RF classification ML model"
#for trait in ${TRAITS[*]}
#do
#echo "Apply RF model to top 5000 features for ${trait}"
#python /mnt/home/seguraab/Shiu_Lab/Project/ML-Pipeline/ML_regression.py -df geno_rf_${trait}.csv -sep , -feat feat_rf_${trait}_top_5000 -test test_rf_${trait}.txt -alg RF -n_jobs 12 -cv_num 5 -save ${trait}_rf -plots t

#echo "Baseline for ${trait}"
# this step only completes about 3 traits in 100 hrs!! That's insane!
#python /mnt/home/seguraab/Shiu_Lab/Project/ML-Pipeline/ML_regression.py -df geno_rf_${trait}.csv -sep , -test /mnt/home/seguraab/Shiu_Lab/Project/Test.txt  -alg RF -n_jobs 12 -cv_num 5 -save ${trait}_rf_baseline -plots t
#done

scontrol show job $SLURM_JOB_ID

