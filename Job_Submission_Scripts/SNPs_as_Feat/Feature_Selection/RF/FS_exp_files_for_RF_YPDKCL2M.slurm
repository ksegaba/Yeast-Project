#!/bin/sh --login                 
#SBATCH --time=10:00:00                 
#SBATCH --ntasks=3                 
#SBATCH --cpus-per-task=1                 
#SBATCH --mem=20G                 
#SBATCH --job-name RF_make_FS_YPDKCL2M                 
#SBATCH --output=%x_%j                 
cd /mnt/gs21/scratch/seguraab/yeast_project/yeast_rf_results/                 
module purge                 
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4                 
trait=YPDKCL2M                 
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline                 
path3=/mnt/gs21/scratch/seguraab/yeast_project/yeast_rf_results
"""data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 2 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_2 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_2 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 4 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_4 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_4 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 8 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_8 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_8 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 16 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_16 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_16 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 32 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_32 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_32 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 64 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_64 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_64 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 128 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_128 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_128 -plots t
"""
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 256 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_256 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_256 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 512 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_512 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_512 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 1024 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_1024 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_1024 -plots t

data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/                         
python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n 2048 -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                         
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_2048 -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_2048 -plots t


scontrol show job $SLURM_JOB_ID