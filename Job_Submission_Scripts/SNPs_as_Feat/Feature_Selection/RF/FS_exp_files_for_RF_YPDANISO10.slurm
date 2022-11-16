#!/bin/sh --login                 
#SBATCH --array=2,4,8,16,32,64,128,256,512,1024                 
#SBATCH --time=3:00:00                 
#SBATCH --ntasks=3                 
#SBATCH --cpus-per-task=1                 
#SBATCH --mem=20G                 
#SBATCH --job-name RF_make_FS_YPDANISO10                 
#SBATCH --output=%x_%j.out                 
cd /mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/                 
module purge                 
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4                 
trait=YPDANISO10                 
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline                 
path2=/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results                 
path3=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
python ${path1}/Feature_Selection.py -df ${path3}/geno.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
python ${path1}/ML_regression.py -df ${path3}/geno.csv -df2 ${path3}/pheno.csv -y_name ${trait} -sep , -feat ${path2}/feat_exp_rf_${trait}_top_$SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -alg RF -n_jobs 12 -n 10 -cv_num 5 -save ${trait}_exp_rf_$SLURM_ARRAY_TASK_ID -plots t


scontrol show job $SLURM_JOB_ID