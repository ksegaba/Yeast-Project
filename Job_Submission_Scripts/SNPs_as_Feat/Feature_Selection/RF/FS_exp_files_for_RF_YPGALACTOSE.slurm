#!/bin/sh --login                 
#SBATCH --array=2,4,8,16,32,64,128,256,512,1024,2048                 
#SBATCH --time=3:00:00                 
#SBATCH --ntasks=3                 
#SBATCH --cpus-per-task=1                 
#SBATCH --mem=20G                 
#SBATCH --job-name RF_make_FS_YPGALACTOSE                 
#SBATCH --output=%x_%j                 
cd /mnt/gs21/scratch/seguraab/yeast_project/yeast_rf_results/                 
module purge                 
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4                 
trait=YPGALACTOSE                 
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline                 
path3=/mnt/gs21/scratch/seguraab/yeast_project/yeast_rf_results
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/
                     python ${path1}/Feature_Selection.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -alg randomforest -n $SLURM_ARRAY_TASK_ID -test ${data}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T
                     python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_$SLURM_ARRAY_TASK_ID -test ${data}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_$SLURM_ARRAY_TASK_ID -plots t


scontrol show job $SLURM_JOB_ID