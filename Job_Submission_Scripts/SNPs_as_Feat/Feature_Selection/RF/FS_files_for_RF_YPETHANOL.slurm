#!/bin/sh --login                 
#SBATCH --array=1000-40000:5000                 
#SBATCH --time=100:00:00                 
#SBATCH --ntasks=3                 
#SBATCH --cpus-per-task=1                 
#SBATCH --mem=20G                 
#SBATCH --job-name RF_make_FS_YPETHANOL                 
#SBATCH --output=%x_%j.out                 
cd /mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/                 
module purge                 
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4                 
trait=YPETHANOL                 
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline                 
path3=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018                 
start_ind=$SLURM_ARRAY_TASK_ID                 
end_ind=$(( $start_ind + 5000 - 1000))                  
echo "this job start from $start_ind to $end_ind"                 
for i in `seq $start_ind 1000 $end_ind`; do
python ${path1}/Feature_Selection.py -df ${path3}/geno.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $i -test ${path3}/Test.txt -sep , -type r -save feat_rf_${trait}_top -scores T
done
scontrol show job $SLURM_JOB_ID