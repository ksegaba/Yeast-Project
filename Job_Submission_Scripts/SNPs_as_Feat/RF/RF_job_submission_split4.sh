#!/bin/bash --login 
#SBATCH --array=0-34
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name baseline_split4_rf
#SBATCH --output=%x_%A_%a.out

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
cd /mnt/scratch/seguraab/yeast_project/

TRAITS=( YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDSDS YPDSODIUMMETAARSENITE YPDNYSTATIN YPDFLUCONAZOLE YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDDMSO YPDANISO50 YPDBENOMYL200 YPDFORMAMIDE4 YPDBENOMYL500 YPDFORMAMIDE5 )
echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]}"

# Baseline models using Test_split4.txt
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} -sep , -test ${data}/Test_split4.txt -alg RF -n_jobs 12 -n 10 -cv_num 5 -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_rf_baseline_split4 -plots t


scontrol show job $SLURM_JOB_ID



