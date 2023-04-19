#!/bin/bash --login 
#SBATCH --array=250-5000:250
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name rf_fs_split1
#SBATCH --output=%x_%A.out

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
cd /mnt/scratch/seguraab/yeast_project/

TRAITS=( YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDSDS YPDSODIUMMETAARSENITE YPDNYSTATIN YPDFLUCONAZOLE YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDDMSO YPDANISO50 YPDBENOMYL200 YPDFORMAMIDE4 YPDBENOMYL500 YPDFORMAMIDE5 )

# Feature selection models using Test_split1.txt
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
path2=/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
start_ind=$SLURM_ARRAY_TASK_ID
end_ind=$(($start_ind+250))
echo "this job starts from $start_ind to $end_ind features"
for trait in ${TRAITS[*]};
do
    echo "$trait"
    #  python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name $trait -sep , -feat ${path2}/feat_rf_${trait}_top_$i -test ${data}/Test_split1.txt -alg RF -n_jobs -1 -n 10 -cv_num 5 -save ${trait}_rf_fs_split1 -plots t
done

scontrol show job $SLURM_JOB_ID



