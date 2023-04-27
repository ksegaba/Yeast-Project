#!/bin/bash --login 
#SBATCH --array=0-34 # 35 array tasks
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name shap_values
#SBATCH --output=%x_%j.out

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4

cd /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis

scratch=/mnt/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
TRAITS=($(awk  -F'\t' '{if (NR!=1) print $1}' ../../Results/RESULTS_RF_SNPs_FS.txt))
SNP_FEATS=($(awk -F'\t' '{if (NR!=1) print $9}' ../../Results/RESULTS_RF_SNPs_FS.txt))

echo "Calculating SHAP values for SNP RF FS models"
echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]} ; ${FEATS[${SLURM_ARRAY_TASK_ID}]}"

### SNP RF FS models top ~20 features (some a little less or a little more)
# python SHAP_training.py \
#     -df ${data}/geno.csv -sep , \
#     -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt \
#     -feat ${scratch}/feat*_rf_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_top_${SNP_FEATS[${SLURM_ARRAY_TASK_ID}]} \
#     -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}*_rf_${SNP_FEATS[${SLURM_ARRAY_TASK_ID}]}_models.pkl \
#     -top 10 -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_snp_rf_${SNP_FEATS[${SLURM_ARRAY_TASK_ID}]}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
scratch=/mnt/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs
TRAITS=($(awk  -F'\t' '{if (NR!=1) print $1}' ../../Results/RESULTS_RF_ORFs_FS.txt))
ORF_FEATS=($(awk -F'\t' '{if (NR!=1) print $9}' ../../Results/RESULTS_RF_ORFs_FS.txt))
echo "Calculating SHAP values for ORF presence/absence RF FS models"
echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]}"

### ORF RF FS
python SHAP_training.py \
    -df ${data}/ORFs_pres_abs.csv -sep , \
    -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -test ${data}/Test.txt \
    -feat ${scratch}/feat_rf_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_top_${ORF_FEATS[${SLURM_ARRAY_TASK_ID}]} \
    -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf_${ORF_FEATS[${SLURM_ARRAY_TASK_ID}]}_models.pkl \
    -top 10 -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf_rf_${ORF_FEATS[${SLURM_ARRAY_TASK_ID}]}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "Calculating SHAP values for ORF copy number RF FS models"
TRAITS=($(awk  -F'\t' '{if (NR!=1) print $1}' ../../Results/RESULTS_RF_CNVs_FS.txt))
CNV_FEATS=($(awk -F'\t' '{if (NR!=1) print $9}' ../../Results/RESULTS_RF_CNVs_FS.txt))

### CNO RF FS
python SHAP_training.py \
    -df ${data}/ORFs_no_NA.csv -sep , \
    -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -test ${data}/Test.txt \
    -feat ${scratch}/feat_rf_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_top_${CNV_FEATS[${SLURM_ARRAY_TASK_ID}]} \
    -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cno_${CNV_FEATS[${SLURM_ARRAY_TASK_ID}]}_models.pkl \
    -top 10 -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cnv_${CNV_FEATS[${SLURM_ARRAY_TASK_ID}]}

scontrol show job $SLURM_JOB_ID
