#!/bin/bash --login 
#SBATCH --array=0-34
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name shap_fs_values
#SBATCH --output=%x_%j.out

# ######################## GET RF MODEL SKLEARN VERSIONS #########################
# ## First, check RF model performances for PCs, SNPs, PAVs, and CNVs
# snp_path=/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results
# for file in ${snp_path}/baseline/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${snp_path}/baseline_test1/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test1" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${snp_path}/baseline_test2/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test2" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${snp_path}/baseline_test3/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test3" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${snp_path}/baseline_test4/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test4" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${snp_path}/baseline_test5/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test5" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${snp_path}/base_rand/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    base_rand" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${snp_path}/fs/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    fs" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done

# orf_path=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results
# for file in ${orf_path}/baseline/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${orf_path}/baseline_test1/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test1" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${orf_path}/baseline_test2/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test2" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${orf_path}/baseline_test3/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test3" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${orf_path}/baseline_test4/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test4" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${orf_path}/baseline_test5/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    baseline_test5" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${orf_path}/base_rand/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    base_rand" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done
# for file in ${orf_path}/fs/*.pkl; do
#     sklearn_v=$(python -m pickletools $file | grep -A 2 -F _sklearn_version | grep -P [0-9].[0-9]+.[0-9] -o)
#     echo "${file}   ${sklearn_v}    fs" >> /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_2/RF_model_versions_v2.tsv
# done

# ############################ CALCULATE SHAP VALUES #############################
shap=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018

## BASELINE MODELS
mapfile -t snp_models < /mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/baseline/snp_best_models.txt
snp_mod=${snp_models[${SLURM_ARRAY_TASK_ID}]}
snp_trait=$(echo $snp_mod | grep -P YP[A-Z0-9]+ -o)
mapfile -t pav_models < /mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline/pav_best_models.txt
pav_mod=${pav_models[${SLURM_ARRAY_TASK_ID}]}
pav_trait=$(echo $pav_mod | grep -P YP[A-Z0-9]+ -o)
mapfile -t cnv_models < /mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline/cnv_best_models.txt
cnv_mod=${cnv_models[${SLURM_ARRAY_TASK_ID}]}
cnv_trait=$(echo $cnv_mod | grep -P YP[A-Z0-9]+ -o)

conda activate /mnt/home/seguraab/miniconda3/envs/shap

### SNP RF baseline models top 20 features (sklearn version 1.2.2)
echo "SNPs RF baseline -----------------"
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP/SNP/baseline
scratch=/mnt/scratch/seguraab/yeast_project/SNP_yeast_RF_results/baseline
python ${shap}/SHAP_training_saving_interaction.py \
    -df ${data}/geno.csv \
    -df2 ${data}/pheno.csv \
    -y_name ${snp_trait} \
    -sep , -test ${data}/Test.txt -feat all \
    -model ${scratch}/${snp_mod} \
    -top 20 -interaction n -interaction_score n \
    -save ${snp_trait}_snp_rf_baseline

### PAV RF baseline models top 20 features (sklearn version 1.2.2)
echo "PAVs RF baseline -----------------"
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP/PAV/baseline
scratch=/mnt/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline
python ${shap}/SHAP_training_saving_interaction.py \
    -df ${data}/ORFs_pres_abs.csv \
    -df2 ${data}/pheno.csv \
    -y_name ${pav_trait} \
    -sep , -test ${data}/Test.txt -feat all \
    -model ${scratch}/${pav_mod} \
    -top 20 -interaction n -interaction_score n \
    -save ${pav_trait}_pav_rf_baseline

### CNV RF baseline models top 20 features (sklearn version 1.2.2)
echo "CNVs RF baseline -----------------"
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP/CNV/baseline
scratch=/mnt/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline
python ${shap}/SHAP_training_saving_interaction.py \
    -df ${data}/ORFs_no_NA.csv \
    -df2 ${data}/pheno.csv \
    -y_name ${cnv_trait} \
    -sep , -test ${data}/Test.txt -feat all \
    -model ${scratch}/${cnv_mod} \
    -top 20 -interaction n -interaction_score n \
    -save ${cnv_trait}_cnv_rf_baseline}

conda deactivate

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## FEATURE SELECTION MODELS
mapfile -t snp_models < /mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs/snp_best_fs_models.txt
snp_mod=${snp_models[${SLURM_ARRAY_TASK_ID}]}
snp_trait=$(echo $snp_mod | grep -P YP[A-Z0-9]+ -o)
snp_feat=$(echo $snp_mod | grep -P '(?<=top_)[0-9]+(?=_models)' -o)
mapfile -t pav_models < /mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs/pav_best_fs_models.txt
pav_mod=${pav_models[${SLURM_ARRAY_TASK_ID}]}
pav_trait=$(echo $pav_mod | grep -P YP[A-Z0-9]+ -o)
pav_feat=$(echo $pav_mod | grep -P '(?<=top_)[0-9]+(?=_models)' -o)
mapfile -t cnv_models < /mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs/cnv_best_fs_models.txt
cnv_mod=${cnv_models[${SLURM_ARRAY_TASK_ID}]}
cnv_trait=$(echo $cnv_mod | grep -P YP[A-Z0-9]+ -o)
cnv_feat=$(echo $cnv_mod | grep -P '(?<=top_)[0-9]+(?=_models)' -o)

conda activate /mnt/home/seguraab/miniconda3/envs/shap

### SNP RF FS models top 20 features (sklearn version 1.2.2)
echo "SNPs RF FS -----------------"
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP/SNP/fs
scratch=/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs
python ${shap}/SHAP_training_saving_interaction.py \
    -df ${data}/geno.csv \
    -df2 ${data}/pheno.csv \
    -y_name ${snp_trait} \
    -sep , -test ${data}/Test.txt \
    -feat ${scratch}/features_${snp_trait}_top_${snp_feat} \
    -model ${scratch}/${snp_mod} \
    -top 20 -interaction n -interaction_score n \
    -save ${snp_trait}_snp_rf_fs_${snp_feat}_top_${snp_feat}

### PAV RF FS models top 20 features (sklearn version 1.2.2)
echo "PAVs RF FS -----------------"
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP/PAV/fs
scratch=/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs
python ${shap}/SHAP_training_saving_interaction.py \
    -df ${data}/ORFs_pres_abs.csv \
    -df2 ${data}/pheno.csv \
    -y_name ${pav_trait} \
    -sep , -test ${data}/Test.txt \
    -feat ${scratch}/features_pav_${pav_trait}_top_${pav_feat} \
    -model ${scratch}/${pav_mod} \
    -top 20 -interaction n -interaction_score n \
    -save ${pav_trait}_pav_rf_fs_${pav_feat}_top_${pav_feat}

### CNV RF FS models top 20 features (sklearn version 1.2.2)
echo "CNVs RF FS -----------------"
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP/CNV/fs
python ${shap}/SHAP_training_saving_interaction.py \
    -df ${data}/ORFs_no_NA.csv \
    -df2 ${data}/pheno.csv \
    -y_name ${cnv_trait} \
    -sep , -test ${data}/Test.txt \
    -feat ${scratch}/features_cnv_${cnv_trait}_top_${cnv_feat} \
    -model ${scratch}/${cnv_mod} \
    -top 20 -interaction n -interaction_score n \
    -save ${cnv_trait}_cnv_rf_fs_${cnv_feat}_top_${cnv_feat}

conda deactivate

scontrol show job $SLURM_JOB_ID
