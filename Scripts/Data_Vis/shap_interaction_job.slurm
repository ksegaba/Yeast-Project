#!/bin/bash --login
#SBATCH --time=6:00:00
# SBATCH --array=0-4 # for submitting the baseline XGB models
# SBATCH --array=0-312 # for submitting FS XGB models 0-312
#SBATCH --array=0-60 # for submitting the 60 fs + lit gene models or the 15 combined lit gene models
# SBATCH --array=0-14 # submit SHAP interaction calculation
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4 # for xgb models
#SBATCH --mem-per-cpu=9G # for xgb models, 6gb is not enough
# SBATCH --cpus-per-task=2 # for shap interaction
# SBATCH --mem-per-cpu=100G # for shap interaction
# SBATCH --account=data-machine
# SBATCH --mem=1TB # use all the data machine memory
# SBATCH --cpus-per-task=2 # num cpus for data machine
# SBATCH --job-name xgb_base
# SBATCH --job-name xgb_fs
# SBATCH --job-name xgb_fs_lit
#SBATCH --job-name xgb_fs_non_lit
# SBATCH --job-name fs_plus_comb
#SBATCH --output=%x_%j.out

TRAITS=( YPDCAFEIN40 YPDCAFEIN50 YPDBENOMYL500 YPDCUSO410MM YPDSODIUMMETAARSENITE )
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
shap=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis
pipe=/mnt/home/seguraab/ML-Pipeline
xgb=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_XGBoost

echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]}"

module purge
conda activate /mnt/home/seguraab/miniconda3/envs/shap
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB

# ############################ Baseline XGBoost Models ###########################
# ## SNP XGB baseline using features with max gini only
# python ${xgb}/xgboost_model.py -X ${data}/geno.csv \
#     -Y ${data}/pheno.csv \
#     -test ${data}/Test.txt \
#     -trait ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -save /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB \
#     -prefix snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -type snp -fold 5 -n 10 \
#     -feat ../RF/Features_max_gini_from_RF_baseline_imp_snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt

# ##PAV XGB baseline using features with max gini only #********* done 3/19/2024
# python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
#     -Y ${data}/pheno.csv \
#     -test ${data}/Test.txt \
#     -trait ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -save /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB \
#     -prefix pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -type pav -fold 5 -n 10 \
#     -feat ../RF/Features_max_gini_from_RF_baseline_imp_pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt

# ## CNV XGB baseline using features with max gini only #********* done 3/19/2024
# python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
#     -Y ${data}/pheno.csv \
#     -test ${data}/Test.txt \
#     -trait ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -save /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB \
#     -prefix cnv_${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -type cnv -fold 5 -n 10 \
#     -feat ../RF/Features_max_gini_from_RF_baseline_imp_cnv_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt
# conda deactivate

################# Feature Selection of Baseline XGBoost Models #################
## Apply feature selection to the above models
# Get feature subsets
# for trait in ${TRAITS[@]}; do
#     ## SNP feature subsets
#     python ${xgb}/feature_selection.py -m snp_${trait}_model_rep_0.pkl \
#         -i snp_${trait}_imp.csv \
#         -start 0 -stop 6500 -step 100 \
#         -o fs/features_snp_${trait}
#     ## PAV feature subsets
#     python ${xgb}/feature_selection.py -m pav_${trait}_model_rep_0.pkl \
#         -i pav_${trait}_imp.csv \
#         -start 0 -stop 6500 -step 100 \
#         -o fs/features_pav_${trait}
#     ## CNV feature subsets
#     python ${xgb}/feature_selection.py -m cnv_${trait}_model_rep_0.pkl \
#         -i cnv_${trait}_imp.csv \
#         -start 0 -stop 6500 -step 100 \
#         -o fs/features_cnv_${trait}
# done

# Get file with xgboost FS model commands to run
# cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs
# for trait in ${TRAITS[@]}; do
#     ## SNP XGB FS models
#     feat_snp=($(ls features_snp_${trait}_top_*))
#     for feat in ${feat_snp[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${xgb}/xgboost_model.py -X ${data}/geno.csv -Y ${data}/pheno.csv -test ${data}/Test.txt -trait ${trait} -save . -prefix snp_${trait}_top_${feat_num} -type snp -fold 5 -n 10 -feat ${feat}" >> XGB_FS_runs.txt
#     done

#     ## PAV XGB FS models
#     feat_pav=($(ls features_pav_${trait}_top_*))
#     for feat in ${feat_pav[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv -Y ${data}/pheno.csv -test ${data}/Test.txt -trait ${trait} -save . -prefix pav_${trait}_top_${feat_num} -type pav -fold 5 -n 10 -feat ${feat}" >> XGB_FS_runs.txt
#     done

#     ## CNV XGB FS models
#     feat_cnv=($(ls features_cnv_${trait}_top_*))
#     for feat in ${feat_cnv[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv -Y ${data}/pheno.csv -test ${data}/Test.txt -trait ${trait} -save . -prefix cnv_${trait}_top_${feat_num} -type cnv -fold 5 -n 10 -feat ${feat}" >> XGB_FS_runs.txt
#     done
# done

# Submit commands in XGB_FS_runs.txt as a SLURM JOB ARRAY
# cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs
# readarray -t joblist < XGB_FS_runs.txt
# ${joblist[${SLURM_ARRAY_TASK_ID}]}
# conda deactivate

########### Get the best FS models to combine with literature genes ############
# Table S11 in Manuscript Tables .py

# ##################### FS + Literature Genes XGBoost Models #####################
# ## Run the new set of models that use the top features and the individual literature gene lists  ##### or the combined literature genes
# ## These models will be used to calculate SHAP Interactions
# cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs
# mapfile -t input_data < fs_plus_lit_genes_models_to_run.txt # individual lit gene lists
# # mapfile -t input_data < fs_plus_comb_lit_genes_models_to_run.txt # combined lit gene lists ##### move to it's own section
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P top_[0-9]+_plus_[a-z0-9_-]+_lit_genes -o)
# if [[ $type == "snp" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/geno.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type snp -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "pav" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type pav -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "cnv" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type cnv -fold 5 -n 10 \
#         -feat $file
# fi
# conda deactivate

# ############# Randomized Label FS + Literature Genes XGBoost Models ############
# ## Only the training set was randomized, the test set was just concatenated to the end of pheno_randomized.csv
# cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/randomized
# mapfile -t input_data < ../fs/fs_plus_comb_lit_genes_models_to_run.txt # combined lit gene lists
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P top_[0-9]+_plus_[a-z0-9_-]+_lit_genes -o)
# for i in {0..100}; do
#     if [[ $type == "snp" ]]; then
#         python ${xgb}/xgboost_model.py -X ${data}/geno.csv \
#             -Y pheno_randomized_${i}.csv \
#             -test ${data}/Test.txt \
#             -trait $trait \
#             -save . \
#             -prefix ${type}_${trait}_${litInfo}_rand_${i} \
#             -type snp -fold 5 -n 10 \
#             -feat $file -plot f
#     fi
#     if [[ $type == "pav" ]]; then
#         python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
#             -Y pheno_randomized_${i}.csv \
#             -test ${data}/Test.txt \
#             -trait $trait \
#             -save . \
#             -prefix ${type}_${trait}_${litInfo}_rand_${i} \
#             -type pav -fold 5 -n 10 \
#             -feat $file -plot f
#     fi
#     if [[ $type == "cnv" ]]; then
#         python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
#             -Y pheno_randomized_${i}.csv \
#             -test ${data}/Test.txt \
#             -trait $trait \
#             -save . \
#             -prefix ${type}_${trait}_${litInfo}_rand_${i} \
#             -type cnv -fold 5 -n 10 \
#             -feat $file -plot f
#     fi
# done
# conda deactivate

################### FS + Non-literature Gene XGBoost Models ####################
## Run the new set of models that use the top features and the individual literature gene lists  ##### or the combined literature genes
## These models will be used to calculate SHAP Interactions
cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs
mapfile -t input_data < fs_plus_unimportant_non_lit_genes_models_to_run.txt # individual gene lists - Top lit + Top non-lit + Unimportant non-lit
file=${input_data[${SLURM_ARRAY_TASK_ID}]}
trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P top_[0-9]+_plus_[a-z0-9_-]+_non_[a-z0-9_-]+_genes -o)
if [[ $type == "snp" ]]; then
    python ${xgb}/xgboost_model.py -X ${data}/geno.csv \
        -Y ${data}/pheno.csv \
        -test ${data}/Test.txt \
        -trait $trait \
        -save . \
        -prefix ${type}_${trait}_${litInfo} \
        -type snp -fold 5 -n 10 \
        -feat $file
fi
if [[ $type == "pav" ]]; then
    python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
        -Y ${data}/pheno.csv \
        -test ${data}/Test.txt \
        -trait $trait \
        -save . \
        -prefix ${type}_${trait}_${litInfo} \
        -type pav -fold 5 -n 10 \
        -feat $file
fi
if [[ $type == "cnv" ]]; then
    python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
        -Y ${data}/pheno.csv \
        -test ${data}/Test.txt \
        -trait $trait \
        -save . \
        -prefix ${type}_${trait}_${litInfo} \
        -type cnv -fold 5 -n 10 \
        -feat $file
fi
conda deactivate

##### SHAP Values & Interactions for FS + Literature Genes XGBoost Models #######
## Get best models first
# python fs/0_get_best_model.py

## Now calculate SHAP interaction scores
# mapfile -t snp_data < /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs/fs_snp_plus_lit_best_models.txt
# mapfile -t pav_data < /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs/fs_pav_plus_lit_best_models.txt
# mapfile -t cnv_data < /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs/fs_cnv_plus_lit_best_models.txt
# mapfile -t comb_models < /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs/fs_plus_comb_lit_best_models.txt # combined lit gene models
# for SLURM_ARRAY_TASK_ID in {0..14}; do
#     file=${comb_models[${SLURM_ARRAY_TASK_ID}]}
#     trait=$(echo ${comb_models[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
#     type=$(echo ${comb_models[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
#     litInfo=$(echo ${comb_models[${SLURM_ARRAY_TASK_ID}]} | grep -P [snpavc]{3}_[A-Z0-9]+_top_[0-9]+_plus_[a-z0-9]+_lit_genes -o)

    ## SNP models
    # if [[ $type == "snp" ]]; then
    #     cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/SNP
    #     python ${shap}/SHAP_training_saving_interaction.py \
    #         -df ${data}/geno.csv \
    #         -df2 ${data}/pheno.csv -sep , \
    #         -y_name ${trait} \
    #         -test ${data}/Test.txt \
    #         -feat ../XGB/fs/${litInfo} \
    #         -model ${file} \
    #         -top 10 -interaction n -interaction_score y \
    #         -save ${litInfo}

    #     echo "Summing SHAP values across isolates for SNP"
    #     python ${shap}/Sum_SHAP_interaction.py \
    #         -path /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/SNP/ \
    #         -y ${trait} \
    #         -dtype snp \
    #         -save shap_interaction_scores_${litInfo}
    # fi
    # if [[ $type == "pav" ]]; then
    #     ## PAV models
    #     cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/PAV
        # python ${shap}/SHAP_training_saving_interaction.py \
        #     -df ${data}/ORFs_pres_abs.csv \
        #     -df2 ${data}/pheno.csv -sep , \
        #     -y_name ${trait} \
        #     -test ${data}/Test.txt \
        #     -feat ../XGB/fs/${litInfo} \
        #     -model ${file} \
        #     -top 10 -interaction n -interaction_score y \
        #     -save ${litInfo}

    #     echo "Summing SHAP values across isolates for PAV"
    #     python ${shap}/Sum_SHAP_interaction.py \
    #         -path /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/PAV/ \
    #         -y ${trait} \
    #         -dtype pav \
    #         -save shap_interaction_scores_${litInfo}
    # fi

    ## CNV models
    # if [[ $type == "cnv" ]]; then
    #     cd /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/CNV
    #     python ${shap}/SHAP_training_saving_interaction.py \
    #         -df ${data}/ORFs_no_NA.csv \
    #         -df2 ${data}/pheno.csv -sep , \
    #         -y_name ${trait} \
    #         -test ${data}/Test.txt \
    #         -feat ../XGB/fs/${litInfo} \
    #         -model ${file} \
    #         -top 10 -interaction n -interaction_score y \
    #         -save ${litInfo}

    #     # echo "Summing SHAP values across isolates for CNV"
    #     python ${shap}/Sum_SHAP_interaction.py \
    #         -path /mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/CNV/ \
    #         -y ${trait} \
    #         -dtype cnv \
    #         -save shap_interaction_scores_${litInfo}
    # fi
# done
# conda deactivate

scontrol show job $SLURM_JOB_ID

