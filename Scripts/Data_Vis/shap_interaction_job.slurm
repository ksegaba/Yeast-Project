#!/bin/bash --login
#SBATCH --time=1:00:00
#SBATCH --array=0 # for submitting the baseline XGB models
# SBATCH --array=0-296 # for submitting FS XGB models 0-312
# SBATCH --array=0-14 # for submitting the 15 fs + lit gene models or the 15 combined lit gene models
# SBATCH --array=0-14 # submit SHAP interaction calculation
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4 # for xgb models
#SBATCH --mem-per-cpu=9G # for xgb models, 6gb is not enough
# SBATCH --cpus-per-task=2 # for shap interaction
# SBATCH --mem-per-cpu=100G # for shap interaction
#SBATCH --job-name xgb_base
# SBATCH --job-name xgb_fs
# SBATCH --job-name xgb_fs_lit_or_non_lit_or_comb_lit
# SBATCH --job-name xgb_fs_lit_or_non_lit_or_comb_lit_rand
# SBATCH --job-name shap_interaction
#SBATCH --output=%x_%j.out

TRAITS=( YPDCAFEIN40 YPDCAFEIN50 YPDBENOMYL500 YPDCUSO410MM YPDSODIUMMETAARSENITE )
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
shap=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis
pipe=/mnt/home/seguraab/ML-Pipeline
xgb=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_XGBoost

# echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]}"

module purge
conda activate /mnt/home/seguraab/miniconda3/envs/shap
cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction

############################ Baseline XGBoost Models ###########################
cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB
## SNP XGB baseline using features with max gini only
python ${xgb}/1b_xgb_regression.py -X ${data}/geno.csv \
    -Y ${data}/pheno.csv \
    -test ${data}/Test.txt \
    -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -save /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB \
    -prefix snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -tag snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]} -fold 5 -n 20 -plot f \
    -feat RF/Features_max_gini_from_RF_baseline_imp_snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt

##PAV XGB baseline using features with max gini only
python ${xgb}/1b_xgb_regression.py -X ${data}/ORFs_pres_abs.csv \
    -Y ${data}/pheno.csv \
    -test ${data}/Test.txt \
    -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -save /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB \
    -prefix pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -tag pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]} -fold 5 -n 20 -plot f \
    -feat RF/Features_max_gini_from_RF_baseline_imp_pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt

## CNV XGB baseline using features with max gini only
python ${xgb}/1b_xgb_regression.py -X ${data}/ORFs_no_NA.csv \
    -Y ${data}/pheno.csv \
    -test ${data}/Test.txt \
    -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -save /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB \
    -prefix cnv_${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -tag cnv_${TRAITS[${SLURM_ARRAY_TASK_ID}]} -fold 5 -n 20 -plot f \
    -feat RF/Features_max_gini_from_RF_baseline_imp_cnv_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt

conda deactivate

# ################# Feature Selection of Baseline XGBoost Models #################
# ## Apply feature selection to the above models
# # Get feature subsets
# for trait in ${TRAITS[@]}; do
#     ## SNP feature subsets
#     python ${xgb}/feature_selection.py -m snp_${trait}_model_rep_0.pkl \
#         -i snp_${trait}_imp.csv \
#         -start 0 -stop 6000 -step 200 \
#         -o fs/features_snp_${trait}
#     ## PAV feature subsets
#     python ${xgb}/feature_selection.py -m pav_${trait}_model_rep_0.pkl \
#         -i pav_${trait}_imp.csv \
#         -start 0 -stop 6000 -step 200 \
#         -o fs/features_pav_${trait}
#     ## CNV feature subsets
#     python ${xgb}/feature_selection.py -m cnv_${trait}_model_rep_0.pkl \
#         -i cnv_${trait}_imp.csv \
#         -start 0 -stop 6000 -step 200 \
#         -o fs/features_cnv_${trait}
# done

# # Get file with xgboost FS model commands to run
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB/fs
# for trait in ${TRAITS[@]}; do
#     ## SNP XGB FS models
#     feat_snp=($(ls features_snp_${trait}_top_*))
#     for feat in ${feat_snp[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${xgb}/1b_xgb_regression.py -X ${data}/geno.csv -Y ${data}/pheno.csv -test ${data}/Test.txt -y_name ${trait} -save . -prefix snp_${trait}_top_${feat_num} -tag snp_${trait}_top_${feat_num} -fold 5 -n 10 -feat ${feat}" >> XGB_FS_runs.txt
#     done

#     ## PAV XGB FS models
#     feat_pav=($(ls features_pav_${trait}_top_*))
#     for feat in ${feat_pav[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${xgb}/1b_xgb_regression.py -X ${data}/ORFs_pres_abs.csv -Y ${data}/pheno.csv -test ${data}/Test.txt -y_name ${trait} -save . -prefix pav_${trait}_top_${feat_num} -tag pav_${trait}_top_${feat_num} -fold 5 -n 10 -feat ${feat}" >> XGB_FS_runs.txt
#     done

#     ## CNV XGB FS models
#     feat_cnv=($(ls features_cnv_${trait}_top_*))
#     for feat in ${feat_cnv[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${xgb}/1b_xgb_regression.py -X ${data}/ORFs_no_NA.csv -Y ${data}/pheno.csv -test ${data}/Test.txt -y_name ${trait} -save . -prefix cnv_${trait}_top_${feat_num} -tag cnv_${trait}_top_${feat_num} -fold 5 -n 10 -feat ${feat}" >> XGB_FS_runs.txt
#     done
# done

# # Submit commands in XGB_FS_runs.txt as a SLURM JOB ARRAY
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB/fs
# readarray -t joblist < XGB_FS_runs.txt
# ${joblist[${SLURM_ARRAY_TASK_ID}]}
# conda deactivate

########### Get the best FS models to combine with literature genes ############
# Table S11 in Manuscript Tables .py

# ##################### FS + Literature Genes XGBoost Models #####################
# ## Run the new set of models that use the top features and the individual literature gene lists  ##### or the combined literature genes
# ## These models will be used to calculate SHAP Interactions
# ## 
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB/fs
# mapfile -t input_data < fs_plus_lit_genes_models_to_run.txt # individual lit gene lists
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
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "pav" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "cnv" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi

# ############# Randomized Label FS + Literature Genes XGBoost Models ############
# ## Only the training set was randomized, the test set was just concatenated to the end of pheno_randomized_{i}.csv
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB/randomized
# mapfile -t input_data < ../fs/fs_plus_lit_genes_models_to_run.txt # individual lit gene lists
# # mapfile -t input_data < fs_plus_unimportant_non_lit_genes_models_to_run.txt
# # mapfile -t input_data < ../fs/fs_plus_comb_lit_genes_models_to_run.txt # combined lit gene lists
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
#             -type ${type}_${trait}_${litInfo}_rand_${i} -fold 5 -n 10 \
#             -feat $file -plot f
#     fi
#     if [[ $type == "pav" ]]; then
#         python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
#             -Y pheno_randomized_${i}.csv \
#             -test ${data}/Test.txt \
#             -trait $trait \
#             -save . \
#             -prefix ${type}_${trait}_${litInfo}_rand_${i} \
#             -type ${type}_${trait}_${litInfo}_rand_${i} -fold 5 -n 10 \
#             -feat $file -plot f
#     fi
#     if [[ $type == "cnv" ]]; then
#         python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
#             -Y pheno_randomized_${i}.csv \
#             -test ${data}/Test.txt \
#             -trait $trait \
#             -save . \
#             -prefix ${type}_${trait}_${litInfo}_rand_${i} \
#             -type ${type}_${trait}_${litInfo}_rand_${i} -fold 5 -n 10 \
#             -feat $file -plot f
#     fi
# done

# ################### FS + Non-literature Gene XGBoost Models ####################
# ## These models will be used to calculate SHAP Interactions
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB/fs
# mapfile -t input_data < fs_plus_unimportant_non_lit_genes_models_to_run.txt # individual gene lists - Top lit + Top non-lit + Unimportant non-lit
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P top_[0-9]+_plus_[a-z0-9_-]+_non_[a-z0-9_-]+_genes -o)
# if [[ $type == "snp" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/geno.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "pav" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "cnv" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi

# ################# FS + Combined Literature Genes XGBoost Models ################
# ## Run the new set of models that use the top features and the combined literature gene list
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/XGB/fs
# mapfile -t input_data < fs_plus_comb_lit_genes_models_to_run.txt # combined lit gene lists
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P top_[0-9]+_plus_comb_lit_genes -o)
# if [[ $type == "snp" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/geno.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "pav" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_pres_abs.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi
# if [[ $type == "cnv" ]]; then
#     python ${xgb}/xgboost_model.py -X ${data}/ORFs_no_NA.csv \
#         -Y ${data}/pheno.csv \
#         -test ${data}/Test.txt \
#         -trait $trait \
#         -save . \
#         -prefix ${type}_${trait}_${litInfo} \
#         -type ${type}_${trait}_${litInfo} -fold 5 -n 10 \
#         -feat $file
# fi
# conda deactivate

# ##### SHAP Values & Interactions for FS + Literature Genes XGBoost Models #######
# ## Get best models first
# # python fs/0_get_best_model.py

# ## Now calculate SHAP interaction scores
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/
# mapfile -t lit_models < XGB/fs/fs_plus_lit_genes_best_models.txt
# mapfile -t non_lit_models < XGB/fs/fs_plus_unimportant_non_lit_genes_best_models.txt
# mapfile -t comb_models < XGB/fs/fs_plus_comb_lit_best_models.txt # combined lit gene models
# for snp_YPDCAFEIN40_top_1400_plus_caffeine_lit_genes_model_rep_3.pkl in {0..14}; do
# # for SLURM_ARRAY_TASK_ID in {3,6,9}; do
#     # file=${lit_models[${SLURM_ARRAY_TASK_ID}]}
#     # trait=$(echo ${lit_models[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
#     # type=$(echo ${lit_models[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
#     # litInfo=$(echo ${lit_models[${SLURM_ARRAY_TASK_ID}]} | grep -P [snpavc]{3}_[A-Z0-9]+_top_[0-9]+_plus_[a-z0-9-_]+_lit_genes -o)

#     # file=${non_lit_models[${SLURM_ARRAY_TASK_ID}]}
#     # trait=$(echo ${non_lit_models[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
#     # type=$(echo ${non_lit_models[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
#     # litInfo=$(echo ${non_lit_models[${SLURM_ARRAY_TASK_ID}]} | grep -P [snpavc]{3}_[A-Z0-9]+_top_[0-9]+_plus_[a-z0-9-_]+_genes -o)

#     file=${comb_models[${SLURM_ARRAY_TASK_ID}]}
#     trait=$(echo ${comb_models[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
#     type=$(echo ${comb_models[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
#     litInfo=$(echo ${comb_models[${SLURM_ARRAY_TASK_ID}]} | grep -P [snpavc]{3}_[A-Z0-9]+_top_[0-9]+_plus_[a-z0-9]+_lit_genes -o)

#     # SNP models
#     if [[ $type == "snp" ]]; then
#         cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/SNP
#         python ${shap}/SHAP_training_saving_interaction.py \
#             -df ${data}/geno.csv \
#             -df2 ${data}/pheno.csv -sep , \
#             -y_name ${trait} \
#             -test ${data}/Test.txt \
#             -feat ../XGB/fs/${litInfo} \
#             -model ../XGB/fs/${file} \
#             -top 10 -interaction n -interaction_score y \
#             -save ${litInfo}

#         echo "Summing SHAP values across isolates for SNP"
#         python ${shap}/Sum_SHAP_interaction.py \
#             -path /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/SNP/ \
#             -y ${trait} \
#             -dtype snp \
#             -save shap_interaction_scores_${litInfo}
#     fi
#     if [[ $type == "pav" ]]; then
#         ## PAV models
#         cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/PAV
#         python ${shap}/SHAP_training_saving_interaction.py \
#             -df ${data}/ORFs_pres_abs.csv \
#             -df2 ${data}/pheno.csv -sep , \
#             -y_name ${trait} \
#             -test ${data}/Test.txt \
#             -feat ../XGB/fs/${litInfo} \
#             -model ../XGB/fs/${file} \
#             -top 10 -interaction n -interaction_score y \
#             -save ${litInfo}

#         echo "Summing SHAP values across isolates for PAV"
#         python ${shap}/Sum_SHAP_interaction.py \
#             -path /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/PAV/ \
#             -y ${trait} \
#             -dtype pav \
#             -save shap_interaction_scores_${litInfo}
#     fi

#     # CNV models
#     if [[ $type == "cnv" ]]; then
#         cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/CNV
#         python ${shap}/SHAP_training_saving_interaction.py \
#             -df ${data}/ORFs_no_NA.csv \
#             -df2 ${data}/pheno.csv -sep , \
#             -y_name ${trait} \
#             -test ${data}/Test.txt \
#             -feat ../XGB/fs/${litInfo} \
#             -model ../XGB/fs/${file} \
#             -top 10 -interaction n -interaction_score y \
#             -save ${litInfo}

#         # echo "Summing SHAP values across isolates for CNV"
#         python ${shap}/Sum_SHAP_interaction.py \
#             -path /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/CNV/ \
#             -y ${trait} \
#             -dtype cnv \
#             -save shap_interaction_scores_${litInfo}
#     fi
# done
# conda deactivate

scontrol show job $SLURM_JOB_ID

