#!/bin/bash --login
#SBATCH --time=4:00:00
# SBATCH --array=0-4 # for submitting the baseline RF models
# SBATCH --array=0-299 # for submitting RF FS models
# SBATCH --array=0-14 # for submitting the 15 only lit or only non-lit gene models
# SBATCH --array=0-24 # for submitting the missing randomized only lit (147) or only top non-lit gene models (25)
#SBATCH --array=11 #0-14 # submit SHAP interaction calculation
#SBATCH --nodes=1
#SBATCH --ntasks=2
# SBATCH --cpus-per-task=12 # for RF models
# SBATCH --mem-per-cpu=3G # for RF models
#SBATCH --cpus-per-task=4 # for shap interaction
#SBATCH --mem-per-cpu=100G # for shap interaction
# SBATCH --job-name rf_base
# SBATCH --job-name rf_fs
# SBATCH --job-name rf_only_lit_genes
# SBATCH --job-name rf_only_top_non_lit_genes
# SBATCH --job-name shap_interaction
#SBATCH --job-name sum_shap_interaction
#SBATCH --output=%x_%j.out


TRAITS=( YPDCAFEIN40 YPDCAFEIN50 YPDBENOMYL500 YPDCUSO410MM YPDSODIUMMETAARSENITE )
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
shap=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis
pipe=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
rf=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF

# ############################ Baseline RF Models ###########################
# module purge
# conda activate ml-pipeline
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models

# ## SNP baseline using features with max gini only
# python ${pipe}/ML_regression.py -df ${data}/geno.csv \
#     -df2 ${data}/pheno.csv -sep , \
#     -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt \
#     -feat ../Features_max_gini_from_RF_baseline_imp_snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt \
#     -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#     -tag snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_max_gini_from_RF_baseline_imp \
#     -save snp_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_max_gini_from_RF_baseline_imp

# ##PAV baseline using features with max gini only
# python ${pipe}/ML_regression.py -df ${data}/ORFs_pres_abs.csv \
#     -df2 ${data}/pheno.csv -sep , \
#     -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt \
#     -feat ../Features_max_gini_from_RF_baseline_imp_pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt \
#     -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#     -tag pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_max_gini_from_RF_baseline_imp \
#     -save pav_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_max_gini_from_RF_baseline_imp

# ## CNV baseline using features with max gini only
# python ${pipe}/ML_regression.py -df ${data}/ORFs_no_NA.csv \
#     -df2 ${data}/pheno.csv -sep , \
#     -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt \
#     -feat ../Features_max_gini_from_RF_baseline_imp_cnv_${TRAITS[${SLURM_ARRAY_TASK_ID}]}.txt \
#     -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#     -tag cnv_TRAITS[${SLURM_ARRAY_TASK_ID}]}_max_gini_from_RF_baseline_imp \
#     -save cnv_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_max_gini_from_RF_baseline_imp

# conda deactivate

# ################# Feature Selection of Baseline XGBoost Models #################
# ## Apply feature selection to the above models
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs
# module purge
# conda activate ml-pipeline

# # Get feature subsets
# for trait in ${TRAITS[@]}; do
#     ## SNP feature subsets
#     python ${rf}/1_feature_selection.py -m ../models/snp_${trait}_max_gini_from_RF_baseline_imp_models_rep_0.pkl \
#         -i ../models/snp_${trait}_max_gini_from_RF_baseline_imp_imp \
#         -start 0 -stop 2000 -step 100 \
#         -o snp_${trait}_max_gini_from_RF_baseline_imp_features
#     ## PAV feature subsets
#     python ${rf}/1_feature_selection.py -m ../models/pav_${trait}_max_gini_from_RF_baseline_imp_models_rep_0.pkl \
#         -i ../models/pav_${trait}_max_gini_from_RF_baseline_imp_imp \
#         -start 0 -stop 2000 -step 100 \
#         -o pav_${trait}_max_gini_from_RF_baseline_imp_features
#     ## CNV feature subsets
#     python ${rf}/1_feature_selection.py -m ../models/cnv_${trait}_max_gini_from_RF_baseline_imp_models_rep_0.pkl \
#         -i ../models/cnv_${trait}_max_gini_from_RF_baseline_imp_imp \
#         -start 0 -stop 2000 -step 100 \
#         -o cnv_${trait}_max_gini_from_RF_baseline_imp_features
# done

# # Get file with RF FS model commands to run
# for trait in ${TRAITS[@]}; do
#     ## SNP RF FS models
#     feat_snp=($(ls snp_${trait}_max_gini_from_RF_baseline_imp_features_top_*))
#     for feat in ${feat_snp[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${pipe}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -sep , -y_name ${trait} -test ${data}/Test.txt -feat snp_${trait}_max_gini_from_RF_baseline_imp_features_top_${feat_num} -alg RF -n_jobs 12 -n 10 -cv_num 5 -plots f -tag snp_${trait}_max_gini_from_RF_baseline_imp_top_${feat_num} -save snp_${trait}_max_gini_from_RF_baseline_imp_top_${feat_num}" >> RF_FS_runs.txt
#     done

#     ## PAV RF FS models
#     feat_pav=($(ls pav_${trait}_max_gini_from_RF_baseline_imp_features_top_*))
#     for feat in ${feat_pav[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${pipe}/ML_regression.py -df ${data}/ORFs_pres_abs.csv -df2 ${data}/pheno.csv -sep , -y_name ${trait} -test ${data}/Test.txt -feat pav_${trait}_max_gini_from_RF_baseline_imp_features_top_${feat_num} -alg RF -n_jobs 12 -n 10 -cv_num 5 -plots f -tag pav_${trait}_max_gini_from_RF_baseline_imp_top_${feat_num} -save pav_${trait}_max_gini_from_RF_baseline_imp_top_${feat_num}" >> RF_FS_runs.txt
#     done

#     ## CNV RF FS models
#     feat_cnv=($(ls cnv_${trait}_max_gini_from_RF_baseline_imp_features_top_*))
#     for feat in ${feat_cnv[@]}; do
#         feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
#         echo "python ${pipe}/ML_regression.py -df ${data}/ORFs_no_NA.csv -df2 ${data}/pheno.csv -sep , -y_name ${trait} -test ${data}/Test.txt -feat cnv_${trait}_max_gini_from_RF_baseline_imp_features_top_${feat_num} -alg RF -n_jobs 12 -n 10 -cv_num 5 -plots f -tag cnv_${trait}_max_gini_from_RF_baseline_imp_top_${feat_num} -save cnv_${trait}_max_gini_from_RF_baseline_imp_top_${feat_num}" >> RF_FS_runs.txt
#     done
# done

# # Submit commands in XGB_FS_runs.txt as a SLURM JOB ARRAY
# readarray -t joblist < RF_FS_runs.txt
# ${joblist[${SLURM_ARRAY_TASK_ID}]}
# conda deactivate

########### Get the best FS models to combine with literature genes ############
# Table S11 in 0_Table_S10-11_make_benchmark_gene_models.py

# ############## Important & unimportant literature genes RF Models ##############
# ## These models will be used to calculate SHAP Interactions
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/only_lit_genes_models
# mapfile -t input_data < /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs/only_lit_genes_models_to_run.txt # individual lit gene lists
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P [a-z4_-]+_lit_genes -o)

# module purge
# conda activate ml-pipeline

# if [[ $type == "snp" ]]; then
#     python ${pipe}/ML_regression.py -df ${data}/geno.csv \
#         -df2 ${data}/pheno.csv -sep , \
#         -y_name ${trait} \
#         -test ${data}/Test.txt \
#         -feat ../models_fs/${type}_${trait}${litInfo} \
#         -alg RF -n_jobs 12 -n 10 -cv_num 5 -plots f \
#         -tag ${type}_${trait}${litInfo} \
#         -save ${type}_${trait}${litInfo}
# fi
# if [[ $type == "pav" ]]; then
#     python ${pipe}/ML_regression.py -df ${data}/ORFs_pres_abs.csv \
#         -df2 ${data}/pheno.csv -sep , \
#         -y_name ${trait} \
#         -test ${data}/Test.txt \
#         -feat ../models_fs/${type}_${trait}${litInfo} \
#         -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#         -tag ${type}_${trait}${litInfo} \
#         -save ${type}_${trait}${litInfo}
# fi
# if [[ $type == "cnv" ]]; then
#     python ${pipe}/ML_regression.py -df ${data}/ORFs_no_NA.csv \
#         -df2 ${data}/pheno.csv -sep , \
#         -y_name ${trait} \
#         -test ${data}/Test.txt \
#         -feat ../models_fs/${type}_${trait}${litInfo} \
#         -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#         -tag ${type}_${trait}${litInfo} \
#         -save ${type}_${trait}${litInfo}
# fi

# conda deactivate

# ### Randomized Important & unimportant literature genes RF Models
# ## Only the training set was randomized, the test set was just concatenated to the end of pheno_randomized_{i}.csv
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/only_lit_genes_models

# # Get list of missing randomized models
# # python shap_interaction_missing_randomized.py

# mapfile -t missing_runs < ../only_lit_genes_models_missing_runs.txt

# module purge
# conda activate ml-pipeline

# run=${missing_runs[${SLURM_ARRAY_TASK_ID}]}
# run='/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs/cnv_YPDBENOMYL500_benomyl_lit_genes'
# run='/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs/cnv_YPDCUSO410MM_cuso4_lit_genes'
# type=$(echo ${run} | grep -P snp\|pav\|cnv -o)
# trait=$(echo ${run} | grep -P YPD[A-Z0-9]+ -o)
# litInfo=$(echo ${run} | grep -P [a-z4_-]+_lit_genes -o)
# i=$(echo ${run} | grep -P [0-9]+$ -o)

# for i in {0..99}; do
    # if [[ $type == "snp" ]]; then
    #     python ${pipe}/ML_regression.py -df ${data}/geno.csv \
    #         -df2 /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/pheno_randomized_${i}.csv -sep , \
    #         -y_name ${trait} \
    #         -test ${data}/Test.txt \
    #         -feat ../../models_fs/${type}_${trait}${litInfo} \
    #         -alg RF -n_jobs 12 -n 10 -cv_num 5 -plots f \
    #         -tag ${type}_${trait}${litInfo}_randomized_${i} \
    #         -save ${type}_${trait}${litInfo}_randomized_${i}
    # fi
    # if [[ $type == "pav" ]]; then
    #     python ${pipe}/ML_regression.py -df ${data}/ORFs_pres_abs.csv \
    #         -df2 /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/pheno_randomized_${i}.csv -sep , \
    #         -y_name ${trait} \
    #         -test ${data}/Test.txt \
    #         -feat ../../models_fs/${type}_${trait}${litInfo} \
    #         -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
    #         -tag ${type}_${trait}${litInfo}_randomized_${i} \
    #         -save ${type}_${trait}${litInfo}_randomized_${i}
    # fi
#     if [[ $type == "cnv" ]]; then
#         python ${pipe}/ML_regression.py -df ${data}/ORFs_no_NA.csv \
#             -df2 /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/pheno_randomized_${i}.csv -sep , \
#             -y_name ${trait} \
#             -test ${data}/Test.txt \
#             -feat ../../models_fs/${type}_${trait}${litInfo} \
#             -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#             -tag ${type}_${trait}${litInfo}_randomized_${i} \
#             -save ${type}_${trait}${litInfo}_randomized_${i}
#     fi
# done

# conda deactivate


# #################### Only Top Non-Literature Genes RF Models ###################
# ## These models will be used to calculate SHAP Interactions
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/only_top_non_lit_genes_models
# mapfile -t input_data < /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs/only_top_non_lit_genes_models_to_run.txt # individual lit gene lists
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P top_[a-z4_-]+_genes -o)

# module purge
# conda activate ml-pipeline

# if [[ $type == "snp" ]]; then
#     python ${pipe}/ML_regression.py -df ${data}/geno.csv \
#         -df2 ${data}/pheno.csv -sep , \
#         -y_name ${trait} \
#         -test ${data}/Test.txt \
#         -feat ../models_fs/${type}_${trait}_${litInfo} \
#         -alg RF -n_jobs 12 -n 10 -cv_num 5 -plots f \
#         -tag ${type}_${trait}_${litInfo} \
#         -save ${type}_${trait}_${litInfo}
# fi
# if [[ $type == "pav" ]]; then
#     python ${pipe}/ML_regression.py -df ${data}/ORFs_pres_abs.csv \
#         -df2 ${data}/pheno.csv -sep , \
#         -y_name ${trait} \
#         -test ${data}/Test.txt \
#         -feat ../models_fs/${type}_${trait}_${litInfo} \
#         -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#         -tag ${type}_${trait}_${litInfo} \
#         -save ${type}_${trait}_${litInfo}
# fi
# if [[ $type == "cnv" ]]; then
#     python ${pipe}/ML_regression.py -df ${data}/ORFs_no_NA.csv \
#         -df2 ${data}/pheno.csv -sep , \
#         -y_name ${trait} \
#         -test ${data}/Test.txt \
#         -feat ../models_fs/${type}_${trait}_${litInfo} \
#         -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
#         -tag ${type}_${trait}_${litInfo} \
#         -save ${type}_${trait}_${litInfo}
# fi

# ### Randomized only top non-literature genes RF Models
# ## Only the training set was randomized, the test set was just concatenated to the end of pheno_randomized_{i}.csv
# cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/only_top_non_lit_genes_models

# # Get list of missing randomized models
# python shap_interaction_missing_randomized.py

# mapfile -t missing_runs < ../only_top_non_lit_genes_models_missing_runs.txt

# module purge
# conda activate ml-pipeline

# for i in {0..99}; do
# run=${missing_runs[${SLURM_ARRAY_TASK_ID}]}
# type=$(echo ${run} | grep -P snp\|pav\|cnv -o)
# trait=$(echo ${run} | grep -P YPD[A-Z0-9]+ -o)
# litInfo=$(echo ${run} | grep -P top_[a-z4_-]+_genes -o)
# i=$(echo ${run} | grep -P [0-9]+$ -o)

    # if [[ $type == "snp" ]]; then
    #     python ${pipe}/ML_regression.py -df ${data}/geno.csv \
    #         -df2 /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/pheno_randomized_${i}.csv -sep , \
    #         -y_name ${trait} \
    #         -test ${data}/Test.txt \
    #         -feat /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs/${type}_${trait}_${litInfo} \
    #         -alg RF -n_jobs 12 -n 10 -cv_num 5 -plots f \
    #         -tag ${type}_${trait}_${litInfo}_randomized_${i} \
    #         -save ${type}_${trait}_${litInfo}_randomized_${i}
    # fi
    # if [[ $type == "pav" ]]; then
    #     python ${pipe}/ML_regression.py -df ${data}/ORFs_pres_abs.csv \
    #         -df2 /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/pheno_randomized_${i}.csv -sep , \
    #         -y_name ${trait} \
    #         -test ${data}/Test.txt \
    #         -feat /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs/${type}_${trait}_${litInfo} \
    #         -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
    #         -tag ${type}_${trait}_${litInfo}_randomized_${i} \
    #         -save ${type}_${trait}_${litInfo}_randomized_${i}
    # fi
    # if [[ $type == "cnv" ]]; then
    #     python ${pipe}/ML_regression.py -df ${data}/ORFs_no_NA.csv \
    #         -df2 /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/pheno_randomized_${i}.csv -sep , \
    #         -y_name ${trait} \
    #         -test ${data}/Test.txt \
    #         -feat /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs/${type}_${trait}_${litInfo} \
    #         -alg RF -n_jobs 12 -n 20 -cv_num 5 -plots f \
    #         -tag ${type}_${trait}_${litInfo}_randomized_${i} \
    #         -save ${type}_${trait}_${litInfo}_randomized_${i}
    # fi
    # done

# conda deactivate

##### SHAP Values & Interactions for the optimized RF models #######
## Get best RF FS (optimized) models first
# see 0_Table_S10-11_make_benchmark_gene_models.py line 153

## Now calculate SHAP interaction scores
module purge
conda activate shap

cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction

## best optimized RF model (from feature selection)
mapfile -t input_data < RF/models_fs/best_rf_fs_models.txt
file=${input_data[${SLURM_ARRAY_TASK_ID}]}
trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -oP 'YPD[A-Z0-9]+')
type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -oP snp\|pav\|cnv)
litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -oP 'max_[a-zRF0-9_]+_models_rep_[0-9]{1,2}(?=\.pkl)')
featFile=$(echo $litInfo | sed -E 's/imp/imp_features/' | sed -E 's/_models_rep_[0-9]+$//')

# ## only lit gene models (idk if I should calculate SHAP interactions for these, would need to make a shap_best_models.py)
# mapfile -t input_data < RF/models_fs/only_lit_genes_models_to_run.txt # individual lit gene lists
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P [a-z4_-]+_lit_genes -o)

# ## only top non-lit gene models (idk if I should calculate SHAP interactions for these, would need to make a shap_best_models.py)
# mapfile -t input_data < RF/models_fs/only_top_non_lit_genes_models_to_run.txt # individual lit gene lists
# file=${input_data[${SLURM_ARRAY_TASK_ID}]}
# trait=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P YPD[A-Z0-9]+ -o)
# type=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P snp\|pav\|cnv -o)
# litInfo=$(echo ${input_data[${SLURM_ARRAY_TASK_ID}]} | grep -P top_[a-z4_-]+_genes -o)

if [[ $type == "snp" ]]; then
    cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/SNP/best_rf_fs
    # python ${shap}/SHAP_training_saving_interaction.py \
    #     -df ${data}/geno.csv \
    #     -df2 ${data}/pheno.csv -sep , \
    #     -y_name ${trait} \
    #     -test ${data}/Test.txt \
    #     -feat ../../RF/models_fs/${type}_${trait}_${featFile} \
    #     -model ${file} \
    #     -top 10 -interaction n -interaction_score y \
    #     -save ${type}_${trait}_${litInfo}

    echo "Summing SHAP values across isolates for SNP"
    python ${shap}/Sum_SHAP_interaction.py \
        -path /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/SNP/best_rf_fs/ \
        -y ${trait} \
        -dtype snp \
        -save shap_interaction_scores_${type}_${trait}_${litInfo}
fi
if [[ $type == "pav" ]]; then
    ## PAV models
    cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/PAV/best_rf_fs
    # python ${shap}/SHAP_training_saving_interaction.py \
    #     -df ${data}/ORFs_pres_abs.csv \
    #     -df2 ${data}/pheno.csv -sep , \
    #     -y_name ${trait} \
    #     -test ${data}/Test.txt \
    #     -feat ../../RF/models_fs/${type}_${trait}_${featFile} \
    #     -model ${file} \
    #     -top 10 -interaction n -interaction_score y \
    #     -save ${type}_${trait}_${litInfo} #-gpu n

    echo "Summing SHAP values across isolates for PAV"
    python ${shap}/Sum_SHAP_interaction.py \
        -path /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/PAV/best_rf_fs/ \
        -y ${trait} \
        -dtype pav \
        -save shap_interaction_scores_${type}_${trait}_${litInfo}
fi


if [[ $type == "cnv" ]]; then
    ## CNV models
    cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/CNV/best_rf_fs
    # python ${shap}/SHAP_training_saving_interaction.py \
    #     -df ${data}/ORFs_no_NA.csv \
    #     -df2 ${data}/pheno.csv -sep , \
    #     -y_name ${trait} \
    #     -test ${data}/Test.txt \
    #     -feat ../../RF/models_fs/${type}_${trait}_${featFile} \
    #     -model ${file} \
    #     -top 10 -interaction n -interaction_score y \
    #     -save ${type}_${trait}_${litInfo} #-gpu n

    echo "Summing SHAP values across isolates for CNV"
    python ${shap}/Sum_SHAP_interaction.py \
        -path /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/CNV/best_rf_fs/ \
        -y ${trait} \
        -dtype cnv \
        -save shap_interaction_scores_${type}_${trait}_${litInfo}
fi

conda deactivate

scontrol show job $SLURM_JOB_ID

