#!/bin/bash --login 
#SBATCH --array=0-34 # 35 array tasks
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name shap_values
#SBATCH --output=%x_%j.out

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4

cd /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis

scratch=/mnt/scratch/seguraab/yeast_project/SNP_yeast_RF_results
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
TRAITS=( YPACETATE YPD6AU YPD14 YPD40 YPD42 YPDANISO10 YPDANISO20 YPDANISO50 YPDBENOMYL200 YPDBENOMYL500 YPDCAFEIN40 YPDCAFEIN50 YPDCHX1 YPDCHX05 YPDCUSO410MM YPDDMSO YPDETOH YPDFLUCONAZOLE YPDFORMAMIDE4 YPDFORMAMIDE5 YPDHU YPDKCL2M YPDLICL250MM YPDMV YPDNACL1M YPDNACL15M YPDNYSTATIN YPDSDS YPDSODIUMMETAARSENITE YPETHANOL YPGALACTOSE YPGLYCEROL YPRIBOSE YPSORBITOL YPXYLOSE )
FEATS=( 16 1024 256 1024 512 128 1024 256 512 1024 512 1024 512 1024 1024 256 32 256 256 128 128 256 128 256 1024 512 16 1024 256 1024 128 32 128 256 256 )

echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]} ; ${FEATS[${SLURM_ARRAY_TASK_ID}]}"

### SNP RF FS models top ~20 features (some a little less or a little more)
python SHAP_training.py \
    -df ${data}/geno.csv -sep , \
    -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -test ${data}/Test.txt \
    -feat ${scratch}/feat_exp_rf_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_top_${FEATS[${SLURM_ARRAY_TASK_ID}]} \
    -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_exp_rf_${FEATS[${SLURM_ARRAY_TASK_ID}]}_models.pkl \
    -top 20 -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_20

### SNP RF FS models all features in figure
python SHAP_training.py \
    -df ${data}/geno.csv -sep , \
    -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
    -test ${data}/Test.txt \
    -feat ${scratch}/feat_exp_rf_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_top_${FEATS[${SLURM_ARRAY_TASK_ID}]} \
    -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_exp_rf_${FEATS[${SLURM_ARRAY_TASK_ID}]}_models.pkl \
    -top ${FEATS[${SLURM_ARRAY_TASK_ID}]} \
    -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_${FEATS[${SLURM_ARRAY_TASK_ID}]}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# scratch=/mnt/scratch/seguraab/yeast_project/ORF_yeast_RF_results
# data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
# TRAITS=( YPACETATE YPD6AU YPD14 YPD40 YPD42 YPDANISO10 YPDANISO20 YPDANISO50 YPDBENOMYL200 YPDBENOMYL500 YPDCAFEIN40 YPDCAFEIN50 YPDCHX1 YPDCHX05 YPDCUSO410MM YPDDMSO YPDETOH YPDFLUCONAZOLE YPDFORMAMIDE4 YPDFORMAMIDE5 YPDHU YPDKCL2M YPDLICL250MM YPDMV YPDNACL1M YPDNACL15M YPDNYSTATIN YPDSDS YPDSODIUMMETAARSENITE YPETHANOL YPGALACTOSE YPGLYCEROL YPRIBOSE YPSORBITOL YPXYLOSE )
# ORF_FEATS=( 250 1250 250 500 250 500 500 250 750 250 500 500 250 250 500 500 250 250 250 250 250 250 750 250 250 250 500 250 250 250 500 250 250 250 250 )
# CNO_FEATS=( 250 500 250 250 250 500 250 250 500 250 250 250 500 250 500 500 250 250 250 250 250 250 250 250 250 250 250 250 250 250 250 250 500 250 250 )

# echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]}"

### ORF RF baseline
# python SHAP_training.py -df ${data}/ORFs_pres_abs.csv -sep , \
#     -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt \
#     -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf_baseline_models.pkl \
<<<<<<< HEAD
<<<<<<< HEAD
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf_baseline
=======
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
=======
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b

### ORF RF FS
# python SHAP_training.py -df ${data}/ORFs_pres_abs.csv -sep , \
#     -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt \
#     -feat ${scratch}/feat_rf_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_top_${ORF_FEATS[${SLURM_ARRAY_TASK_ID}]} \
#     -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf_${ORF_FEATS[${SLURM_ARRAY_TASK_ID}]}_models.pkl \
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_orf_${ORF_FEATS[${SLURM_ARRAY_TASK_ID}]}

### CNO RF baseline
# python SHAP_training.py -df ${data}/ORFs_no_NA.csv -sep , \
#     -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cno_baseline_models.pkl \
<<<<<<< HEAD
<<<<<<< HEAD
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cno_baselines
=======
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cno
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
=======
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cno
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b

### CNO RF FS
# python SHAP_training.py -df ${data}/ORFs_no_NA.csv -sep , \
#     -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} \
#     -test ${data}/Test.txt \
#     -feat ${scratch}/feat_rf_${TRAITS[${SLURM_ARRAY_TASK_ID}]}_top_${CNO_FEATS[${SLURM_ARRAY_TASK_ID}]} \
#     -model ${scratch}/${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cno_${CNO_FEATS[${SLURM_ARRAY_TASK_ID}]}_models.pkl \
#     -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_cno_${CNO_FEATS[${SLURM_ARRAY_TASK_ID}]}

scontrol show job $SLURM_JOB_ID