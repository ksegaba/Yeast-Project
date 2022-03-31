#!/bin/bash --login 
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name shap_values
#SBATCH --output=%x_%j

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4

cd /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis
scratch=/mnt/scratch/seguraab/yeast_project/yeast_rf_results
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018

# top ~20 features (some a little less or a little more)
# 03/28/2022: I opened multiple terminals and submitted via the command line (not slurm)
python SHAP_training.py -df ${scratch}/geno_rf_YPACETATE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPACETATE_top_16 -model ${scratch}/YPACETATE_exp_rf_16_models.pkl -top 16 -save YPACETATE
python SHAP_training.py -df ${scratch}/geno_rf_YPD6AU.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD6AU_top_1024 -model ${scratch}/YPD6AU_exp_rf_1024_models.pkl -top 20 -save YPD6AU
python SHAP_training.py -df ${scratch}/geno_rf_YPD14.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD14_top_256 -model ${scratch}/YPD14_exp_rf_256_models.pkl -top 20 -save YPD14
python SHAP_training.py -df ${scratch}/geno_rf_YPD40.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD40_top_1024 -model ${scratch}/YPD40_exp_rf_1024_models.pkl -top 20 -save YPD40
python SHAP_training.py -df ${scratch}/geno_rf_YPD42.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD42_top_512 -model ${scratch}/YPD42_exp_rf_512_models.pkl -top 20 -save YPD42
python SHAP_training.py -df ${scratch}/geno_rf_YPDANISO10.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDANISO10_top_128 -model ${scratch}/YPDANISO10_exp_rf_128_models.pkl -top 20 -save YPDANISO10
python SHAP_training.py -df ${scratch}/geno_rf_YPDANISO20.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDANISO20_top_1024 -model ${scratch}/YPDANISO20_exp_rf_1024_models.pkl -top 20 -save YPDANISO20
python SHAP_training.py -df ${scratch}/geno_rf_YPDANISO50.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDANISO50_top_256 -model ${scratch}/YPDANISO50_exp_rf_256_models.pkl -top 20 -save YPDANISO50
python SHAP_training.py -df ${scratch}/geno_rf_YPDBENOMYL200.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDBENOMYL200_top_512 -model ${scratch}/YPDBENOMYL200_exp_rf_512_models.pkl -top 20 -save YPDBENOMYL200
python SHAP_training.py -df ${scratch}/geno_rf_YPDBENOMYL500.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDBENOMYL500_top_1024 -model ${scratch}/YPDBENOMYL500_exp_rf_1024_models.pkl -top 20 -save YPDBENOMYL500

python SHAP_training.py -df ${scratch}/geno_rf_YPDCAFEIN40.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCAFEIN40_top_512 -model ${scratch}/YPDCAFEIN40_exp_rf_512_models.pkl -top 20 -save YPDCAFEIN40
python SHAP_training.py -df ${scratch}/geno_rf_YPDCAFEIN50.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCAFEIN50_top_1024 -model ${scratch}/YPDCAFEIN50_exp_rf_1024_models.pkl -top 20 -save YPDCAFEIN50
python SHAP_training.py -df ${scratch}/geno_rf_YPDCHX1.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCHX1_top_512 -model ${scratch}/YPDCHX1_exp_rf_512_models.pkl -top 20 -save YPDCHX1
python SHAP_training.py -df ${scratch}/geno_rf_YPDCHX05.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCHX05_top_1024 -model ${scratch}/YPDCHX05_exp_rf_1024_models.pkl -top 20 -save YPDCHX05
python SHAP_training.py -df ${scratch}/geno_rf_YPDCUSO410MM.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCUSO410MM_top_1024 -model ${scratch}/YPDCUSO410MM_exp_rf_1024_models.pkl -top 20 -save YPDCUSO410MM
python SHAP_training.py -df ${scratch}/geno_rf_YPDDMSO.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDDMSO_top_256 -model ${scratch}/YPDDMSO_exp_rf_256_models.pkl -top 20 -save YPDDMSO
python SHAP_training.py -df ${scratch}/geno_rf_YPDETOH.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDETOH_top_32 -model ${scratch}/YPDETOH_exp_rf_32_models.pkl -top 32 -save YPDETOH
python SHAP_training.py -df ${scratch}/geno_rf_YPDFLUCONAZOLE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDFLUCONAZOLE_top_256 -model ${scratch}/YPDFLUCONAZOLE_exp_rf_256_models.pkl -top 20 -save YPDFLUCONAZOLE
python SHAP_training.py -df ${scratch}/geno_rf_YPDFORMAMIDE4.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDFORMAMIDE4_top_256 -model ${scratch}/YPDFORMAMIDE4_exp_rf_256_models.pkl -top 20 -save YPDFORMAMIDE4
python SHAP_training.py -df ${scratch}/geno_rf_YPDFORMAMIDE5.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDFORMAMIDE5_top_128 -model ${scratch}/YPDFORMAMIDE5_exp_rf_128_models.pkl -top 20 -save YPDFORMAMIDE5

python SHAP_training.py -df ${scratch}/geno_rf_YPDHU.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDHU_top_128 -model ${scratch}/YPDHU_exp_rf_128_models.pkl -top 20 -save YPDHU
python SHAP_training.py -df ${scratch}/geno_rf_YPDKCL2M.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDKCL2M_top_256 -model ${scratch}/YPDKCL2M_exp_rf_256_models.pkl -top 20 -save YPDKCL2M
python SHAP_training.py -df ${scratch}/geno_rf_YPDLICL250MM.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDLICL250MM_top_128 -model ${scratch}/YPDLICL250MM_exp_rf_128_models.pkl -top 20 -save YPDLICL250MM
python SHAP_training.py -df ${scratch}/geno_rf_YPDMV.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDMV_top_256 -model ${scratch}/YPDMV_exp_rf_256_models.pkl -top 20 -save YPDMV
python SHAP_training.py -df ${scratch}/geno_rf_YPDNACL1M.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDNACL1M_top_1024 -model ${scratch}/YPDNACL1M_exp_rf_1024_models.pkl -top 20 -save YPDNACL1M
python SHAP_training.py -df ${scratch}/geno_rf_YPDNACL15M.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDNACL15M_top_512 -model ${scratch}/YPDNACL15M_exp_rf_512_models.pkl -top 20 -save YPDNACL15M
python SHAP_training.py -df ${scratch}/geno_rf_YPDNYSTATIN.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDNYSTATIN_top_16 -model ${scratch}/YPDNYSTATIN_exp_rf_16_models.pkl -top 16 -save YPDNYSTATIN
python SHAP_training.py -df ${scratch}/geno_rf_YPDSDS.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDSDS_top_1024 -model ${scratch}/YPDSDS_exp_rf_1024_models.pkl -top 20 -save YPDSDS
python SHAP_training.py -df ${scratch}/geno_rf_YPDSODIUMMETAARSENITE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDSODIUMMETAARSENITE_top_256 -model ${scratch}/YPDSODIUMMETAARSENITE_exp_rf_256_models.pkl -top 20 -save YPDSODIUMMETAARSENITE
python SHAP_training.py -df ${scratch}/geno_rf_YPETHANOL.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPETHANOL_top_1024 -model ${scratch}/YPETHANOL_exp_rf_1024_models.pkl -top 20 -save YPETHANOL

python SHAP_training.py -df ${scratch}/geno_rf_YPGALACTOSE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPGALACTOSE_top_128 -model ${scratch}/YPGALACTOSE_exp_rf_128_models.pkl -top 20 -save YPGALACTOSE
python SHAP_training.py -df ${scratch}/geno_rf_YPGLYCEROL.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPGLYCEROL_top_32 -model ${scratch}/YPGLYCEROL_exp_rf_32_models.pkl -top 32 -save YPGLYCEROL
python SHAP_training.py -df ${scratch}/geno_rf_YPRIBOSE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPRIBOSE_top_128 -model ${scratch}/YPRIBOSE_exp_rf_128_models.pkl -top 20 -save YPRIBOSE
python SHAP_training.py -df ${scratch}/geno_rf_YPSORBITOL.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPSORBITOL_top_256 -model ${scratch}/YPSORBITOL_exp_rf_256_models.pkl -top 20 -save YPSORBITOL
python SHAP_training.py -df ${scratch}/geno_rf_YPXYLOSE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPXYLOSE_top_256 -model ${scratch}/YPXYLOSE_exp_rf_256_models.pkl -top 20 -save YPXYLOSE

# SHAP figures with all features
python SHAP_training.py -df ${scratch}/geno_rf_YPD6AU.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD6AU_top_1024 -model ${scratch}/YPD6AU_exp_rf_1024_models.pkl -top 1024 -save YPD6AU
python SHAP_training.py -df ${scratch}/geno_rf_YPD14.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD14_top_256 -model ${scratch}/YPD14_exp_rf_256_models.pkl -top 256 -save YPD14
python SHAP_training.py -df ${scratch}/geno_rf_YPD40.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD40_top_1024 -model ${scratch}/YPD40_exp_rf_1024_models.pkl -top 1024 -save YPD40
python SHAP_training.py -df ${scratch}/geno_rf_YPD42.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPD42_top_512 -model ${scratch}/YPD42_exp_rf_512_models.pkl -top 512 -save YPD42
python SHAP_training.py -df ${scratch}/geno_rf_YPDANISO10.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDANISO10_top_128 -model ${scratch}/YPDANISO10_exp_rf_128_models.pkl -top 128 -save YPDANISO10
python SHAP_training.py -df ${scratch}/geno_rf_YPDANISO20.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDANISO20_top_1024 -model ${scratch}/YPDANISO20_exp_rf_1024_models.pkl -top 1024 -save YPDANISO20
python SHAP_training.py -df ${scratch}/geno_rf_YPDANISO50.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDANISO50_top_256 -model ${scratch}/YPDANISO50_exp_rf_256_models.pkl -top 256 -save YPDANISO50
python SHAP_training.py -df ${scratch}/geno_rf_YPDBENOMYL200.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDBENOMYL200_top_512 -model ${scratch}/YPDBENOMYL200_exp_rf_512_models.pkl -top 512 -save YPDBENOMYL200
python SHAP_training.py -df ${scratch}/geno_rf_YPDBENOMYL500.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDBENOMYL500_top_1024 -model ${scratch}/YPDBENOMYL500_exp_rf_1024_models.pkl -top 1024 -save YPDBENOMYL500

python SHAP_training.py -df ${scratch}/geno_rf_YPDCAFEIN40.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCAFEIN40_top_512 -model ${scratch}/YPDCAFEIN40_exp_rf_512_models.pkl -top 512 -save YPDCAFEIN40
python SHAP_training.py -df ${scratch}/geno_rf_YPDCAFEIN50.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCAFEIN50_top_1024 -model ${scratch}/YPDCAFEIN50_exp_rf_1024_models.pkl -top 1024 -save YPDCAFEIN50
python SHAP_training.py -df ${scratch}/geno_rf_YPDCHX1.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCHX1_top_512 -model ${scratch}/YPDCHX1_exp_rf_512_models.pkl -top 512 -save YPDCHX1
python SHAP_training.py -df ${scratch}/geno_rf_YPDCHX05.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCHX05_top_1024 -model ${scratch}/YPDCHX05_exp_rf_1024_models.pkl -top 1024 -save YPDCHX05
python SHAP_training.py -df ${scratch}/geno_rf_YPDCUSO410MM.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDCUSO410MM_top_1024 -model ${scratch}/YPDCUSO410MM_exp_rf_1024_models.pkl -top 1024 -save YPDCUSO410MM
python SHAP_training.py -df ${scratch}/geno_rf_YPDDMSO.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDDMSO_top_256 -model ${scratch}/YPDDMSO_exp_rf_256_models.pkl -top 256 -save YPDDMSO
python SHAP_training.py -df ${scratch}/geno_rf_YPDFLUCONAZOLE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDFLUCONAZOLE_top_256 -model ${scratch}/YPDFLUCONAZOLE_exp_rf_256_models.pkl -top 256 -save YPDFLUCONAZOLE
python SHAP_training.py -df ${scratch}/geno_rf_YPDFORMAMIDE4.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDFORMAMIDE4_top_256 -model ${scratch}/YPDFORMAMIDE4_exp_rf_256_models.pkl -top 256 -save YPDFORMAMIDE4
python SHAP_training.py -df ${scratch}/geno_rf_YPDFORMAMIDE5.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDFORMAMIDE5_top_128 -model ${scratch}/YPDFORMAMIDE5_exp_rf_128_models.pkl -top 128 -save YPDFORMAMIDE5

python SHAP_training.py -df ${scratch}/geno_rf_YPDHU.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDHU_top_128 -model ${scratch}/YPDHU_exp_rf_128_models.pkl -top 128 -save YPDHU
python SHAP_training.py -df ${scratch}/geno_rf_YPDKCL2M.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDKCL2M_top_256 -model ${scratch}/YPDKCL2M_exp_rf_256_models.pkl -top 256 -save YPDKCL2M
python SHAP_training.py -df ${scratch}/geno_rf_YPDLICL250MM.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDLICL250MM_top_128 -model ${scratch}/YPDLICL250MM_exp_rf_128_models.pkl -top 128 -save YPDLICL250MM
python SHAP_training.py -df ${scratch}/geno_rf_YPDMV.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDMV_top_256 -model ${scratch}/YPDMV_exp_rf_256_models.pkl -top 256 -save YPDMV
python SHAP_training.py -df ${scratch}/geno_rf_YPDNACL1M.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDNACL1M_top_1024 -model ${scratch}/YPDNACL1M_exp_rf_1024_models.pkl -top 1024 -save YPDNACL1M
python SHAP_training.py -df ${scratch}/geno_rf_YPDNACL15M.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDNACL15M_top_512 -model ${scratch}/YPDNACL15M_exp_rf_512_models.pkl -top 512 -save YPDNACL15M
python SHAP_training.py -df ${scratch}/geno_rf_YPDSDS.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDSDS_top_1024 -model ${scratch}/YPDSDS_exp_rf_1024_models.pkl -top 1024 -save YPDSDS
python SHAP_training.py -df ${scratch}/geno_rf_YPDSODIUMMETAARSENITE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPDSODIUMMETAARSENITE_top_256 -model ${scratch}/YPDSODIUMMETAARSENITE_exp_rf_256_models.pkl -top 256 -save YPDSODIUMMETAARSENITE
python SHAP_training.py -df ${scratch}/geno_rf_YPETHANOL.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPETHANOL_top_1024 -model ${scratch}/YPETHANOL_exp_rf_1024_models.pkl -top 1024 -save YPETHANOL

python SHAP_training.py -df ${scratch}/geno_rf_YPGALACTOSE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPGALACTOSE_top_128 -model ${scratch}/YPGALACTOSE_exp_rf_128_models.pkl -top 128 -save YPGALACTOSE
python SHAP_training.py -df ${scratch}/geno_rf_YPRIBOSE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPRIBOSE_top_128 -model ${scratch}/YPRIBOSE_exp_rf_128_models.pkl -top 128 -save YPRIBOSE
python SHAP_training.py -df ${scratch}/geno_rf_YPSORBITOL.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPSORBITOL_top_256 -model ${scratch}/YPSORBITOL_exp_rf_256_models.pkl -top 256 -save YPSORBITOL
python SHAP_training.py -df ${scratch}/geno_rf_YPXYLOSE.csv -sep , -test ${data}/Test.txt -feat ${scratch}/feat_exp_rf_YPXYLOSE_top_256 -model ${scratch}/YPXYLOSE_exp_rf_256_models.pkl -top 256 -save YPXYLOSE


scontrol show job $SLURM_JOB_ID