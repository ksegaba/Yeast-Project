#!/usr/bin/env python3
""" Manuscript supplementary tables code """

import os
import re
import tqdm
import pickle
import swifter
import glob
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
from sklearn.metrics import r2_score
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import spearmanr

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

# Isolate growth condition labels; will be used throughout the script
mapping = {"YPACETATE":"YP Acetate 2%", "YPD14":"YPD 14ºC", "YPD40":"YPD 40ºC",
		   "YPD42":"YPD 42ºC", "YPD6AU":"YPD 6-Azauracile 600 µg/ml",
		   "YPDANISO10":"YPD Anisomycin 10 µg/ml", "YPDANISO20":"YPD Anisomycin 20 µg/ml",
		   "YPDANISO50":"YPD Anisomycin 50 µg/ml", "YPDBENOMYL200":"YPD Benomyl 200 µg/ml",
		   "YPDBENOMYL500":"YPD Benomyl 500 µg/ml", "YPDCAFEIN40":"YPD Caffeine 40 mM",
		   "YPDCAFEIN50":"YPD Caffeine 50 mM", "YPDCHX05":"YPD Cycloheximide 0.5 µg/ml",
		   "YPDCHX1":"YPD Cycloheximide 1 µg/ml", "YPDCUSO410MM":"YPD CuSO4 10 mM",
		   "YPDDMSO":"YPD DMSO 6%", "YPDETOH":"YPD Ethanol 15%",
		   "YPDFLUCONAZOLE":"YPD Fluconazole 20 µg/ml", "YPDFORMAMIDE4":"YPD Formamide 4%",
		   "YPDFORMAMIDE5":"YPD Formamide 5%", "YPDHU":"YPD Hydroxyurea 30 mg/ml",
		   "YPDKCL2M":"YPD KCL 2 M", "YPDLICL250MM":"YPD LiCl 250 mM",
		   "YPDMV":"YPD Methylviologen 20 mM", "YPDNACL15M":"YPD NaCl 1.5 M",
		   "YPDNACL1M":"YPD NaCl 1 M", "YPDNYSTATIN":"YPD Nystatin 10 µg/ml",
		   "YPDSDS":"YPD SDS 0.2%", "YPDSODIUMMETAARSENITE":"YPD Sodium metaarsenite 2.5 mM",
		   "YPETHANOL":"YP Ethanol 2%", "YPGALACTOSE":"YP Galactose 2%",
		   "YPRIBOSE":"YP Ribose 2%", "YPGLYCEROL":"YP Glycerol 2%",
		   "YPXYLOSE":"YP Xylose 2%", "YPSORBITOL":"YP Sorbitol 2%"}

################################################################################ #********* done 3/19/2024 for RF PAV and CNV; need RF PC and SNP; need to re-do FS for all alg
### TABLE S1
################################################################################
## RANDOM FOREST (RF) PREDICTION PERFORMANCES (BASELINE USING ALL FEATURES)
## PC Models
path = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/baseline/RESULTS_reg.txt"
rf_pc = pd.read_csv(path, sep="\t")
rf_pc.insert(0, "new_cond", rf_pc.replace({"Y": mapping})["Y"]) # add full condition names
rf_pc = rf_pc.sort_values(by="r2_test", ascending=False)
cond_order = rf_pc["new_cond"]
rf_pc.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t", index=False) #********* done 3/19/2024

## SNP Models
path = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results"
rf_base = pd.read_csv(f"{path}/baseline/RESULTS_SNP_reg.txt", sep="\t")
rf_base.insert(0, "new_cond", rf_base.replace({"Y": mapping})["Y"]) # add full condition names
rf_base.set_index("new_cond", inplace=True)
rf_base = rf_base.loc[cond_order] # sort by PC order
rf_base.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_baseline.txt", sep="\t") #********* do when sds and nystatin finish running

## PAV Models
path = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
rf_orf = pd.read_csv(f"{path}/baseline/RESULTS_reg.txt", sep="\t")
rf_pav = rf_orf.loc[rf_orf.ID.str.contains("_pav_baseline")]
rf_pav.insert(0, "new_cond", rf_pav.replace({"Y":mapping})["Y"])
rf_pav.set_index("new_cond", inplace=True)
rf_pav = rf_pav.loc[cond_order]
rf_pav.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_baseline.txt", sep="\t") #********* done 3/19/2024

## CNV Models
rf_cnv = rf_orf[rf_orf.ID.str.contains("_cnv_baseline")]
rf_cnv.insert(0, "new_cond", rf_cnv.replace({"Y":mapping})["Y"])
rf_cnv.set_index("new_cond", inplace=True)
rf_cnv = rf_cnv.loc[cond_order]
rf_cnv.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_baseline.txt", sep="\t") #********* done 3/19/2024

## Merge them into one table
combined = pd.concat([rf_pc, rf_base.reset_index(), rf_pav.reset_index(), rf_cnv.reset_index()])
combined.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_ALL_baseline.txt", sep="\t", index=None)

## BASELINE PREDICTION PERFORMANCES OF ALL OTHER ALGORITHMS
## PCs and SNPs
d = "/mnt/gs21/scratch/seguraab/yeast_project/"
xgb_snp = pd.read_csv(d + "SNP_yeast_XGBoost_results/baseline/RESULTS_xgboost.txt", sep="\t") # XGBoost results
rrblup_snp = pd.read_csv(d + "SNP_yeast_rrBLUP_results/baseline/RESULTS_rrblup.txt", sep="\t") # rrBLUP results
bl_snp = pd.read_csv(d + "SNP_yeast_BL_results/baseline/RESULTS_BL.txt", sep="\t") # Bayesian LASSO results
bayesc_snp = pd.read_csv(d + "SNP_yeast_BayesC_results/baseline/RESULTS_BayesC.txt", sep="\t") # BayesC results

# Keep only relevant columns and insert missing information
xgb_snp.rename(columns={"R2_val            ":"R2_val", "R2_test            ":"R2_test"}, inplace=True)
xgb_snp = xgb_snp[["Data", "Trait", "R2_val", "R2_val_sd", "PCC_val", "PCC_val_sd", 
				   "R2_test", "R2_test_sd", "PCC_test", "PCC_test_sd"]]
xgb_snp.insert(4, "R2_val_se", xgb_snp.R2_val_sd/np.sqrt(10))
xgb_snp.insert(7, "PCC_val_se", xgb_snp.PCC_val_sd/np.sqrt(10))
xgb_snp.insert(10, "R2_test_se", xgb_snp.R2_test_sd/np.sqrt(10))
xgb_snp.insert(13, "PCC_test_se", xgb_snp.PCC_test_sd/np.sqrt(10))
xgb_snp.insert(2, "NumInstances", 625)
xgb_snp.insert(3, "CV-Fold", 5)
xgb_snp.insert(4, "NumRepetitions", 10)
xgb_snp.insert(3, "NumFeatures", np.repeat(5,35).tolist()+np.repeat(118382,35).tolist())
xgb_snp.insert(0, "Alg", "XGBoost")

rrblup_snp.NumInstances = 625
rrblup_snp.insert(0, "Alg", "rrBLUP")
rrblup_snp.insert(1, "Data", np.repeat("SNPs",35).tolist()+np.repeat("PCs_tassel",35).tolist())
rrblup_snp.drop(columns=["Date", "ID"], inplace=True)

bl_snp.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_snp.insert(0, "Data", np.repeat("PCs_tassel",35).tolist()+np.repeat("SNPs",35).tolist())
bl_snp.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					   "FeatureNum":"NumFeatures"}, inplace=True)
bl_snp.columns = bl_snp.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

bayesc_snp.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_snp.insert(0, "Data", np.repeat("PCs_tassel",35).tolist()+np.repeat("SNPs",35).tolist())
bayesc_snp.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					   "FeatureNum":"NumFeatures"}, inplace=True)
bayesc_snp.columns = bayesc_snp.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

## PAVs and CNVs
d = "/mnt/gs21/scratch/seguraab/yeast_project/"
xgb_orf = pd.read_csv(d + "ORF_yeast_XGBoost_results/baseline/RESULTS_xgboost.txt", sep="\t") # XGBoost results
rrblup_orf = pd.read_csv(d + "ORF_yeast_rrBLUP_results/baseline/RESULTS_rrblup.txt", sep="\t") # rrBLUP results
bl_orf = pd.read_csv(d + "ORF_yeast_BL_results/baseline/RESULTS_BL.txt", sep="\t") # Bayesian LASSO results
bayesc_orf = pd.read_csv(d + "ORF_yeast_BayesC_results/baseline/RESULTS_BayesC.txt", sep="\t") # BayesC results

# Keep only relevant columns and insert missing information
xgb_orf.rename(columns={"R2_val            ":"R2_val", "R2_test            ":"R2_test"}, inplace=True)
xgb_orf = xgb_orf[["Data", "Trait", "R2_val", "R2_val_sd", "PCC_val", "PCC_val_sd", 
				   "R2_test", "R2_test_sd", "PCC_test", "PCC_test_sd"]]
xgb_orf.insert(4, "R2_val_se", xgb_orf.R2_val_sd/np.sqrt(10))
xgb_orf.insert(7, "PCC_val_se", xgb_orf.PCC_val_sd/np.sqrt(10))
xgb_orf.insert(10, "R2_test_se", xgb_orf.R2_test_sd/np.sqrt(10))
xgb_orf.insert(13, "PCC_test_se", xgb_orf.PCC_test_sd/np.sqrt(10))
xgb_orf.insert(2, "NumInstances", 625)
xgb_orf.insert(3, "CV-Fold", 5)
xgb_orf.insert(4, "NumRepetitions", 10)
xgb_orf.insert(3, "NumFeatures", 7708)
xgb_orf.insert(0, "Alg", "XGBoost")

rrblup_orf.NumInstances = 625
rrblup_orf.insert(0, "Alg", "rrBLUP")
rrblup_orf["Data"] = rrblup_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1)
rrblup_orf.drop(columns=["Date", "ID"], inplace=True)

bl_orf.insert(0, "Data", bl_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))
bl_orf.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_orf.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					   "FeatureNum":"NumFeatures"}, inplace=True)
bl_orf.columns = bl_orf.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

bayesc_orf = bayesc_orf.loc[~bayesc_orf.ID.str.contains("__"),:]
bayesc_orf = bayesc_orf.loc[~bayesc_orf.ID.duplicated(),:]
bayesc_orf.insert(0, "Data", bayesc_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))
bayesc_orf.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_orf.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					   "FeatureNum":"NumFeatures"}, inplace=True)
bayesc_orf.columns = bayesc_orf.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

# Merge them into one table
combined.insert(0, "Data", combined.apply(lambda x: "PCs_tassel" if "PC" in x.ID else \
								 ("PAV" if "orf" in x.ID else \
								 ("CNV" if "cno" in x.ID else "SNP")), axis=1))
combined.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se", "Tag", "new_cond", \
						 "EVS_val", "EVS_val_sd", "EVS_val_se", "EVS_test", "EVS_test_sd", \
						 "EVS_test_se"], inplace=True)
combined.rename(columns={"Y":"Trait", "FeatureNum":"NumFeatures", "CVfold":"CV-Fold", \
						 "CV_rep":"NumRepetitions"}, inplace=True)
combined.columns = combined.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

all_alg = pd.concat([combined, xgb_snp, rrblup_snp, bl_snp, bayesc_snp, \
					 xgb_orf, rrblup_orf, bl_orf, bayesc_orf])
all_alg.insert(1, "new_cond", all_alg.replace({"Trait":mapping})["Trait"])
all_alg.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t", index=False)

## PREDICTION PERFORMANCES OF ALL ALGORITHMS USING RF FEATURE SELECTION FEATURES
# First, determine optimal feature subsets based on RF model performances
d = "/mnt/gs21/scratch/seguraab/yeast_project"

## SNPs
rf_snp = pd.read_csv(f"{d}/SNP_yeast_RF_results/fs/RESULTS_reg.txt", sep="\t") # Random Forest results
rf_snp.insert(0, "new_cond", rf_snp.replace({"Y": mapping})["Y"]) # add full condition names
rf_fs = pd.DataFrame(columns=rf_snp.columns) # Optimal FS features
for env in mapping.keys():
	tmp = rf_snp.loc[rf_snp["ID"].str.match(f"{env}_(exp_rf|rf|top)_[0-9]+"),:]
	rf_fs = pd.concat([rf_fs, tmp.loc[tmp.r2_val==max(tmp.r2_val),:]])
	print(tmp.loc[tmp.r2_val==max(tmp.r2_val),"ID"].values)
rf_fs.set_index("new_cond", inplace=True)
rf_fs = rf_fs.loc[cond_order,] # set to PC env order
rf_fs.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t", index=True)

## PAVs and CNVs
rf_orf = pd.read_csv(f"{d}/ORF_yeast_RF_results/RESULTS_reg.txt", sep="\t") # Random Forest results
rf_orf = rf_orf.loc[rf_orf.ID.str.match("[A-Z0-9]+_(orf|cno)_[0-9]+"),:]
rf_orf.insert(0, "new_cond", rf_orf.replace({"Y": mapping})["Y"]) # add full condition names
rf_pav = pd.DataFrame(columns=rf_orf.columns) # Optimal FS features
rf_cnv = pd.DataFrame(columns=rf_orf.columns)
for env in mapping.keys():
	tmp_pav = rf_orf.loc[rf_orf["ID"].str.match(f"{env}_orf_[0-9]+"),:]
	rf_pav = pd.concat([rf_pav, tmp_pav.loc[tmp_pav.r2_val==max(tmp_pav.r2_val),:]])
	print(tmp_pav.loc[tmp_pav.r2_val==max(tmp_pav.r2_val),"ID"].values)
	tmp_cnv = rf_orf.loc[rf_orf["ID"].str.match(f"{env}_cno_[0-9]+"),:]
	rf_cnv = pd.concat([rf_cnv, tmp_cnv.loc[tmp_cnv.r2_val==max(tmp_cnv.r2_val),:]])
	print(tmp_cnv.loc[tmp_cnv.r2_val==max(tmp_cnv.r2_val),"ID"].values)
rf_pav.set_index("new_cond", inplace=True)
rf_cnv.set_index("new_cond", inplace=True)
rf_pav = rf_pav.loc[cond_order,] # set to PC env order
rf_cnv = rf_cnv.loc[cond_order,]
rf_pav.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t", index=True)
rf_cnv.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t", index=True)

# Load in results from models from the other algorithms run with the RF feature subsets
## SNPs
xgb_snp = pd.read_csv(d + "/SNP_yeast_XGBoost_results/fs/RESULTS_xgboost.txt", sep="\t") # XGBoost results
rrblup_snp = pd.read_csv(d + "/SNP_yeast_rrBLUP_results/fs_based_on_rf/RESULTS_rrblup.txt", sep="\t")
bl_snp = pd.read_csv(d + "/SNP_yeast_BL_results/fs/RESULTS_BL.txt", sep="\t") # Bayesian LASSO results
bayesc_snp = pd.read_csv(d + "/SNP_yeast_BayesC_results/fs/RESULTS_BayesC.txt", sep="\t") # BayesC results

# Keep only relevant columns and insert missing information
xgb_snp.rename(columns={"R2_val            ":"R2_val",
				"R2_test            ":"R2_test", "NumFeatures            ":"NumFeatures"},
				inplace=True)
xgb_snp = xgb_snp[["Data", "Trait", "NumInstances", "NumFeatures", "CV-Fold",
				   "NumRepetitions", "R2_val", "R2_val_sd", "PCC_val", "PCC_val_sd",
				   "R2_test", "R2_test_sd", "PCC_test", "PCC_test_sd"]]
xgb_snp.insert(8, "R2_val_se", xgb_snp.R2_val_sd/np.sqrt(10))
xgb_snp.insert(11, "PCC_val_se", xgb_snp.PCC_val_sd/np.sqrt(10))
xgb_snp.insert(14, "R2_test_se", xgb_snp.R2_test_sd/np.sqrt(10))
xgb_snp.insert(17, "PCC_test_se", xgb_snp.PCC_test_sd/np.sqrt(10))
xgb_snp.insert(0, "Alg", "XGBoost")

rrblup_snp.NumInstances = 625
rrblup_snp.insert(0, "Alg", "rrBLUP")
rrblup_snp.insert(1, "Data", "SNP")
rrblup_snp.drop(columns=["Date", "ID"], inplace=True)

bl_snp.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_snp.rename(columns={"FeatureNum":"NumFeatures"}, inplace=True)
bl_snp.insert(4, "CV-Fold", 5)
bl_snp.insert(5, "NumRepetitions", 10)
bl_snp.insert(0, "Data", "SNP")
bl_snp.columns = bl_snp.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

bayesc_snp.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_snp.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					   "FeatureNum":"NumFeatures"}, inplace=True)
bayesc_snp.insert(0, "Data", "SNP")
bayesc_snp.columns = bayesc_snp.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

rf_fs.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se", "EVS_val", "EVS_val_sd",
				"EVS_val_se", "EVS_test", "EVS_test_sd", "EVS_test_se", "Tag"], inplace=True)
rf_fs.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					"FeatureNum":"NumFeatures", "Y":"Trait"}, inplace=True)
rf_fs.insert(0, "Data", "SNP")
rf_fs.columns = rf_fs.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

## PAVs and CNVs
xgb_orf = pd.read_csv(d + "/ORF_yeast_XGBoost_results/fs/RESULTS_xgboost.txt", sep="\t")
rrblup_orf = pd.read_csv(d + "/ORF_yeast_rrBLUP_results/fs/RESULTS_rrblup.txt", sep="\t")
bl_orf = pd.read_csv(d + "/ORF_yeast_BL_results/fs/RESULTS_BL.txt", sep="\t")
bayesc_orf = pd.read_csv(d + "/ORF_yeast_BayesC_results/fs/RESULTS_BayesC.txt", sep="\t")

# Keep only relevant columns and insert missing information
xgb_orf.rename(columns={"R2_val            ":"R2_val",
				"R2_test            ":"R2_test", "NumFeatures            ":"NumFeatures"},
				inplace=True)
xgb_orf = xgb_orf[["Data", "Trait", "NumInstances", "NumFeatures", "CV-Fold",
				   "NumRepetitions", "R2_val", "R2_val_sd", "PCC_val", "PCC_val_sd",
				   "R2_test", "R2_test_sd", "PCC_test", "PCC_test_sd"]]
xgb_orf.insert(8, "R2_val_se", xgb_orf.R2_val_sd/np.sqrt(10))
xgb_orf.insert(11, "PCC_val_se", xgb_orf.PCC_val_sd/np.sqrt(10))
xgb_orf.insert(14, "R2_test_se", xgb_orf.R2_test_sd/np.sqrt(10))
xgb_orf.insert(17, "PCC_test_se", xgb_orf.PCC_test_sd/np.sqrt(10))
xgb_orf.insert(0, "Alg", "XGBoost")

rrblup_orf.NumInstances = 625
rrblup_orf.insert(0, "Alg", "rrBLUP")
rrblup_orf.insert(1, "Data", rrblup_orf.apply(lambda x: "PAV" if \
	"pav" in x.ID else "CNV", axis=1))
rrblup_orf.drop(columns=["Date", "ID"], inplace=True)

bl_orf.rename(columns={"FeatureNum":"NumFeatures"}, inplace=True)
bl_orf.insert(7, "CV-Fold", 5)
bl_orf.insert(8, "NumRepetitions", 10)
bl_orf.insert(0, "Data", bl_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))
bl_orf.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_orf.columns = bl_orf.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

bayesc_orf.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					   "FeatureNum":"NumFeatures"}, inplace=True)
bayesc_orf.insert(0, "Data", bayesc_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))
bayesc_orf.drop(columns=["Date", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_orf.columns = bayesc_orf.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

rf_pav.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se", "EVS_val", "EVS_val_sd",
				"EVS_val_se", "EVS_test", "EVS_test_sd", "EVS_test_se", "Tag"], inplace=True)
rf_pav.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					"FeatureNum":"NumFeatures", "Y":"Trait"}, inplace=True)
rf_pav.insert(0, "Data", "PAV")
rf_pav.columns = rf_pav.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

rf_cnv.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se", "EVS_val", "EVS_val_sd",
				"EVS_val_se", "EVS_test", "EVS_test_sd", "EVS_test_se", "Tag"], inplace=True)
rf_cnv.rename(columns={"CVfold":"CV-Fold", "CV_rep":"NumRepetitions",
					"FeatureNum":"NumFeatures", "Y":"Trait"}, inplace=True)
rf_cnv.insert(0, "Data", "CNV")
rf_cnv.columns = rf_cnv.apply(lambda x: x.name if "R2" in x.name else x.name.replace("r2", "R2"), axis=0)

## Merge all single-environment models together (baseline and FS)
all_baseline = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t")
all_single_env = pd.concat([rf_fs, xgb_snp, rrblup_snp, bl_snp, bayesc_snp, rf_pav, rf_cnv,
			xgb_orf, rrblup_orf, bl_orf, bayesc_orf, all_baseline], axis=0)
all_single_env.drop(columns=["CVfold", "CV_rep"], inplace=True)
all_single_env["new_cond"] = all_single_env.replace({"Trait": mapping})["Trait"]
all_single_env.set_index("new_cond", inplace=True)
all_single_env.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_SINGLE_ENV.txt", sep="\t")

all_fs = pd.concat([rf_fs, xgb_snp, rrblup_snp, bl_snp, bayesc_snp, rf_pav, rf_cnv,
			xgb_orf, rrblup_orf, bl_orf, bayesc_orf], axis=0)
all_fs.drop(columns=["CVfold", "CV_rep"], inplace=True)
all_fs["new_cond"] = all_fs.replace({"Trait": mapping})["Trait"]
all_fs.set_index("new_cond", inplace=True)
all_fs.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_FS.txt", sep="\t")

################################################################################
### TABLE S2 ---- perhaps S2
################################################################################
## linear regression of model performance on fitness-related factors
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)  # fitness data
var = pheno.var(axis=0) ; var.name = "var"
med = pheno.median(axis=0) ; med.name = "median"

h2 = pd.read_csv("Scripts/Data_Vis/Section_2/Heritability_h2_H2_sommer.csv") # trait heritabilities
h2.set_index("Conditions", inplace=True)
X = pd.concat([var, med, h2.h2], ignore_index=False, axis=1) # fitness factors
X_s = (X-X.mean())/X.std() # center and scale

for data_type in ["SNPs", "PAVs", "CNVs", "PCs"]:
	if data_type == "PCs": # model performance results
		Y = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")
	else:
		Y = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_%s_FS.txt" % data_type, sep="\t")
	Y.set_index("Y", inplace=True)
	## Regress model performance on all factors at once
	df = pd.concat([X_s, Y.r2_test], axis=1)
	mod = smf.ols(formula=f"r2_test ~ {' * '.join(df.columns[:3])}", data=df)
	res = mod.fit()
	yhats = res.predict()
	r2_scores = r2_score(df.r2_test, yhats, multioutput=None) # double check sm.ols R2
	with open(f"Scripts/Data_Vis/Section_3/factors_ols_{data_type}_results.txt", "w") as out:
		out.write(res.summary().as_text())
		out.write(f"\nR-sq: {r2_scores}")
	# vars(res) # attributes
	pickle.dump(mod, open(f"Scripts/Data_Vis/Section_3/factors_ols_{data_type}_model.pkl", 'wb')) # save the model
	yhats=pd.Series(yhats)
	yhats.index = df.index
	yhats.name = 'y_pred'
	pd.concat([Y.r2_test, yhats], ignore_index=False, axis=1).\
		to_csv(f"Scripts/Data_Vis/Section_3/factors_ols_{data_type}_preds.csv")
	del df, Y, mod

## GxE models linear regression
# model performance regressed on degree of correlation among the correlated environments

################################################################################ #********* done 3/19/2024 for SNP, PAV, and CNV baseline only
### SUPPLEMENTARY DATA FILE 10
###############################################################################
## Combine RF FS & baseline model gini importance values for SNPs, PAVs, and CNVs individually
# read feature to gene map files
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)
map_snps.merge(map_orfs, how="inner", on="gene").gene.nunique() # 5356 shared genes
map_snps["gene_with_intergenic"] = map_snps.apply(lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1)

# paths to feature importance score files
dir = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs"
snp_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t")  # SNP FS results
snp_fs_files = [os.path.join(dir, f"{x}_imp") for x in snp_rf_res['ID']]
dir = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/baseline"
snp_baseline_files = [os.path.join(dir, f"{x}_rf_baseline_imp") for x in mapping.keys()]
dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs"
pav_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t")  # ORF pres/abs FS results
pav_fs_files = [os.path.join(dir, f"{x}_imp") for x in pav_rf_res['ID']]
cnv_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t")  # CNV FS results
cnv_fs_files = [os.path.join(dir, f"{x}_imp") for x in cnv_rf_res['ID']]
dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline"
pav_baseline_files = [os.path.join(dir, f"{x}_pav_baseline_imp") for x in mapping.keys()]
cnv_baseline_files = [os.path.join(dir, f"{x}_cnv_baseline_imp") for x in mapping.keys()]

## combine gini importance for all envs per data type
def combine_imp_indiv(imp_files, map=map_snps, dtype="snp", save="", mapping=mapping):
	for i,env in enumerate(mapping.keys()):
		print(env)
		# Read gini importance file
		file = [f for f in imp_files if env in f]
		print(len(file)) # should be 1
		imp = dt.fread(file[0]).to_pandas()
		imp.set_index(imp.iloc[:,0], inplace=True) # feature names as index
		imp = imp.loc[:,"mean_imp"] # use mean gini importances
		imp.rename(env, inplace=True)
		imp = pd.DataFrame(imp)
		if dtype != "snp":
			imp.index = imp.apply(lambda x: re.sub("^X", "", x.name), axis=1) # rename index
			imp.index = imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
		if i == 0:
			merged = imp.copy(deep=True)
		else:
			merged = pd.concat([merged, imp], axis=1, ignore_index=False) # add to dictionary
		del imp
	print(merged.shape)
	# map to genes
	if dtype == "snp":
		merged = map_snps[["snp","gene","gene_with_intergenic"]].merge(merged, how="right", left_on="snp", right_index=True)
	else:
		merged = map_orfs[["orf","gene"]].merge(merged, how="right", left_on="orf", right_index=True)
	merged.to_csv(save, sep="\t", index=False)
	return merged

combine_imp_indiv(snp_fs_files, save="Scripts/Data_Vis/Section_4/RF_FS_imp_snp.tsv")
combine_imp_indiv(pav_fs_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_FS_imp_pav.tsv")
combine_imp_indiv(cnv_fs_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_FS_imp_cnv.tsv")
combine_imp_indiv(snp_baseline_files, save="Scripts/Data_Vis/Section_4/RF_baseline_imp_snp.tsv")
combine_imp_indiv(pav_baseline_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_baseline_imp_pav.tsv")
combine_imp_indiv(cnv_baseline_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_baseline_imp_cnv.tsv")

## Combine feature importance information from different environments into one table including all data types
# def combine_imp(imp_files, merged={}, merged_50={}, dtype="snp", mapping=mapping):
	# for env in mapping.keys():
	#     if env not in merged.keys():
	#         merged[env]={}
	#     if env not in merged_50.keys():
	#         merged_50[env]={}
	#     file = [f for f in imp_files if env in f] # Read feature importance score files
	#     imp = pd.read_csv(file[0], sep="\t", skiprows=1, names=["X", "mean_imp"])
	#     if dtype == "snp": # Map markers to genes and drop intergenic SNPs
	#         imp = map_snps[["gene", "snp"]].merge(imp, left_on="snp", right_on="X")
	#         imp = imp[imp["gene"] != "intergenic"]
	#     elif dtype != "snp": # Map ORFs to genes and drop those with no mapping
	#         imp.X = imp.X.apply(lambda x: re.sub("^X", "", x))
	#         imp.X = imp.X.apply(lambda x: re.sub("\.", "-", x))
	#         imp = map_orfs.merge(imp, left_on="orf", right_on="X")
	#     imp = imp.groupby("gene")["mean_imp"].max().reset_index()  # Get max importance score per gene
	#     imp.rename(columns={"mean_imp":"max_mean_imp"}, inplace=True)
	#     imp = imp.sort_values(by="max_mean_imp", ascending=False) # Rank genes
	#     imp["rank"] = imp.max_mean_imp.rank(method="average", ascending=False)
	#     if dtype == "snp": # Calculate rank percentiles
	#         imp["rank_per"] = imp["rank"] / num_genes_snps ### note 1/22/2024: I need to change to 1 - n/N because my bins below identify >=.95 as top genes. Right now, I was identifying low rank genes as top
	#     elif dtype != "snp":
	#         imp["rank_per"] = imp["rank"] / num_genes_orfs
	#     imp["rank_per_bin"] = pd.cut(imp["rank_per"],
	#         bins=[0, 0.01, 0.05, 0.1, 1],
	#         labels=["(0.99, 1]", "(0.95, 0.99]", "(0.9, 0.95]", "(0, 0.90]"]) # Bin percentiles into categories
	#     if dtype not in merged[env].keys():
	#         merged[env][dtype] = imp[["gene", "max_mean_imp", "rank_per", "rank_per_bin"]] # Merge feature importance measures
	#     if dtype not in merged_50[env].keys():
	#         imp_50 = imp.iloc[:50,:] # Top 50 genes
	#         new_idx = imp_50.gene.str.split(" // ", expand=True).stack().reset_index(level=1, drop=True).reset_index()
	#         new_idx.columns = ["index", "gene"]
	#         imp_50 = new_idx.merge(imp_50, left_on="gene", right_on="gene")
	#         imp_50.drop(columns=["index", "rank"], inplace=True)
	#         merged_50[env][dtype] = imp_50
#     return merged, merged_50


# def dict2df(nested_dict, idx_col="gene"):
#     """
#     Return a pandas dataframe with a column multi-index and and a row index
#     Arguments:
#             nested_dict: nested dictionary with pandas dataframes in the
#                             innermost values. The keys will become the column
#                             multi-index.
#             idx_col: column to set as row index.
#     """
#     new_dict = {}
#     for outer_key, inner_dict in nested_dict.items():# Iterate through the nested dictionary
#         for inner_key, df in inner_dict.items(): # Append the current dataframe and corresponding index information
#             new_dict[(outer_key, inner_key)] = df.set_index(idx_col)
#     result_df = pd.concat(new_dict, axis=1, ignore_index=False) # Concatenate the dataframes along the new multi-index
#     return result_df

# # combine feature importance data
# snp_baseline_merged, snp_base_merged_50 = combine_imp(snp_baseline_files, merged={}, merged_50={})
# orf_baseline_merged, orf_base_merged_50 = combine_imp(orf_baseline_files, merged=snp_baseline_merged, merged_50=snp_base_merged_50, dtype="orf")
# cnv_baseline_merged, cnv_base_merged_50 = combine_imp(cnv_baseline_files, merged=orf_baseline_merged, merged_50=orf_base_merged_50, dtype="cnv")
# snp_merged, snp_merged_50 = combine_imp(snp_files, merged={}, merged_50={})
# orf_merged, orf_merged_50 = combine_imp(orf_files, merged=snp_merged, merged_50=snp_merged_50, dtype="orf")
# cnv_merged, cnv_merged_50 = combine_imp(cnv_files, merged=orf_merged, merged_50=orf_merged_50, dtype="cnv")
# fs_combined = dict2df(cnv_merged) ; fs_combined_50 = dict2df(cnv_merged_50)
# baseline_combined = dict2df(cnv_baseline_merged) ; base_combined_50 = dict2df(cnv_base_merged_50)
# fs_combined.to_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types.tsv", sep="\t")
# fs_combined_50.to_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types_top50.tsv", sep="\t")
# baseline_combined.to_csv("Scripts/Data_Vis/RF_baseline_imp_all_data_types.tsv", sep="\t")
# base_combined_50.to_csv("Scripts/Data_Vis/RF_baseline_imp_all_data_types_top50.tsv", sep="\t")

# # correlations between baseline feature importances for different datatypes (how to do this for FS models?)
# imp = pd.read_csv("Scripts/Data_Vis/RF_baseline_imp_all_data_types.tsv", sep="\t", index_col=0, header=[0,1,2])
# corr_df = {"env":[], "comparison":[], "gene_num":[], "r":[], "p-val":[]}
# for env in imp.columns.get_level_values(0).unique():
#     snp_vs_pav = pd.concat([imp[env, "snp", "max_mean_imp"], imp[env, "orf", "max_mean_imp"]], axis=1).dropna()
#     snp_vs_cnv = pd.concat([imp[env, "snp", "max_mean_imp"], imp[env, "cnv", "max_mean_imp"]], axis=1).dropna()
#     pav_vs_cnv = pd.concat([imp[env, "orf", "max_mean_imp"], imp[env, "cnv", "max_mean_imp"]], axis=1).dropna()
#     r, p = pearsonr(snp_vs_pav.iloc[:,0], snp_vs_pav.iloc[:,1])
#     corr_df["env"].append(env)
#     corr_df["comparison"].append("snp_vs_pav")
#     corr_df["gene_num"].append(snp_vs_pav.shape[0])
#     corr_df["r"].append(r)
#     corr_df["p-val"].append(p)
#     r, p = pearsonr(snp_vs_cnv.iloc[:,0], snp_vs_cnv.iloc[:,1])
#     corr_df["env"].append(env)
#     corr_df["comparison"].append("snp_vs_cnv")
#     corr_df["gene_num"].append(snp_vs_cnv.shape[0])
#     corr_df["r"].append(r)
#     corr_df["p-val"].append(p)
#     r, p = pearsonr(pav_vs_cnv.iloc[:,0], pav_vs_cnv.iloc[:,1])
#     corr_df["env"].append(env)
#     corr_df["comparison"].append("pav_vs_cnv")
#     corr_df["gene_num"].append(pav_vs_cnv.shape[0])
#     corr_df["r"].append(r)
#     corr_df["p-val"].append(p)
#     del snp_vs_pav, snp_vs_cnv, pav_vs_cnv

# pd.DataFrame.from_dict(corr_df).sort_values('r', ascending=False).\
#     to_csv("Scripts/Data_Vis/RF_baseline_imp_data_type_correlations.tsv", sep="\t", index=False)


################################################################################
### ALSO IN SUPPLEMENTARY DATA FILE 10
################################################################################
## Combine RF FS & baseline model SHAP values for SNPs, PAVs, and CNVs individually
# paths to feature SHAP value files
dir = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP"
snp_shap_files = [f"{dir}/SNP/fs/{file}" for file in os.listdir(dir + "/SNP/fs") if file.startswith("SHAP_values_sorted_average_Y")]
snp_shap_baseline_files = [f"{dir}/SNP/baseline/{file}" for file in os.listdir(dir + "/SNP/baseline") if file.startswith("SHAP_values_sorted_average_Y")]
pav_shap_files = [f"{dir}/PAV/fs/{file}" for file in os.listdir(dir + "/PAV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
pav_shap_baseline_files = [f"{dir}/PAV/baseline/{file}" for file in os.listdir(dir + "/PAV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_files = [f"{dir}/CNV/fs/{file}" for file in os.listdir(dir + "/CNV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_baseline_files = [f"{dir}/CNV/baseline/{file}" for file in os.listdir(dir + "/CNV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]

# read feature to gene map files
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)
map_snps["gene_with_intergenic"] = map_snps.apply(lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1) # label intergenic snps

# combine shap values for all envs per data type
def combine_shap_indiv(shap_files, map=map_snps, merged={}, dtype="snp", save="", mapping=mapping):
	for i,env in enumerate(mapping.keys()):
		print(env)
		# Read SHAP file
		file = [f for f in shap_files if env in f]
		print(len(file)) # should be 1
		shap = dt.fread(file[0]).to_pandas()
		shap.set_index(shap.iloc[:,0], inplace=True)
		shap = shap.iloc[:,1:]
		shap.rename(columns={"C1":env}, inplace=True)
		if dtype != "snp":
			shap.index = shap.apply(lambda x: re.sub("^X", "", x.name), axis=1) # rename index
			shap.index = shap.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
		if i == 0:
			merged = shap.copy(deep=True)
		else:
			merged = pd.concat([merged, shap], axis=1, ignore_index=False) # add to dictionary
		del shap
	print(merged.shape)
	# map to genes
	if dtype == "snp":
		merged = map_snps[["snp","gene","gene_with_intergenic"]].merge(merged, how="right", left_on="snp", right_index=True)
	else:
		merged = map_orfs[["orf","gene"]].merge(merged, how="right", left_on="orf", right_index=True)
	merged.to_csv(save, sep="\t", index=False)
	return merged

combine_shap_indiv(snp_shap_files, save="Scripts/Data_Vis/Section_4/RF_FS_shap_snp.tsv")
combine_shap_indiv(pav_shap_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_FS_shap_pav.tsv")
combine_shap_indiv(cnv_shap_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_FS_shap_cnv.tsv")
combine_shap_indiv(snp_shap_baseline_files, save="Scripts/Data_Vis/Section_4/RF_baseline_shap_snp.tsv")
combine_shap_indiv(pav_shap_baseline_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_baseline_shap_pav.tsv")
combine_shap_indiv(cnv_shap_baseline_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_baseline_shap_cnv.tsv")

################################################################################
### TABLE S4 (CORRELATIONS) & SUPPLEMENTARY DATA FILE 11 (RANK PERCENTILE MATRICES)
################################################################################
## SPEARMAN'S RHO CORRELATIONS BETWEEN IMPORTANCE MEASURES
res = [["Model Type", "Data Type", "Env", "rho", "pval"]]
for mod_type in ["baseline", "FS"]:
	for data_type in ["snp", "pav", "cnv"]:
		gini = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_imp_{data_type}.tsv").to_pandas()
		shap = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_shap_{data_type}.tsv").to_pandas()
		gini.index = gini.iloc[:,0]
		gini = gini.iloc[:,1:]
		shap.index = shap.iloc[:,0]
		shap = shap.iloc[:,1:]
		## Convert data to rank percentiles
		for env in mapping.keys():
			gini[env] = gini[env].rank(axis=0, method="average", numeric_only=True, pct=True) # ranks as percentiles (1= most important)
			shap[env] = shap[env].abs().rank(axis=0, method="average", numeric_only=True, pct=True) # based on absolute shap value
		gini.to_csv(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_imp_{data_type}_rank_per.tsv")
		shap.to_csv(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_shap_{data_type}_rank_per.tsv")
		## Calculate spearman's rho
		for env in mapping.keys():
			df = pd.concat([gini.loc[:,env], shap.loc[:,env]], ignore_index=False, axis=1).dropna()
			rho = df.corr(method=lambda x, y: spearmanr(x,y).statistic)
			pval = df.corr(method=lambda x, y: spearmanr(x,y).pvalue)
			res.append([mod_type, data_type, env, rho.iloc[0,1], pval.iloc[0,1]])
res = pd.DataFrame(res)
res.columns = res.iloc[0,:]
res = res.iloc[1:,:]
res.sort_values(by='rho', ascending=False, inplace=True)
res.to_csv("Scripts/Data_Vis/Section_4/Table_S_gini_vs_shap_rank_per_corr.tsv", sep="\t", index=False)

## Hypothesis test to prove ranking distributions are not random



## Commented out 02/28/2024
# ## Scatterplots gini vs imp per data type
# snp_imp = dt.fread("Scripts/Data_Vis/RF_baseline_imp_snp_rank_per.tsv").to_pandas()
# pav_imp = dt.fread("Scripts/Data_Vis/RF_baseline_imp_pav_rank_per.tsv").to_pandas()
# cnv_imp = dt.fread("Scripts/Data_Vis/RF_baseline_imp_cnv_rank_per.tsv").to_pandas()
# snp_shap = dt.fread("Scripts/Data_Vis/RF_baseline_shap_snp_rank_per.tsv").to_pandas()
# pav_shap = dt.fread("Scripts/Data_Vis/RF_baseline_shap_pav_rank_per.tsv").to_pandas()
# cnv_shap = dt.fread("Scripts/Data_Vis/RF_baseline_shap_cnv_rank_per.tsv").to_pandas()
# snp_imp.set_index("snp", inplace=True)
# snp_shap.set_index("snp", inplace=True)
# pav_imp.set_index("orf", inplace=True)
# pav_shap.set_index("orf", inplace=True)
# cnv_imp.set_index("orf", inplace=True)
# cnv_shap.set_index("orf", inplace=True)
# def make_density_plot(df1, df2):
#     fig, axes = plt.subplots(7, 5, figsize=(25, 28))
#     for i, env in enumerate(tqdm.tqdm(mapping.keys())):
#         values = np.vstack([df1[env], df2[env]])
#         kernel = stats.gaussian_kde(values)(values) # point kernel density
#         row_idx, col_idx = divmod(i, 5) # subplot coordinates
#         ax = axes[row_idx, col_idx] # set subplot
#         sns.scatterplot(x=df1[env], y=df2[env], c=kernel, edgecolor = 'none',
#                         size=-5, alpha=.3, cmap="viridis", legend=False, ax=ax) # density scatter plot
#         ax.set_title(f'{mapping[env]}')
#         norm = plt.Normalize(kernel.min(), kernel.max()) # normalize kde values
#         cbar = plt.cm.ScalarMappable(norm=norm, cmap="viridis") # colorbar
#         ax.figure.colorbar(cbar, ax=ax) # set colorbar
#         del values, kernel, row_idx, col_idx, ax, norm , cbar
#     return fig, axes
# pav_fig, pav_axes = make_density_plot(pav_imp, pav_shap)
# plt.tight_layout()
# plt.savefig("Scripts/Data_Vis/Density_gini_vs_shap_pav_baseline.pdf")
# plt.close()
# cnv_fig, cnv_axes = make_density_plot(cnv_imp, cnv_shap)
# plt.tight_layout()
# plt.savefig("Scripts/Data_Vis/Density_gini_vs_shap_cnv_baseline.pdf")
# plt.close()
# pav_vs_cnv_fig, pav_vs_cnv_axes = make_density_plot(pav_imp, cnv_imp)
# plt.tight_layout()
# plt.savefig("Scripts/Data_Vis/Density_pav_vs_cnv_imp_baseline.pdf")
# plt.close()
# pav_vs_cnv_fig, pav_vs_cnv_axes = make_density_plot(pav_shap, cnv_shap)
# plt.tight_layout()
# plt.savefig("Scripts/Data_Vis/Density_pav_vs_cnv_shap_baseline.pdf")
# plt.close()
# snp_fig, snp_axes = make_density_plot(snp_imp, snp_shap) # Note 1/18/24 takes a long time, would need to submit a job
# plt.tight_layout()
# plt.savefig("Scripts/Data_Vis/Density_gini_vs_imp_snp_baseline.pdf")
# plt.close()

################################################################################
### TABLE S4 TOO???
################################################################################
## SPEARMAN'S RHO CORRELATIONS BEWTEEN DATA TYPES
res = [["Model Type", "Importance Type", "Env", "Comparison", "rho", "pval"]]
for mod_type in ["baseline", "FS"]:
	for imp_type in ["imp", "shap"]:
		df_snp = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_snp_rank_per.tsv").to_pandas()
		df_pav = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_pav_rank_per.tsv").to_pandas()
		df_cnv = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_cnv_rank_per.tsv").to_pandas()
		## Calculate spearman's rho
		if imp_type == "imp":
			imp_score_type = "mean gini importance"
		else:
			imp_score_type = "mean absolute SHAP value"
		# prep data for SNP to PAV/CNV comparisons
		tmp_snp = df_snp.loc[~df_snp.gene_with_intergenic.str.contains('intergenic'),:] # drop intergenic snps
		tmp_pav = df_pav.loc[df_pav.gene != "",:] # drop ORFs with no gene mapping
		tmp_cnv = df_cnv.loc[df_cnv.gene != "",:]
		tmp_snp.drop(columns="snp", inplace=True) # drop snp column
		tmp_pav.drop(columns=["orf"], inplace=True) # drop orf and organism columns
		tmp_cnv.drop(columns=["orf"], inplace=True)
		tmp_snp = tmp_snp.groupby("gene_with_intergenic").max() # get the max rank percentile per gene
		tmp_pav = tmp_pav.groupby("gene").max()
		tmp_cnv = tmp_cnv.groupby("gene").max()
		for env in mapping.keys():
			## SNP vs PAV
			try:
				df = pd.concat([tmp_snp.loc[:,env], tmp_pav.loc[:,env]], ignore_index=False, axis=1).dropna()
				rho = df.corr(method=lambda x, y: spearmanr(x,y).statistic)
				pval = df.corr(method=lambda x, y: spearmanr(x,y).pvalue)
				res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs", \
							rho.iloc[0,1], pval.iloc[0,1]])
				del df
			except:
				res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs", \
					None, None])
			## SNP vs CNV
			try:
				df = pd.concat([tmp_snp.loc[:,env], tmp_cnv.loc[:,env]], ignore_index=False, axis=1).dropna()
				rho = df.corr(method=lambda x, y: spearmanr(x,y).statistic)
				pval = df.corr(method=lambda x, y: spearmanr(x,y).pvalue)
				res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs", \
							rho.iloc[0,1], pval.iloc[0,1]])
				del df
			except:
				res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs", \
					None, None])
			## PAV vs CNV (keep all ORFs regardless of gene mapping)
			try:
				df = pd.concat([df_pav.set_index("orf").loc[:,env], \
								df_cnv.set_index("orf").loc[:,env]], ignore_index=False, \
								axis=1).dropna()
				rho = df.corr(method=lambda x, y: spearmanr(x,y).statistic)
				pval = df.corr(method=lambda x, y: spearmanr(x,y).pvalue)
				res.append([mod_type, imp_score_type, env, "PAVs vs CNVs", \
							rho.iloc[0,1], pval.iloc[0,1]])
				del df
			except:
				res.append([mod_type, imp_score_type, env, "PAVs vs CNVs", \
					None, None])
		del df_snp, df_pav, df_cnv, tmp_snp, tmp_pav, tmp_cnv
res = pd.DataFrame(res)
res.columns = res.iloc[0,:]
res = res.iloc[1:,:]
res.to_csv("Scripts/Data_Vis/Section_4/Table_S_rank_per_corr_btwn_data_types.tsv", sep="\t", index=False)

### Don't think I need below anymore 12/18/2023; at least I know it needs to be modified
## Combine RF FS and baseline model SHAP values for SNPs, PAVs, and CNVs by gene
# combine SHAP values into one table
# def combine_shap(shap_files, merged={}, merged_50={}, dtype="snp", drop_intergenic:bool=True, baseline:bool=False, mapping=mapping):
#     for env in mapping.keys():
#         print(env)
#         if env not in merged.keys():
#             merged[env]={}
#         if env not in merged_50.keys():
#             merged_50[env]={}
#         file = [f for f in shap_files if env in f] # Read SHAP files
#         if not file:
#             print(f"No file found for {env}")
#             return
#         else:
#             if baseline: # Read in data
#                 shap = dt.fread(file[0]).to_pandas()
#                 shap.set_index(shap.iloc[:,0], inplace=True)
#                 shap = shap.iloc[:,1:]
#                 shap = shap.T
#             else:
#                 shap = pd.read_csv(file[0], sep="\t", header=0, index_col=0).T
#             if dtype == "snp": # Map markers to genes and drop intergenic SNPs
#                 shap = map_snps[["gene", "snp"]].merge(shap, left_on="snp", right_index=True)
#                 if drop_intergenic:
#                     shap = shap[shap["gene"] != "intergenic"]
#                 else: # Rename intergenic snps to intergenic_<snp>
#                     shap["gene"] = shap.apply(lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1)
#                 shap.drop(columns="snp", inplace=True)
#             if dtype != "snp": # Map ORFs to genes and drop those with no mapping
#                 shap.index = shap.apply(lambda x: re.sub("^X", "", x.name), axis=1)
#                 shap.index = shap.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
#                 shap = map_orfs.merge(shap, left_on="orf", right_index=True)
#                 shap.drop(columns=["organism", "orf"], inplace=True)
#             shap.set_index("gene", inplace=True)
#             shap = shap.mean(axis=1)
#             shap = shap.groupby("gene").max().reset_index() # Get max importance score per gene
#             shap.rename(columns={0:"max_mean_shap"}, inplace=True)
#             shap = shap.sort_values(by="max_mean_shap", ascending=False)
#             shap["rank"] = shap.max_mean_shap.rank(method="average", ascending=False) # Rank genes
#             if dtype == "snp": # Calculate rank percentiles
#                 if drop_intergenic:
#                     shap["rank_per"] = shap["rank"] / num_genes_snps ### note 1/22/2024: I need to change to 1 - n/N because my bins below identify >=.95 as top genes. Right now, I was identifying low rank genes as top
#                 else:
#                     shap["rank_per"] = shap["rank"] / num_genes_snps_with_intergenic
#             elif dtype != "snp":
#                 shap["rank_per"] = shap["rank"] / num_genes_orfs
#             shap["rank_per_bin"] = pd.cut(shap["rank_per"],
#                 bins=[0, 0.01, 0.05, 0.1, 1],
#                 labels=["(0.99, 1]", "(0.95, 0.99]", "(0.9, 0.95]", "(0, 0.90]"]) # Bin percentiles into categories
#             if dtype not in merged[env].keys():
#                 merged[env][dtype] = shap.drop(columns="rank")
#             if dtype not in merged_50[env].keys():
#                 shap_50 = shap.iloc[:50,:] # Top 50 genes
#                 new_idx = shap_50.gene.str.split(" // ", expand=True).stack().reset_index(level=1, drop=True).reset_index()
#                 new_idx.columns = ["index", "gene"]
#                 shap_50 = new_idx.merge(shap_50, left_on="gene", right_on="gene")
#                 shap_50.drop(columns=["index", "rank"], inplace=True)
#                 merged_50[env][dtype] = shap_50
#     return merged, merged_50

# # Merge single-env RF FS shap results into a master table
# snp_shap_merged, snp_shap_merged_50 = combine_shap(snp_shap_files, merged={}, merged_50={}, dtype="snp")
# pav_shap_merged, pav_shap_merged_50 = combine_shap(pav_shap_files, merged=snp_shap_merged, merged_50=snp_shap_merged_50, dtype="pav")
# cnv_shap_merged, cnv_shap_merged_50 = combine_shap(cnv_shap_files, merged=pav_shap_merged, merged_50=pav_shap_merged_50, dtype="cnv")
# dict2df(cnv_shap_merged).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types.tsv", sep="\t") # function on line 252
# dict2df(cnv_shap_merged_50).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types_top50.tsv", sep="\t")
# snp_shap_merged_with_intergenic, snp_shap_merged_50_with_intergenic = combine_shap(snp_shap_files, merged={}, merged_50={}, dtype="snp", drop_intergenic=False)
# pav_shap_merged_with_intergenic, pav_shap_merged_50_with_intergenic = combine_shap(pav_shap_files, merged=snp_shap_merged_with_intergenic, merged_50=snp_shap_merged_50_with_intergenic, dtype="pav", drop_intergenic=False)
# cnv_shap_merged_with_intergenic, cnv_shap_merged_50_with_intergenic = combine_shap(cnv_shap_files, merged=pav_shap_merged_with_intergenic, merged_50=pav_shap_merged_50_with_intergenic, dtype="cnv", drop_intergenic=False)
# dict2df(cnv_shap_merged_with_intergenic).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types_with_intergenic.tsv", sep="\t")
# dict2df(cnv_shap_merged_50_with_intergenic).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types_with_intergenic_top50.tsv", sep="\t")

# # Merge single-env baseline RF shap results into a master table
# snp_shap_merged_base, snp_shap_merged_base_50 = combine_shap(snp_shap_baseline_files, merged={}, merged_50={}, dtype="snp", baseline=True)
# pav_shap_merged_base, pav_shap_merged_base_50 = combine_shap(pav_shap_baseline_files, merged=snp_shap_merged_base, merged_50=snp_shap_merged_base_50, dtype="pav", baseline=True)
# cnv_shap_merged_base, cnv_shap_merged_base_50 = combine_shap(cnv_shap_baseline_files, merged=pav_shap_merged_base, merged_50=pav_shap_merged_base_50, dtype="cnv", baseline=True)
# dict2df(cnv_shap_merged_base).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types.tsv", sep="\t")
# dict2df(cnv_shap_merged_base_50).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types_top50.tsv", sep="\t")
# snp_shap_merged_base_with_intergenic, snp_shap_merged_base_50_with_intergenic = combine_shap(snp_shap_baseline_files, merged={}, merged_50={}, dtype="snp", drop_intergenic=False, baseline=True)
# pav_shap_merged_base_with_intergenic, pav_shap_merged_base_50_with_intergenic = combine_shap(pav_shap_baseline_files, merged=snp_shap_merged_base_with_intergenic, merged_50=snp_shap_merged_base_50_with_intergenic, dtype="pav", drop_intergenic=False, baseline=True)
# cnv_shap_merged_base_with_intergenic, cnv_shap_merged_base_50_with_intergenic = combine_shap(cnv_shap_baseline_files, merged=pav_shap_merged_base_with_intergenic, merged_50=pav_shap_merged_base_50_with_intergenic, dtype="cnv", drop_intergenic=False, baseline=True)
# dict2df(cnv_shap_merged_base_with_intergenic).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types_with_intergenic.tsv", sep="\t")
# dict2df(cnv_shap_merged_base_50_with_intergenic).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types_with_intergenic_top50.tsv", sep="\t")

################################################################################
### TABLE S7 perhaps???
################################################################################
# Call script to create GO enrichment and pathway enrichment excel files
os.system("python Scripts/Data_Vis/combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs ^ORA_SHAP_values_sorted_average_ SNP go Scripts/Data_Vis")
os.system("python Scripts/Data_Vis/combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs ^ORA_SHAP_values_sorted_average_[A-Z0-9]+_orf ORF go Scripts/Data_Vis")
os.system("python Scripts/Data_Vis/combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs ^ORA_SHAP_values_sorted_average_[A-Z0-9]+_cno CNV go Scripts/Data_Vis")
os.system("python Scripts/Data_Vis/combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs PWY_ORA_SHAP_values_sorted_average_ SNP pwy Scripts/Data_Vis")
os.system("python Scripts/Data_Vis/combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs PWY_ORA_SHAP_values_sorted_average_[A-Z0-9]+_orf ORF pwy Scripts/Data_Vis")
os.system("python Scripts/Data_Vis/combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs PWY_ORA_SHAP_values_sorted_average_[A-Z0-9]+_cno CNV pwy Scripts/Data_Vis")
np_merged, snp_merged_50 = combine_imp(snp_files)
# Some GO terms had no BP, MF, or CC annotations so I added them manually by searching on AmiGO 2

################################################################################
### TABLE S3 --- perhaps S7 ---- Perhaps only do for lit genes and move to a later sheet, so maybe S10 or something like that
################################################################################
## SPEARMAN's RHO CORRELATIONS BETWEEN ENVIRONMENTS
for i,mod_type in enumerate(["baseline", "FS"]):
	for j,data_type in enumerate(["snp", "pav", "cnv"]):
		for k,imp_type in enumerate(["imp", "shap"]):
			print(f"Calculating spearman's rho for Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_{data_type}_rank_per.tsv")
			df = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_{data_type}_rank_per.tsv").to_pandas()
			df.set_index(df.columns[0], inplace=True)
			# calculate spearman's rho
			if data_type != "snp":
				rho = df.iloc[:,1:].corr(method=lambda x, y: spearmanr(x, y).statistic)
				pval = df.iloc[:,1:].corr(method=lambda x, y: spearmanr(x, y).pvalue)
				print(df.iloc[:,1:].shape)
				print(rho.isna().sum())
				print(pval.isna().sum())
			else:
				rho = df.iloc[:,2:].corr(method=lambda x, y: spearmanr(x, y).statistic)
				pval = df.iloc[:,2:].corr(method=lambda x, y: spearmanr(x, y).pvalue)
				print(df.iloc[:,2:].shape)
			# df[["YPACETATE", "YPSORBITOL"]].dropna().corr(method=lambda x, y: spearmanr(x,y).statistic) # check to make sure correlations are calculated correctly for FS files; looks good!
			# prepare table to save
			print(rho.shape, pval.shape)
			# print(rho)
			# print(pval)
			# print(rho.where(np.triu(np.ones(rho.shape), k=1).astype(bool)).stack())
			out = rho.where(np.triu(rho, k=1).astype(bool)).stack()
			out = pd.concat([out, pval.where(np.triu(pval, k=1).astype(bool)).stack()], ignore_index=False, axis=1)
			out.columns = ["rho", "pval"]
			out.pval.fillna(0.0, inplace=True) # 0.0 pvalues become NaNs with astype(bool)
			out.index.set_names(["Env1", "Env2"], inplace=True)
			print(out.shape)
			# out.sort_values(by="rho", ascending=False).\
			#     to_csv(f"Scripts/Data_Vis/Section_4/Table_S_RF_{mod_type}_{imp_type}_{data_type}_rank_per_corr_btwn_envs.tsv", sep="\t")
			if (i==0 and j==0 and k==0):
				res = out.copy(deep=True)
				res.columns = [f"{imp_type}_{data_type}_{mod_type}_rho", \
							   f"{imp_type}_{data_type}_{mod_type}_pvalue"]
			else:
				out.columns = [f"{imp_type}_{data_type}_{mod_type}_rho", \
							   f"{imp_type}_{data_type}_{mod_type}_pvalue"]
				res = pd.concat([res, out], ignore_index=False, axis=1)
			del df, rho, pval, out
res.to_csv("Scripts/Data_Vis/Section_4/Table_S_rank_per_corr_btwn_envs.tsv", sep="\t")

# Note: For FS models, the feature subsets may not overlap and hence no rho or pval can be calculated.
# This is why for SNP FS models there are some missing env pairs for gini and shap.

## Hypothesis test to prove ranking distributions are not random


################################################################################
### TABLE S8
################################################################################
## ENRICHMENT OF BENCHMARK STRESS-RESPONSE GENES IN RF FS MODELS
# Target phenotype and mutant information annotations
phenotypes = ["resistance to chemicals: decreased ", "viability: decreased ",
			  "resistance to chemicals: decreased", "metal resistance: decreased",
			  "metal resistance: decreased ", "oxidative stress resistance: decreased",
			  "respiratory growth: decreased", "stress resistance: decreased"] # sensitive mutants
mutants = ["null Allele", "reduction of function", "reduction of function Allele"]

## Genes lists from SGD
benomyl_sgd = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations.txt", skiprows=8, sep="\t")
benomyl_sgd[["Mutant Type", "Mutant Description"]] = benomyl_sgd["Mutant Information"].str.split(":", expand=True)
# benomyl_sgd.groupby("Mutant Type").count() # see mutant types
ben_genes = benomyl_sgd.loc[(benomyl_sgd.Phenotype.isin(phenotypes)) & \
							(benomyl_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
ben_genes.to_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

caffeine_sgd = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations.txt", skip_to_line=9, sep="\t").to_pandas()
caffeine_sgd[["Mutant Type", "Mutant Description"]] = caffeine_sgd["Mutant Information"].str.split(":", expand=True)
caf_genes = caffeine_sgd.loc[(caffeine_sgd.Phenotype.isin(phenotypes)) & \
							(caffeine_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
caf_genes.to_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

cuso4_sgd = dt.fread("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations.txt", skip_to_line=9, sep="\t").to_pandas()
cuso4_sgd[["Mutant Type", "Mutant Description"]] = cuso4_sgd["Mutant Information"].str.split(":", expand=True)
cu_genes  = cuso4_sgd.loc[(cuso4_sgd.Phenotype.isin(phenotypes)) & \
							(cuso4_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
cu_genes.to_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

sodmetars_sgd = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations.txt", skiprows=8, sep="\t")
sodmetars_sgd[["Mutant Type", "Mutant Description"]] = sodmetars_sgd["Mutant Information"].str.split(":", expand=True)
sma_genes = sodmetars_sgd.loc[(sodmetars_sgd.Phenotype.isin(phenotypes)) & \
							(sodmetars_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
sma_genes.to_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

## Enrichment of literature genes within the feature selection models' features
# Remove literature genes that are not found within our datasets
map_snps = pd.read_csv("Scripts/Data_Vis/biallelic_snps_diploid_and_S288C_genes.txt",
					   header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)
# map_orfs.to_csv("~/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", sep="\t", index=False) ### I have to regenerate the snp map because I accidentally overwrote it :'(
map_orfs.to_csv("~/Shiu_Lab/Project/Data/Peter_2018/final_map_orf_to_gene_16_removed.txt", sep="\t", index=False)

ben_genes = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf_genes = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt").to_pandas()
cu_genes = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma_genes = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")
curated = pd.read_csv("Data/SGD_Experiment_Genes/manually_curated_genes.txt", sep="\t")

ben_snp = ben_genes.loc[ben_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
caf_snp = caf_genes.loc[caf_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
cu_snp = cu_genes.loc[cu_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
sma_snp = sma_genes.loc[sma_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
curated_snp = curated.loc[curated["gene"].isin(map_snps.gene), "gene"].unique()
pd.Series(ben_snp).to_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(caf_snp).to_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(cu_snp).to_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(sma_snp).to_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(curated_snp).to_csv("Data/SGD_Experiment_Genes/manually_curated_genes_snps.txt", sep="\t", index=False)

ben_orf = ben_genes.loc[ben_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
caf_orf = caf_genes.loc[caf_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
cu_orf = cu_genes.loc[cu_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
sma_orf = sma_genes.loc[sma_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
curated_orf = curated.loc[curated.gene.isin(map_orfs.gene), "gene"].unique()
pd.Series(ben_orf).to_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(caf_orf).to_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(cu_orf).to_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(sma_orf).to_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(curated_orf).to_csv("Data/SGD_Experiment_Genes/manually_curated_genes_orfs.txt", sep="\t", index=False)

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"] # most predictive envs
lit_gene_envs = ["benomyl", "caffeine", "copper(II) sulfate", "sodium meta-arsenite", "manually curated"]

# Create contingency table
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
def enrichment_direction(k, n, C, G):
    """determine direction of enrichment
    if >= 1: + (overrepresented)
    if < 1: - (underrepresented)
    k: number of genes in top 5% and are known literature genes
    n: total number of genes in top 5%
    C: total number of genes (in top 5% + background) and are known literature genes
    G: total number of genes (in top 5% + background)"""
    return((k/C)/(n/G))

# Check for enrichment at different ranking percentiles
for percentile in [0.99, 0.95, 0.90, 0.85, 0.80, 0.75]:
	# Contingency table header
	enrich_res = [] # collect results from enrichment test
	for data_type in ["snp", "pav", "cnv"]:
		for imp_type in ["imp", "shap"]:
			print(percentile, data_type, imp_type)
			# Read and prepare dataset
			df_fs = dt.fread(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}_rank_per.tsv").to_pandas()
			df_base = dt.fread(f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_{data_type}_rank_per.tsv").to_pandas()
			df = pd.concat([])
			df = df.iloc[1:,]
			#
			# Count number of literature genes identified by each model
			for env in target_envs:
				env_df = df[["gene", env]].dropna()
				env_df_rank = env_df.loc[env_df[env]!=0.0,:] # Exclude genes with 0 importance
				env_df_rank["rank_per"] = env_df_rank.rank(method="average", numeric_only=True, pct=True)
				env_df_rank.sort_values(by="rank_per", ascending=False, inplace=True)
				env_df_genes = env_df_rank.gene.unique()
				#
				# Variable 1: How many genes are known stress-response genes?
				if data_type=="snp":
					num_ben = np.intersect1d(env_df_genes, ben_snp)
					num_caf = np.intersect1d(env_df_genes, caf_snp)
					num_cu = np.intersect1d(env_df_genes, cu_snp)
					num_sma = np.intersect1d(env_df_genes, sma_snp)
					num_man = np.intersect1d(env_df_genes, curated_snp)
				else:
					num_ben = np.intersect1d(env_df_genes, ben_orf)
					num_caf = np.intersect1d(env_df_genes, caf_orf)
					num_cu = np.intersect1d(env_df_genes, cu_orf)
					num_sma = np.intersect1d(env_df_genes, sma_orf)
					num_man = np.intersect1d(env_df_genes, curated_orf)
				#
				# Variable 2: How many genes are NOT known stress-response genes?
				if data_type=="snp":
					num_not_ben = np.setdiff1d(env_df_genes, ben_snp, assume_unique=True)
					num_not_caf = np.setdiff1d(env_df_genes, caf_snp, assume_unique=True)
					num_not_cu = np.setdiff1d(env_df_genes, cu_snp, assume_unique=True)
					num_not_sma = np.setdiff1d(env_df_genes, sma_snp, assume_unique=True)
					num_not_man = np.setdiff1d(env_df_genes, curated_snp, assume_unique=True)
				else:
					num_not_ben = np.setdiff1d(env_df_genes, ben_orf, assume_unique=True)
					num_not_caf = np.setdiff1d(env_df_genes, caf_orf, assume_unique=True)
					num_not_cu = np.setdiff1d(env_df_genes, cu_orf, assume_unique=True)
					num_not_sma = np.setdiff1d(env_df_genes, sma_orf, assume_unique=True)
					num_not_man = np.setdiff1d(env_df_genes, curated_orf, assume_unique=True)
				#
				# Variable 3: How many genes are ranked in the top 5%?
				env_df_rank_top = env_df_rank.loc[env_df_rank.rank_per > percentile,:]
				num_top = len(env_df_rank_top)
				#
				# Variable 4: How many genes are NOT ranked in the top 5%?
				num_background = len(env_df_rank) - len(env_df_rank_top)
				#
				## Create contingency table
				# for benomyl lit genes
				for i,var in enumerate([num_ben, num_caf, num_cu, num_sma, num_man]):
					a = len(env_df_rank_top.loc[env_df_rank_top.gene.isin(var),:]) # gene is a known literature gene AND is in top 5%
					b = len(env_df_rank_top.loc[~env_df_rank_top.gene.isin(var),:]) # gene is not known AND is in top 5%
					c = len(env_df_rank.loc[(env_df_rank.rank_per <= percentile) & (env_df_rank.gene.isin(var)),:]) # gene is known AND is not in top 5%
					d = len(env_df_rank.loc[(env_df_rank.rank_per <= percentile) & (~env_df_rank.gene.isin(var)),:]) # gene is not known AND is not in top 5%
					odds, pval = fisher_exact([[a,c],[b,d]]) # Run fisher's exact test
					#
					# Direction of enrichment
					if enrichment_direction(k=a, n=a+b, C=a+c, G=a+b+c+d) >= 1:
						direction = "+"
					else:
						direction = "-"
					#
					# Save results
					enrich_res.append(["RF_single-env_baseline", data_type, imp_type, env, lit_gene_envs[i], a,\
					b, c, d, odds, pval, np.log2(odds), np.log10(pval), direction])
	#
	enrich_res = pd.DataFrame(enrich_res)
	enrich_res.columns = ["Model Type", "DataType", "ImpType", "Env", "LitGeneList",
					f"Known & above {int(percentile*100)}%",
					f"Not known & above {int(percentile*100)}%",
					f"Known & not above {int(percentile*100)}%",
					f"Not known & not above {int(percentile*100)}%", "Odds Ratio",
					"p-value", "log2(odds ratio)", "log10(p-value)", "direction"]
	enrich_res_sub = enrich_res.loc[enrich_res["p-value"] != 1.0,:]
	q_values = false_discovery_control(enrich_res_sub["p-value"], method="bh") # Adjust the p-values
	enrich_res_sub["q-values"] = q_values
	enrich_res_sub["log10(q-values)"] = enrich_res_sub.apply(lambda x: np.log10(x["q-values"]) if x.direction=="-" else -np.log10(x["q-values"]), axis=1)
	enrich_res_sub.to_csv(f"Scripts/Data_Vis/Section_4/Enrichment_of_literature_genes_in_above_{int(percentile*100)}%_RF_baseline.csv", index=False)
	#
	# Visualize the significant enrichment values
	significant = enrich_res_sub.loc[enrich_res_sub["q-values"] <= 0.05, ["DataType", "ImpType", "Env", "LitGeneList", "log2(odds ratio)"]]
	if significant.shape[0] > 0:
		for i_type in significant.ImpType.unique():
			significant_sub = significant.loc[significant.ImpType==i_type,:]
			print(significant)
			print(significant_sub)
			for dat_type in significant_sub.DataType.unique():
				significant_sub2 = significant_sub.loc[significant_sub.DataType==dat_type, \
					["Env", "LitGeneList", "log2(odds ratio)"]].\
					pivot(index="Env", columns="LitGeneList", values="log2(odds ratio)") # pivot wider
				significant_sub2.replace([np.inf, -np.inf], np.nan, inplace=True) # set -Inf to NaN
				plt.figure(figsize=(7, 11))
				sns.heatmap(data=significant_sub2, cmap="RdBu_r", center=0, mask=significant_sub2.isna())
				plt.tight_layout()
				plt.savefig(f"Scripts/Data_Vis/Section_4/Enrichment_of_literature_genes_in_above_{int(percentile*100)}_RF_{dat_type}_baseline_{i_type}.pdf")
				plt.close()
			del significant_sub2, significant_sub

################################################################################
### TABLE S8 (this is outdated, see 0_TablS8_kstest.py for full code) -- perhaps S9
################################################################################
## Compare ranking distribution of literature genes to random chance
target_envs = {"YPDCAFEIN40":"Caffeine", "YPDCAFEIN50":"Caffeine",
			   "YPDBENOMYL500":"Benomyl", "YPDCUSO410MM":"CuSO4",
			   "YPDSODIUMMETAARSENITE":"Sodium (meta)arsenite"} # most predictive envs

# Remove literature genes that are not found within our datasets
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
					   header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)

num_genes_snps = map_snps['gene'].nunique() - 1  # don't count intergenic
num_genes_orfs = map_orfs['gene'].nunique()

# lit_genes = pd.read_csv("Data/Peter_2018/literature_genes.txt", sep="\t") # target genes from the literature
# lit_genes_snps = lit_genes.loc[lit_genes.Gene.isin(map_snps.gene),:].reset_index().iloc[:,1:]
# lit_genes_orfs = lit_genes.loc[lit_genes.Gene.isin(map_orfs.gene),:].reset_index().iloc[:,1:]

# List to collect results
permut_res = [["Model Type", "DataType", "ImpType", "Env", "NumGenes", \
			   "NumLitGenesFound", "AvgRankPerSampDist", "AvgRankPerLitGenes", \
			   "Z-score", "p-value"]] # collect results from permutation test

# Compare max rank percentile of literature genes identified by the model to random chance
n = 1000000 # number of repetitions
seeds = 3*2*5*n # random seed counter
for data_type in ["snp"]:#, "pav", "cnv"]:
	for imp_type in ["imp"]:#, "shap"]:
		# Read and prepare dataset
		df = dt.fread(f"Scripts/Data_Vis/RF_baseline_{imp_type}_{data_type}.tsv").to_pandas()
		if data_type == "snp":
			df.set_index(["snp", "gene"], inplace=True)
			df = df.iloc[:,1:] # remove "gene_with_intergenic" column; ranking will not includer intergenic snps
			df = df.loc[df.index.get_level_values("gene") != "intergenic",:] # drop intergenic snps
		else:
			df.set_index(["orf", "gene"], inplace=True)
		if imp_type == "shap":
			df = df.abs()
		## Compare ranking distribution of literature genes to random chance
		for env in target_envs.keys():
			print(data_type, imp_type, env)
			# Rank the genes based on gini importance or shap value
			env_imp = df.loc[:,env].dropna().sort_values(ascending=False) # all the genic snps for this env, excluding intergenic snps
			max_env_ranks = env_imp.groupby("gene").max() # take the max importance (either gini or shap) per gene
			max_env_ranks.sort_values(ascending=False, inplace=True) # sort from highest importance to lowest
			max_env_ranks = pd.concat([max_env_ranks, max_env_ranks.rank(method="average", pct=True)], axis=1) # assign ordinal rank
			max_env_ranks.columns = ["importance", "rank_per"]
			max_lit_ranks = max_env_ranks.loc[max_env_ranks.index.get_level_values("gene")\
								.isin(lit_genes.loc[lit_genes.Env==target_envs[env]].Gene),:]
			## Begin random sampling
			if max_lit_ranks.shape[0] > 1: # needed for calculating std of t-test and z-test
				rand_ranks = [] # distribution of mean rank of the literature genes found in the top 5% ranking
				for i in tqdm.tqdm(range(n)):
					rand_ranks.append(max_env_ranks.sample(n=max_lit_ranks.shape[0],
									random_state=seeds, ignore_index=True).rank_per.mean()) # average max rank percentile of random sample
					seeds -= 1
				# Figure to check for normality
				fig, ax = plt.subplots(nrows=1, ncols=2)
				ax[0].hist(max_lit_ranks["rank_per"], alpha=0.5, label="original")
				ax[1].hist(rand_ranks, alpha=0.5, label="random sample")
				ax[1].vlines(x=np.mean(rand_ranks), ymin=0, ymax=2000, colors="black", linestyles="dashed", label="avg random")
				ax[1].vlines(x=max_lit_ranks.rank_per.mean(), ymin=0, ymax=2000, colors="red", linestyles="dashed", label="avg lit genes")
				ax[0].legend(loc="upper right") ; ax[1].legend(loc="upper right")
				plt.savefig(f"Scripts/Data_Vis/Samp_dist_rank_per_means_{data_type}_{imp_type}_{env}.png")
				plt.close()
				## Right-tailed Z-test
				z_score = (max_lit_ranks.rank_per.mean() - np.mean(rand_ranks)) / \
					np.std(rand_ranks)
				# z_score = (np.mean(rand_ranks) - max_lit_ranks.rank_per.mean()) / \
				#     max_lit_ranks.rank_per.std()
				z_critical = stats.norm.ppf(1-0.05) # for alpha = 0.05, the z_critical = 1.6448536269514722
				# note: if z_score >= z_critical, reject the null hypothesis
				z_p_val = 1 - stats.norm.cdf(np.abs(z_score))
				# Save results
				permut_res.append(["RF FS", data_type, imp_type, env, len(max_env_ranks), \
				len(max_lit_ranks), np.mean(rand_ranks), max_lit_ranks.rank_per.mean(), \
				z_score, z_p_val])
			else:
				permut_res.append(["RF FS", data_type, imp_type, env, len(max_env_ranks), \
				len(max_lit_ranks), len(max_lit_ranks), np.mean(rand_ranks), None, None])

pd.DataFrame(permut_res).to_csv("Scripts/Data_Vis/Table_S8_mean_lit_gene_rank_per_vs_random.tsv",
	sep="\t", header=False, index=False)


################################################################################ #********* done 3/19/2024 for snp, pav, and cnv
### TABLE S10
################################################################################
"""Generate the maximum gini importance feature subsets from SNP, PAV, CNV
RF baseline models that will be used to build a new set of XGBoost models, which
only contains one feature per gene, for SHAP interaction score analysis"""

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/RF"

for data_type in ["snp", "pav", "cnv"]:
	# Read and prepare dataset
	df = dt.fread(f"Scripts/Data_Vis/Section_4/RF_baseline_imp_{data_type}.tsv").to_pandas()
	#
	if data_type == "snp":
		df.set_index(["snp", "gene"], inplace=True)
		df = df.iloc[:,1:] # remove "gene_with_intergenic" column; ranking will not include intergenic snps
		df = df.loc[df.index.get_level_values("gene") != "intergenic",:] # drop intergenic snps
	else:
		df.gene = df.apply(lambda x: x.orf if x.gene=="" else x.gene, axis=1) # add orfs with no gene matches to gene column
		df.set_index(["orf", "gene"], inplace=True)
	#
	for env in target_envs:
		# Take the max importance (either gini or shap) per gene
		env_imp = df.loc[:,env].dropna().sort_values(ascending=False)
		env_imp = env_imp.loc[env_imp != 0.0,:] # remove features with 0 gini importance
		max_features = pd.concat([env_imp.groupby("gene").idxmax().apply(lambda x: x[0]),
								  env_imp.groupby("gene").max()],
								ignore_index=False, axis=1) 
		#
		if data_type=="snp":
			max_features.columns = ["snp", "max_imp"]
		else:
			max_features.columns = ["orf", "max_imp"]
			max_features.orf = max_features.apply(lambda x: "X" + x.orf, axis=1)
			max_features.orf = max_features.apply(lambda x: re.sub("-", ".", x.orf), axis=1)
		#
		max_features.to_csv(
			f"{d}/Feature_map_max_gini_from_RF_baseline_imp_{data_type}_{env}.csv")
		max_features.iloc[:,0].to_csv(
			f"{d}/Features_max_gini_from_RF_baseline_imp_{data_type}_{env}.txt",
			index=False, header=False)

################################################################################
### TABLE S11
################################################################################
"""After running features selection on the XGBoost models used for SHAP interaction
analysis, concatenate the literature genes to the top features and run a new set
of XGBoost models. Also run the models with the label randomized. The purpose of
these is to obtain SHAP gene-gene interactions between top genes and literature
genes and also the feature SHAP values. 
First kind: top FS features plus either benomyl, caffeine, cuso4, or sma
Second kind: top FS features plus all lit genes combined"""

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs"
res = pd.read_csv(f"{d}/RESULTS_xgboost.txt", sep="\t") # XGB FS performance results

# Genes lists from SGD
ben_genes = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf_genes = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt").to_pandas()
cu_genes = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma_genes = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")
lit_envs = ["benomyl", "caffeine", "cuso4", "sma"]

# SNP and ORF gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
					   header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_16_removed.txt", sep="\t")

# Get literature genes found in SNP features
ben_snp = map_snps.swifter.apply(lambda x: (x.snp,x.gene) if \
	(x.gene in ben_genes["Gene Systematic Name"].values) else (None,None), axis=1)
ben_snp = pd.DataFrame(ben_snp.swifter.apply(pd.Series)).dropna()
caf_snp = map_snps.swifter.apply(lambda x: (x.snp,x.gene) if \
	(x.gene in caf_genes["Gene Systematic Name"].values) else (None,None), axis=1)
caf_snp = pd.DataFrame(caf_snp.swifter.apply(pd.Series)).dropna()
cu_snp = map_snps.swifter.apply(lambda x: (x.snp,x.gene) if \
	(x.gene in cu_genes["Gene Systematic Name"].values) else (None,None), axis=1)
cu_snp = pd.DataFrame(cu_snp.swifter.apply(pd.Series)).dropna()
sma_snp = map_snps.swifter.apply(lambda x: (x.snp,x.gene) if \
	(x.gene in sma_genes["Gene Systematic Name"].values) else (None,None), axis=1)
sma_snp = pd.DataFrame(sma_snp.swifter.apply(pd.Series)).dropna()

# Get literature genes found in ORF features
ben_orf = map_orfs.swifter.apply(lambda x: (x.orf,x.gene) if \
	(x.gene in ben_genes["Gene Systematic Name"].values) else (None,None), axis=1)
ben_orf = pd.DataFrame(ben_orf.swifter.apply(pd.Series)).dropna()
caf_orf = map_orfs.swifter.apply(lambda x: (x.orf,x.gene) if \
	(x.gene in caf_genes["Gene Systematic Name"].values) else (None,None), axis=1)
caf_orf = pd.DataFrame(caf_orf.swifter.apply(pd.Series)).dropna()
cu_orf = map_orfs.swifter.apply(lambda x: (x.orf,x.gene) if \
	(x.gene in cu_genes["Gene Systematic Name"].values) else (None,None), axis=1)
cu_orf = pd.DataFrame(cu_orf.swifter.apply(pd.Series)).dropna()
sma_orf = map_orfs.swifter.apply(lambda x: (x.orf,x.gene) if \
	(x.gene in sma_genes["Gene Systematic Name"].values) else (None,None), axis=1)
sma_orf = pd.DataFrame(sma_orf.swifter.apply(pd.Series)).dropna()

# Now create the feature lists of top features + individual literature gene lists
d2 = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB"
models_to_run_a = open(f"{d}/fs_plus_lit_genes_models_to_run.txt", "a") # list of models to run
models_to_run_b = open(f"{d}/fs_plus_unimportant_non_lit_genes_models_to_run.txt", "a")
for data_type in ["snp", "pav", "cnv"]:
	for env in target_envs:
		for lit_env in lit_envs:
			print(data_type, env, lit_env)
			#
			# Determine which model has the top score
			res_env = res.loc[(res.Trait.str.contains(f"{env}") & \
				res.Data.str.contains(f"{data_type}")),:]
			top = res_env.loc[res_env["R2_val            "] == res_env["R2_val            "].max(),"FeatureNum            "].values[0]
			#
			# Read in the top model feature importances
			mod = pickle.load(open(f"{d}/{data_type}_{env}_top_{top}_model_rep_0.pkl", "rb"))
			imp = pd.read_csv(f"{d}/{data_type}_{env}_top_{top}_imp.csv", index_col=0)
			imp.index = mod.feature_names_in_
			imp["mean_imp"] = imp.mean(axis=1)
			imp = imp.loc[imp.mean_imp != 0.0,:] # drop unimportant features
			#
			# Get literature genes that map to baseline model
			imp_base = pd.read_csv(f"{d2}/{data_type}_{env}_imp.csv", index_col=0)
			mod_base = pickle.load(open(f"{d2}/{data_type}_{env}_model_rep_0.pkl", "rb"))
			imp_base.index = mod_base.feature_names_in_
			imp_base["mean_imp"] = imp_base.mean(axis=1)
			#
			if data_type != "snp":
				imp.index = imp.apply(lambda x: re.sub("^X", "", x.name), axis=1)
				imp.index = imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
				imp_base.index = imp_base.apply(lambda x: re.sub("^X", "", x.name), axis=1)
				imp_base.index = imp_base.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
			#
			if lit_env == "benomyl":
				if data_type=="snp":
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(ben_snp.iloc[:,0]),:] # snp features
					T_l = T_l["mean_imp"].to_frame().merge(ben_snp, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "snp", "gene"]
					T_nl = imp.loc[~imp.index.isin(ben_snp.iloc[:,0]),:] # snp features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					#
					# Identify the unimportant model features found in ben_snp and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(ben_snp.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(ben_snp.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_b = pd.concat([T_l, T_nl, B_nl])
						
				else:
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(ben_orf.iloc[:,0]),:] # orf features
					T_l = T_l["mean_imp"].to_frame().merge(ben_orf, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "orf", "gene"]
					T_nl = imp.loc[~imp.index.isin(ben_orf.iloc[:,0]),:] # orf features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					#
					# Identify the unimportant model features found in ben_orf and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(ben_orf.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(ben_orf.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_b = pd.concat([T_l, T_nl, B_nl])
					#
					# Replace missing gene values with the orf ID
					features_a.gene = features_a.gene.fillna(features_a.orf)
					features_b.gene = features_b.gene.fillna(features_b.orf)
				#
				# Ensure that the genes are unique
				if (features_a.gene.nunique() == len(features_a)) & \
					(features_b.gene.nunique() == len(features_b)):
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_benomyl_lit_genes_info", features_a, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_benomyl_lit_genes", features_a.iloc[:,1], fmt="%s")
					models_to_run_a.write(f"{d}/{data_type}_{env}_top_{top}_plus_benomyl_lit_genes\n")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_benomyl_genes_info", features_b, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_benomyl_genes", features_b.iloc[:,1], fmt="%s")
					models_to_run_b.write(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_benomyl_genes\n")
				#
			elif lit_env == "caffeine":
				if data_type=="snp":
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(caf_snp.iloc[:,0]),:] # snp features
					T_l = T_l["mean_imp"].to_frame().merge(caf_snp, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "snp", "gene"]
					T_nl = imp.loc[~imp.index.isin(caf_snp.iloc[:,0]),:] # snp features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					#
					# Identify the unimportant model features found in caf_snp and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(caf_snp.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(caf_snp.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_b = pd.concat([T_l, T_nl, B_nl])
						
				else:
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(caf_orf.iloc[:,0]),:] # orf features
					T_l = T_l["mean_imp"].to_frame().merge(caf_orf, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "orf", "gene"]
					T_nl = imp.loc[~imp.index.isin(caf_orf.iloc[:,0]),:] # orf features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					#
					# Identify the unimportant model features found in caf_orf and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(caf_orf.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(caf_orf.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_b = pd.concat([T_l, T_nl, B_nl])
					#
					# Replace missing gene values with the orf ID
					features_a.gene = features_a.gene.fillna(features_a.orf)
					features_b.gene = features_b.gene.fillna(features_b.orf)
				#
				# Ensure that the genes are unique
				if (features_a.gene.nunique() == len(features_a)) & \
					(features_b.gene.nunique() == len(features_b)):
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_caffeine_lit_genes_info", features_a, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_caffeine_lit_genes", features_a.iloc[:,1], fmt="%s")
					models_to_run_a.write(f"{d}/{data_type}_{env}_top_{top}_plus_caffeine_lit_genes\n")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_caffeine_genes_info", features_b, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_caffeine_genes", features_b.iloc[:,1], fmt="%s")
					models_to_run_b.write(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_caffeine_genes\n")
				#
			elif lit_env == "cuso4":
				if data_type=="snp":
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(cu_snp.iloc[:,0]),:] # snp features
					T_l = T_l["mean_imp"].to_frame().merge(cu_snp, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "snp", "gene"]
					T_nl = imp.loc[~imp.index.isin(cu_snp.iloc[:,0]),:] # snp features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					#
					# Identify the unimportant model features found in cu_snp and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(cu_snp.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(cu_snp.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_b = pd.concat([T_l, T_nl, B_nl])
						
				else:
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(cu_orf.iloc[:,0]),:] # orf features
					T_l = T_l["mean_imp"].to_frame().merge(cu_orf, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "orf", "gene"]
					T_nl = imp.loc[~imp.index.isin(cu_orf.iloc[:,0]),:] # orf features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					#
					# Identify the unimportant model features found in cu_orf and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(cu_orf.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(cu_orf.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_b = pd.concat([T_l, T_nl, B_nl])
					#
					# Replace missing gene values with the orf ID
					features_a.gene = features_a.gene.fillna(features_a.orf)
					features_b.gene = features_b.gene.fillna(features_b.orf)
				#
				# Ensure that the genes are unique
				if (features_a.gene.nunique() == len(features_a)) & \
					(features_b.gene.nunique() == len(features_b)):
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_cuso4_lit_genes_info", features_a, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_cuso4_lit_genes", features_a.iloc[:,1], fmt="%s")
					models_to_run_a.write(f"{d}/{data_type}_{env}_top_{top}_plus_cuso4_lit_genes\n")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_cuso4_genes_info", features_b, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_cuso4_genes", features_b.iloc[:,1], fmt="%s")
					models_to_run_b.write(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_cuso4_genes\n")
				#
			elif lit_env == "sma":
				if data_type=="snp":
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(sma_snp.iloc[:,0]),:] # snp features
					T_l = T_l["mean_imp"].to_frame().merge(sma_snp, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "snp", "gene"]
					T_nl = imp.loc[~imp.index.isin(sma_snp.iloc[:,0]),:] # snp features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					#
					# Identify the unimportant model features found in sma_snp and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(sma_snp.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(sma_snp.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
					features_b = pd.concat([T_l, T_nl, B_nl])
						
				else:
					# Determine literature and non-literature genes in the top features
					T_l = imp.loc[imp.index.isin(sma_orf.iloc[:,0]),:] # orf features
					T_l = T_l["mean_imp"].to_frame().merge(sma_orf, how="left", left_index=True, right_on=0)
					T_l.columns = ["mean_imp", "orf", "gene"]
					T_nl = imp.loc[~imp.index.isin(sma_orf.iloc[:,0]),:] # orf features
					T_nl = T_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					#
					# Identify the unimportant model features found in sma_orf and combine with imp
					B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(imp_base.index.isin(sma_orf.iloc[:,0])),:]
					B_l = B_l["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_a = pd.concat([T_l, T_nl, B_l])
					B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
						(~imp_base.index.isin(sma_orf.iloc[:,0])),:]
					B_nl = B_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
					features_b = pd.concat([T_l, T_nl, B_nl])
					#
					# Replace missing gene values with the orf ID
					features_a.gene = features_a.gene.fillna(features_a.orf)
					features_b.gene = features_b.gene.fillna(features_b.orf)
				#
				# Ensure that the genes are unique
				if (features_a.gene.nunique() == len(features_a)) & \
					(features_b.gene.nunique() == len(features_b)):
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_sodium_meta-arsenite_lit_genes_info", features_a, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_sodium_meta-arsenite_lit_genes", features_a.iloc[:,1], fmt="%s")
					models_to_run_a.write(f"{d}/{data_type}_{env}_top_{top}_plus_sodium_meta-arsenite_lit_genes\n")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_sodium_meta-arsenite_genes_info", features_b, fmt="%s")
					np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_sodium_meta-arsenite_genes", features_b.iloc[:,1], fmt="%s")
					models_to_run_b.write(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_sodium_meta-arsenite_genes\n")
				#
			del res_env, top, imp, imp_base, mod, mod_base

models_to_run_a.close()
models_to_run_b.close()

# Now create the feature lists of top features + combined literature gene list
models_to_run = open(f"{d}/fs_plus_comb_lit_genes_models_to_run.txt", "a") # list of models to run
for data_type in ["snp", "pav", "cnv"]:
	for env in target_envs:
		print(data_type, env)
		# Determine which model has the top score
		res_env = res.loc[(res.Trait.str.contains(f"{env}") & \
			res.Data.str.contains(f"{data_type}")),:]
		top = res_env.loc[res_env["R2_val            "] == res_env["R2_val            "].max(),"FeatureNum            "].values[0]
		# Get the mean feature importances
		mod = pickle.load(open(f"{d}/{data_type}_{env}_top_{top}_model_rep_0.pkl", "rb"))
		imp = pd.read_csv(f"{d}/{data_type}_{env}_top_{top}_imp.csv", index_col=0)
		imp.index = mod.feature_names_in_
		imp["mean_imp"] = imp.mean(axis=1)
		imp = imp.loc[imp.mean_imp != 0.0,:] # drop unimportant features
		# Get literature genes that map to baseline model
		imp_base = pd.read_csv(f"{d2}/{data_type}_{env}_imp.csv", index_col=0)
		mod_base = pickle.load(open(f"{d2}/{data_type}_{env}_model_rep_0.pkl", "rb"))
		imp_base.index = mod_base.feature_names_in_
		if data_type=="snp":
			# Identify the baseline model features found in ben_snp and combine with imp
			lit_ben_snp = imp_base.loc[imp_base.index.isin(ben_snp.iloc[:,0]),:]
			lit_caf_snp = imp_base.loc[imp_base.index.isin(caf_snp.iloc[:,0]),:]
			lit_cu_snp = imp_base.loc[imp_base.index.isin(cu_snp.iloc[:,0]),:]
			lit_sma_snp = imp_base.loc[imp_base.index.isin(sma_snp.iloc[:,0]),:]
			features = pd.concat([imp, lit_ben_snp, lit_caf_snp, lit_cu_snp,\
								lit_sma_snp]).index.unique()
		else:
			lit_ben_orf = imp_base.loc[imp_base.index.isin(ben_orf.iloc[:,0]),:]
			lit_caf_orf = imp_base.loc[imp_base.index.isin(caf_orf.iloc[:,0]),:]
			lit_cu_orf = imp_base.loc[imp_base.index.isin(cu_orf.iloc[:,0]),:]
			lit_sma_orf = imp_base.loc[imp_base.index.isin(sma_orf.iloc[:,0]),:]
			features = pd.concat([imp, lit_ben_orf, lit_caf_orf, lit_cu_orf,\
								lit_sma_orf]).index.unique()
		np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_comb_lit_genes", features, fmt="%s")
		models_to_run.write(f"{d}/{data_type}_{env}_top_{top}_plus_comb_lit_genes\n")
		del res_env, top, imp, imp_base, mod, mod_base
models_to_run.close()

# Randomize the label and add it to the top features + combined literature gene list datasets
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
pheno_train = pheno.loc[~pheno.index.isin(test[0]),:]
dr = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/randomized"
for i in range(100):
	randomized = pd.DataFrame(index=pheno_train.index)
	for env in target_envs:
		randomized[env] = pheno_train[env].sample(frac=1, random_state=i).reset_index(drop=True).values
	randomized = pd.concat([randomized, pheno.loc[pheno.index.isin(test[0]),target_envs]]) # add the test set back
	randomized.to_csv(f"{dr}/pheno_randomized_{i}.csv")
	del randomized

########################### Model performance figure ###########################
## These are the 15 models with the FS features and all lit genes combined
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/fs"
res = pd.read_csv(f"{d}/RESULTS_xgboost.txt", sep="\t") # Model results
res.rename(columns={"R2_test            ":"R2_test"}, inplace=True)
snp = res.loc[res.Data=="snp", ["Trait", "R2_test", "R2_test_sd"]].sort_values("Trait")
pav = res.loc[res.Data=="pav", ["Trait", "R2_test", "R2_test_sd"]].sort_values("Trait")
cnv = res.loc[res.Data=="cnv", ["Trait", "R2_test", "R2_test_sd"]].sort_values("Trait")

# Randomized model results
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
res_rand = pd.DataFrame() # Randomized label model results
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None) # test set instances
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # label
pheno_test = pheno.loc[test[0],:]
dr = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction/XGB/randomized"
for env in target_envs:
	for data in ["snp", "pav", "cnv"]:
		print(env, data)
		scores = []
		scores_sd = []
		for i in range(100): # 100 randomizations per environment
			try: # to ensure no file is missing
				file = glob.glob(f"{dr}/{data}_{env}_top_*_plus_comb_lit_genes_rand_{i}_preds.csv")[0]
				preds = pd.read_csv(file, index_col=0) # predicted labels
				preds = preds.loc[test[0],:]
				r2_test = preds.apply(lambda x: r2_score(pheno_test.loc[x.index,env], x))
				scores.append(r2_test.mean()) # Average R2_test across reps
				scores_sd.append(r2_test.std())
				del file, preds, r2_test
			except IndexError:
				print(f"Missing {data}_{env}_top_*_plus_comb_lit_genes_rand_{i}_preds.csv")
		try:
			res_rand = res_rand._append(pd.Series([data, env, np.mean(scores), np.mean(scores_sd)]), ignore_index=True)
			del scores, scores_sd
		except:
			pass

res_rand.columns = ["Data", "Trait", "R2_test", "R2_test_sd"]
res_rand.to_csv(f"{dr}/RESULTS_xgboost_summary.txt", index=False)

# Plot model performances
fig, ax = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(8,3.8))
sns.barplot(snp, x="Trait", y="R2_test", errorbar="sd", yerr=snp["R2_test_sd"], ax=ax[0])
sns.barplot(pav, x="Trait", y="R2_test", errorbar="sd", yerr=snp["R2_test_sd"], ax=ax[1])
sns.barplot(cnv, x="Trait", y="R2_test", errorbar="sd", yerr=snp["R2_test_sd"], ax=ax[2])
ax[0].set_xticklabels(snp["Trait"], size=7, rotation=50, ha="right")
ax[1].set_xticklabels(snp["Trait"], size=7, rotation=50, ha="right")
ax[2].set_xticklabels(snp["Trait"], size=7, rotation=50, ha="right")
snp_rand = res_rand.loc[res_rand.Data=="snp",:].set_index('Trait')
pav_rand = res_rand.loc[res_rand.Data=="pav",:].set_index('Trait')
cnv_rand = res_rand.loc[res_rand.Data=="cnv",:].set_index('Trait')
sns.scatterplot(snp, x="Trait", y=snp_rand.loc[snp.Trait,"R2_test"].values,
	marker="o", color="black", ax=ax[0]) # Add randomization models average R-sq
ax[0].errorbar(x=snp.Trait, y=snp_rand.loc[snp.Trait,"R2_test"].values,
	yerr=snp_rand.loc[snp.Trait,"R2_test_sd"].values)
sns.scatterplot(pav, x="Trait", y=pav_rand.loc[pav.Trait,"R2_test"].values,
	marker="o", color="black", ax=ax[1])
ax[1].errorbar(x=pav.Trait, y=pav_rand.loc[pav.Trait,"R2_test"].values,
	yerr=pav_rand.loc[pav.Trait,"R2_test_sd"].values)
sns.scatterplot(cnv, x="Trait", y=cnv_rand.loc[cnv.Trait,"R2_test"].values,
	marker="o", color="black", ax=ax[2])
ax[2].errorbar(cnv.Trait, y=cnv_rand.loc[cnv.Trait,"R2_test"].values,
	yerr=cnv_rand.loc[cnv.Trait,"R2_test_sd"].values)
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/XGB_comb_lit_gene_results.pdf")
plt.close()

################################################################################
### TABLE S12
################################################################################
# For each data type:
		# total number of unique genes in gene1 and gene2 columns
			# overlap in gene1 and gene2
		# total number of interactions
		# how many literature genes were identified
			# how many interactions were literature genes found to have
			# how many of the known GIs of these literature genes were identified
		# are there "novel" interactions that make sense
	# Comparison of data types:
		# how many interactions overlap between data types

## Map summed SHAP interaction scores to genes
# Feature to gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
					   header=None, names=["snp", "chr", "pos", "gene"]) # snp to gene map
map_snps.set_index("snp", inplace=True)
map_snps_dict = map_snps["gene"].to_dict()

map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t") # orf to gene map
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)
map_orfs.set_index("orf", inplace=True)
map_orfs_dict = map_orfs["gene"].to_dict()

# Literature gene lists
d = "Data/SGD_Experiment_Genes"
beno = pd.read_csv(f"{d}/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf = dt.fread(f"{d}/caffeine_phenotype_annotations_sensitive_genes.txt").to_pandas()
cu = pd.read_csv(f"{d}/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma = pd.read_csv(f"{d}/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")

# Add whether the SNPs and ORFs are literature genes to map files
map_snps["benomyl_gene"] = map_snps.apply(lambda x: 1 if \
	x.gene in beno["Gene Systematic Name"].values else 0, axis=1)
map_snps["caffeine_gene"] = map_snps.apply(lambda x: 1 if \
	x.gene in caf["Gene Systematic Name"].values else 0, axis=1)
map_snps["cuso4_gene"] = map_snps.apply(lambda x: 1 if \
	x.gene in cu["Gene Systematic Name"].values else 0, axis=1)
map_snps["sodium_meta-arsenite_gene"] = map_snps.apply(lambda x: 1 if \
	x.gene in sma["Gene Systematic Name"].values else 0, axis=1)
map_orfs["benomyl_gene"] = map_orfs.apply(lambda x: 1 if \
	x.gene in beno["Gene Systematic Name"].values else 0, axis=1)
map_orfs["caffeine_gene"] = map_orfs.apply(lambda x: 1 if \
	x.gene in caf["Gene Systematic Name"].values else 0, axis=1)
map_orfs["cuso4_gene"] = map_orfs.apply(lambda x: 1 if \
	x.gene in cu["Gene Systematic Name"].values else 0, axis=1)
map_orfs["sodium_meta-arsenite_gene"] = map_orfs.apply(lambda x: 1 if \
	x.gene in sma["Gene Systematic Name"].values else 0, axis=1)
map_snps.to_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_with_lit.txt",
	sep="\t")
map_orfs.to_csv("Data/Peter_2018/final_map_orf_to_gene_with_lit.txt", sep="\t")

################################################################################
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction"
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_with_lit.txt", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_with_lit.txt", sep="\t")
map_snps.drop(columns=["chr", "pos"], inplace=True)
map_orfs.drop(columns="organism", inplace=True)

################################# SHAP Values ##################################
top_features = {'snp':{}, 'pav':{}, 'cnv':{}} # for scatterplots in next section
kinship = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/kinship.csv", index_col=0)
test = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", sep="\t", header=None)
for env in target_envs:
	print(env)
	snp_shap = dt.fread(glob.glob(f"{d}/SNP/SHAP_values_sorted_average_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap = dt.fread(glob.glob(f"{d}/PAV/SHAP_values_sorted_average_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	cnv_shap = dt.fread(glob.glob(f"{d}/CNV/SHAP_values_sorted_average_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	# map to genes and lit genes
	snp_shap.rename(columns={"C0":"snp", "0":"mean_shap"}, inplace=True)
	pav_shap.rename(columns={"C0":"orf", "0":"mean_shap"}, inplace=True)
	cnv_shap.rename(columns={"C0":"orf", "0":"mean_shap"}, inplace=True)
	snp_shap = snp_shap.loc[snp_shap.mean_shap != 0.0,:] # drop unimportant features
	pav_shap = pav_shap.loc[pav_shap.mean_shap != 0.0,:]
	cnv_shap = cnv_shap.loc[cnv_shap.mean_shap != 0.0,:]
	pav_shap["orf"] = pav_shap.apply(lambda x: re.sub("^X", "", x.orf), axis=1) # fix orf IDs
	pav_shap["orf"] = pav_shap.apply(lambda x: re.sub("\.", "-", x.orf), axis=1)
	cnv_shap["orf"] = cnv_shap.apply(lambda x: re.sub("^X", "", x.orf), axis=1)
	cnv_shap["orf"] = cnv_shap.apply(lambda x: re.sub("\.", "-", x.orf), axis=1)
	snp_shap = snp_shap.set_index("snp").merge(map_snps.set_index("snp"), how="left", left_index=True, right_index=True) # mapping
	pav_shap = pav_shap.set_index("orf").merge(map_orfs.set_index("orf"), how="left", left_index=True, right_index=True)
	cnv_shap = cnv_shap.set_index("orf").merge(map_orfs.set_index("orf"), how="left", left_index=True, right_index=True)
	snp_shap["abs_mean_shap"] = snp_shap.mean_shap.abs() # absolute value of mean shap value
	pav_shap["abs_mean_shap"] = pav_shap.mean_shap.abs()
	cnv_shap["abs_mean_shap"] = cnv_shap.mean_shap.abs()
	snp_shap.sort_values("abs_mean_shap", ascending=False, inplace=True) # sort shap values
	pav_shap.sort_values("abs_mean_shap", ascending=False, inplace=True)
	cnv_shap.sort_values("abs_mean_shap", ascending=False, inplace=True)
	top_features['snp'][env] = snp_shap.index[0]
	top_features['pav'][env] = pav_shap.index[0]
	top_features['cnv'][env] = cnv_shap.index[0]
	# fill in NaNs in gene column with orf name
	pav_shap["gene"] = pav_shap.apply(lambda x: x.name if x.gene is np.nan else x.gene, axis=1)
	cnv_shap["gene"] = cnv_shap.apply(lambda x: x.name if x.gene is np.nan else x.gene, axis=1)
	pav_shap.set_index("gene", inplace=True)
	cnv_shap.set_index("gene", inplace=True)
	snp_shap.set_index("gene", inplace=True)
	# create a heatmap of shap values
	pav_shap.fillna(0, inplace=True)
	cnv_shap.fillna(0, inplace=True)
	# fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6,10)) ## SNP Figure
	# sns.heatmap(snp_shap.iloc[:20,:1], ax=ax[0], square=True, cbar_kws={"orientation":"horizontal"})
	# sns.heatmap(snp_shap.iloc[:20,1:], ax=ax[1], vmin=0, vmax=1, square=True, cbar_kws={"orientation":"horizontal"})
	# plt.tight_layout()
	# plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_snp_{env}_top_20_plus_comb_lit_genes_training.pdf")
	# plt.close()
	# fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6, 10)) ## PAV Figure
	# sns.heatmap(pav_shap.iloc[:20,:1], ax=ax[0], square=True, cbar_kws={"orientation":"horizontal"})
	# sns.heatmap(pav_shap.iloc[:20,1:], ax=ax[1], vmin=0, vmax=1, square=True, cbar_kws={"orientation":"horizontal"})
	# plt.tight_layout()
	# plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_pav_{env}_top_20_plus_comb_lit_genes_training.pdf")
	# plt.close()
	# fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6, 10)) ## CNV Figure
	# sns.heatmap(cnv_shap.iloc[:20,:1], ax=ax[0], square=True, cbar_kws={"orientation":"horizontal"})
	# sns.heatmap(cnv_shap.iloc[:20,1:], ax=ax[1], vmin=0, vmax=1, square=True, cbar_kws={"orientation":"horizontal"})
	# plt.tight_layout()
	# plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_cnv_{env}_top_20_plus_comb_lit_genes_training.pdf")
	# plt.close()
	# Plot the SHAP value of the top genes vs kinship with W303
	snp_shap_all = pd.read_csv(glob.glob(f"{d}/SNP/SHAP_values_sorted_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
	pav_shap_all = pd.read_csv(glob.glob(f"{d}/PAV/SHAP_values_sorted_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
	cnv_shap_all = pd.read_csv(glob.glob(f"{d}/CNV/SHAP_values_sorted_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
	pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix orf IDs
	pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	cnv_shap_all.columns = cnv_shap_all.apply(lambda x: re.sub("^X", "", x.name), axis=0)
	cnv_shap_all.columns = cnv_shap_all.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	x = kinship.loc[~kinship.index.isin(test[0]), "SACE_GAV"]
	fig, ax = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(8.5,3.5))
	sns.regplot(x=x, y=snp_shap_all.loc[x.index, top_features['snp'][env]], ci=None, ax=ax[0])
	m, b, r, p, sd = stats.linregress(x, snp_shap_all.loc[x.index, top_features['snp'][env]])
	ax[0].annotate(f"Slope: {m:.2f}\nIntercept: {b:.2f}\nPCC: {r:.2f}\nP-value: {p:.2f}\nStd. Err.: {sd:.2f}",
		xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
	ax[0].set_xlabel("Kinship with W303")
	ax[0].set_ylabel(f"SHAP value of top gene {top_features['snp'][env]}")
	#
	sns.regplot(x=x, y=pav_shap_all.loc[x.index, top_features['pav'][env]], ci=None, ax=ax[1])
	m, b, r, p, sd = stats.linregress(x, pav_shap_all.loc[x.index, top_features['pav'][env]])
	ax[1].annotate(f"Slope: {m:.2f}\nIntercept: {b:.2f}\nPCC: {r:.2f}\nP-value: {p:.2f}\nStd. Err.: {sd:.2f}",
		xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
	ax[1].set_xlabel("Kinship with W303")
	ax[1].set_ylabel(f"SHAP value of top gene {top_features['pav'][env]}")
	#
	sns.regplot(x=x, y=cnv_shap_all.loc[x.index, top_features['cnv'][env]], ci=None, ax=ax[2])
	m, b, r, p, sd = stats.linregress(x, cnv_shap_all.loc[x.index, top_features['cnv'][env]])
	ax[2].annotate(f"Slope: {m:.2f}\nIntercept: {b:.2f}\nPCC: {r:.2f}\nP-value: {p:.2f}\nStd. Err.: {sd:.2f}",
		xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
	ax[2].set_xlabel("Kinship with W303")
	ax[2].set_ylabel(f"SHAP value of top gene {top_features['cnv'][env]}")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/{env}_top_gene_shap_vs_w303_kinship.pdf")
	plt.close()

# plot with the top gene for w303
	snp_shap = dt.fread(glob.glob(f"{d}/SNP/SHAP_values_sorted_average_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap = dt.fread(glob.glob(f"{d}/PAV/SHAP_values_sorted_average_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	cnv_shap = dt.fread(glob.glob(f"{d}/CNV/SHAP_values_sorted_average_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()

pd.DataFrame.from_dict(top_features).to_csv("Scripts/Data_Vis/Section_5/shap_value_top_features_per_data_type.csv")

###################### SHAP value trends of select genes #######################
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)
geno_snp = dt.fread("Data/Peter_2018/geno.csv").to_pandas().set_index("ID")
geno_pav = pd.read_csv("Data/Peter_2018/ORFs_pres_abs.csv", index_col=0)
geno_cnv = pd.read_csv("Data/Peter_2018/ORFs_no_NA.csv", index_col=0)
for env in target_envs:
	snp_shap = dt.fread(glob.glob(f"{d}/SNP/SHAP_values_sorted_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap = dt.fread(glob.glob(f"{d}/PAV/SHAP_values_sorted_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	cnv_shap = dt.fread(glob.glob(f"{d}/CNV/SHAP_values_sorted_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap.columns = pav_shap.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix orf IDs
	pav_shap.columns = pav_shap.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	cnv_shap.columns = cnv_shap.apply(lambda x: re.sub("^X", "", x.name), axis=0)
	cnv_shap.columns = cnv_shap.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	# calculate average and get the feature with the max average shap
	snp_top_most_feature = top_features['snp'][env]
	pav_top_most_feature = top_features['pav'][env]
	cnv_top_most_feature = top_features['cnv'][env]
	pav_top_most_feature2 = "X" + re.sub("-", ".", pav_top_most_feature)
	cnv_top_most_feature2 = "X" + re.sub("-", ".", cnv_top_most_feature)
	# plot data type values vs fitness scatterplot colored by shap value
	plt.figure(figsize=(3,3))
	sns.scatterplot(x=geno_snp.loc[snp_shap.ID, snp_top_most_feature].values,\
					y=pheno.loc[snp_shap.ID, env].values,\
					hue=snp_shap.loc[:,snp_top_most_feature].values,\
					palette="RdBu", alpha=0.5)
	plt.title(map_snps.loc[map_snps.snp==snp_top_most_feature,'gene'].values[0])
	plt.legend(frameon=False)
	plt.xlabel("Genotype")
	plt.ylabel("Fitness")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_snp_{env}_top_1_scatterplot.pdf") ## SNP Figure
	plt.close()
	plt.figure(figsize=(3,3))
	sns.scatterplot(x=geno_pav.loc[pav_shap.ID, pav_top_most_feature2].values,\
					y=pheno.loc[pav_shap.ID, env].values,\
					hue=pav_shap.loc[:,pav_top_most_feature].values,\
					palette="RdBu", alpha=0.5)
	try:
		plt.title(map_orfs.loc[map_orfs.orf==pav_top_most_feature,'gene'].values[0])
	except:
		plt.title(pav_top_most_feature)
	plt.legend(frameon=False)
	plt.xlabel("Presence (1) or Absence (0)")
	plt.ylabel("Fitness")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_pav_{env}_top_1_scatterplot.pdf") ## PAV Figure
	plt.close()
	plt.figure(figsize=(3,3))
	sns.scatterplot(x=geno_cnv.loc[cnv_shap.ID, cnv_top_most_feature2].values,\
					y=pheno.loc[cnv_shap.ID, env].values,\
					hue=cnv_shap.loc[:,cnv_top_most_feature].values,\
					palette="RdBu", alpha=0.5)
	try:
		plt.title(map_orfs.loc[map_orfs.orf==cnv_top_most_feature,'gene'].values[0])
	except:
		plt.title(cnv_top_most_feature)
	plt.legend(frameon=False)
	plt.xlabel("Copy number")
	plt.ylabel("Fitness")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_cnv_{env}_top_1_scatterplot.pdf") ## CNV Figure
	plt.close()

########################### SHAP Interaction Scores ############################
# Violin plots of interactions cores
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction"
for env in target_envs:	
	print(env)
	try:
		snp = dt.fread(glob.glob(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	except:
		pass # missing benomyl snp summed file
	pav = dt.fread(glob.glob(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	cnv = dt.fread(glob.glob(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	# violin plot of interaction scores
	try:
		sns.violinplot(pd.DataFrame({"snp":snp.Interaction, "pav":pav.Interaction, "cnv":cnv.Interaction}), fill=False)
	except:
		sns.violinplot(pd.DataFrame({"pav":pav.Interaction, "cnv":cnv.Interaction}), fill=False)
	plt.ylim(0,2) # for YPDSODIUMMETAARSENITE
	plt.ylabel("SHAP Interaction Scores")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_interaction_{env}_violin.pdf") #_ylim.pdf for YPDSODIUMMETAARSENITE
	plt.close()
	try: # for violin plot debugging
		del snp, pav, cnv
	except:
		del pav, cnv

# Count the number of experimentally verified lit genes were identified by SHAP
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt")
map_snps.columns = ["snp", "chr", "pos", "gene"]
map_snps.set_index("snp", inplace=True)
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs.set_index("orf", inplace=True)
for env in target_envs:
	try:
		snp = dt.fread(glob.glob(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	except:
		pass # missing benomyl snp summed file
	pav = dt.fread(glob.glob(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	cnv = dt.fread(glob.glob(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	## map feature names to genes
	try:
		snp["Gene1"] = snp.apply(lambda x: map_snps["gene"][x.Feature1], axis=1)
		snp["Gene2"] = snp.apply(lambda x: map_snps["gene"][x.Feature2], axis=1)
	except:
		pass
	# print(snp[["Gene1", "Gene2", "Interaction"]].groupby(["Gene1", "Gene2"]).max().shape)
	# print('snp genes', snp.Gene1.nunique(), snp.Gene2.nunique())
	#
	pav["Feature1"] = pav.apply(lambda x: re.sub("^X", "", x.Feature1), axis=1)
	pav["Feature1"] = pav.apply(lambda x: re.sub("\.", "-", x.Feature1), axis=1)
	pav["Feature2"] = pav.apply(lambda x: re.sub("^X", "", x.Feature2), axis=1)
	pav["Feature2"] = pav.apply(lambda x: re.sub("\.", "-", x.Feature2), axis=1)
	pav["Gene1"] = pav.apply(lambda x: map_orfs["gene"][x.Feature1] \
							if x.Feature1 in map_orfs["gene"].keys() else x.Feature1, axis=1)
	pav["Gene2"] = pav.apply(lambda x: map_orfs["gene"][x.Feature2] \
							if x.Feature2 in map_orfs["gene"].keys() else x.Feature2, axis=1)
	# print(pav[["Gene1", "Gene2", "Interaction"]].groupby(["Gene1", "Gene2"]).max().shape)
	# print('gene1 genes', pd.Series(pav.Gene1.unique()).isin(map_orfs["gene"].values()).sum())
	# print('gene2 genes', pd.Series(pav.Gene2.unique()).isin(map_orfs["gene"].values()).sum())
	# print('gene1 orfs', pav.Gene1.nunique() - pd.Series(pav.Gene1.unique()).isin(map_orfs["gene"].values()).sum())
	# print('gene2 orfs', pav.Gene2.nunique() - pd.Series(pav.Gene2.unique()).isin(map_orfs["gene"].values()).sum())
	cnv["Feature1"] = cnv.apply(lambda x: re.sub("^X", "", x.Feature1), axis=1)
	cnv["Feature1"] = cnv.apply(lambda x: re.sub("\.", "-", x.Feature1), axis=1)
	cnv["Feature2"] = cnv.apply(lambda x: re.sub("^X", "", x.Feature2), axis=1)
	cnv["Feature2"] = cnv.apply(lambda x: re.sub("\.", "-", x.Feature2), axis=1)
	cnv["Gene1"] = cnv.apply(lambda x: map_orfs["gene"][x.Feature1] \
							if x.Feature1 in map_orfs["gene"].keys() else x.Feature1, axis=1)
	cnv["Gene2"] = cnv.apply(lambda x: map_orfs["gene"][x.Feature2] \
							if x.Feature2 in map_orfs["gene"].keys() else x.Feature2, axis=1)
	#
	## number of literature genes identified
	try:
		snp["Gene1_benomyl"] = snp.apply(lambda x: 1 if x.Gene1 in beno["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_benomyl"] = snp.apply(lambda x: 1 if x.Gene2 in beno["Gene Systematic Name"] else 0, axis=1)
		snp["Gene1_caffeine"] = snp.apply(lambda x: 1 if x.Gene1 in caf["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_caffeine"] = snp.apply(lambda x: 1 if x.Gene2 in caf["Gene Systematic Name"] else 0, axis=1)
		snp["Gene1_cuso4"] = snp.apply(lambda x: 1 if x.Gene1 in cu["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_cuso4"] = snp.apply(lambda x: 1 if x.Gene2 in cu["Gene Systematic Name"] else 0, axis=1)
		snp["Gene1_sma"] = snp.apply(lambda x: 1 if x.Gene1 in sma["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_sma"] = snp.apply(lambda x: 1 if x.Gene2 in sma["Gene Systematic Name"] else 0, axis=1)
		print('snp gene1 beno', pd.Series(snp.Gene1.unique()).isin(beno["Gene Systematic Name"]).sum(),\
			'snp gene2 beno', pd.Series(snp.Gene2.unique()).isin(beno["Gene Systematic Name"]).sum())
		print('snp gene1 caf', pd.Series(snp.Gene1.unique()).isin(caf["Gene Systematic Name"]).sum(),\
			'snp gene2 caf', pd.Series(snp.Gene2.unique()).isin(caf["Gene Systematic Name"]).sum())
		print('snp gene1 cu', pd.Series(snp.Gene1.unique()).isin(cu["Gene Systematic Name"]).sum(),\
			'snp gene2 cu', pd.Series(snp.Gene2.unique()).isin(cu["Gene Systematic Name"]).sum())
		print('snp gene1 sma', pd.Series(snp.Gene1.unique()).isin(sma["Gene Systematic Name"]).sum(),\
			'snp gene2 sma', pd.Series(snp.Gene2.unique()).isin(sma["Gene Systematic Name"]).sum())
	except:
		pass
	#
	pav["Gene1_benomyl"] = pav.apply(lambda x: 1 if x.Gene1 in beno["Gene Systematic Name"] else 0, axis=1) 
	pav["Gene2_benomyl"] = pav.apply(lambda x: 1 if x.Gene2 in beno["Gene Systematic Name"] else 0, axis=1)
	pav["Gene1_caffeine"] = pav.apply(lambda x: 1 if x.Gene1 in caf["Gene Systematic Name"] else 0, axis=1) 
	pav["Gene2_caffeine"] = pav.apply(lambda x: 1 if x.Gene2 in caf["Gene Systematic Name"] else 0, axis=1)
	pav["Gene1_cuso4"] = pav.apply(lambda x: 1 if x.Gene1 in cu["Gene Systematic Name"] else 0, axis=1) 
	pav["Gene2_cuso4"] = pav.apply(lambda x: 1 if x.Gene2 in cu["Gene Systematic Name"] else 0, axis=1)
	pav["Gene1_sma"] = pav.apply(lambda x: 1 if x.Gene1 in sma["Gene Systematic Name"] else 0, axis=1)
	pav["Gene2_sma"] = pav.apply(lambda x: 1 if x.Gene2 in sma["Gene Systematic Name"] else 0, axis=1)
	print('pav gene1 beno', pd.Series(pav.Gene1.unique()).isin(beno["Gene Systematic Name"]).sum(),\
		'pav gene2 beno', pd.Series(pav.Gene2.unique()).isin(beno["Gene Systematic Name"]).sum())
	print('pav gene1 caf', pd.Series(pav.Gene1.unique()).isin(caf["Gene Systematic Name"]).sum(),\
		'pav gene2 caf', pd.Series(pav.Gene2.unique()).isin(caf["Gene Systematic Name"]).sum())
	print('pav gene1 cu', pd.Series(pav.Gene1.unique()).isin(cu["Gene Systematic Name"]).sum(),\
		'pav gene2 cu', pd.Series(pav.Gene2.unique()).isin(cu["Gene Systematic Name"]).sum())
	print('pav gene1 sma', pd.Series(pav.Gene1.unique()).isin(sma["Gene Systematic Name"]).sum(),\
		'pav gene2 sma', pd.Series(pav.Gene2.unique()).isin(sma["Gene Systematic Name"]).sum())
	#
	cnv["Gene1_benomyl"] = cnv.apply(lambda x: 1 if x.Gene1 in beno["Gene Systematic Name"] else 0, axis=1) 
	cnv["Gene2_benomyl"] = cnv.apply(lambda x: 1 if x.Gene2 in beno["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene1_caffeine"] = cnv.apply(lambda x: 1 if x.Gene1 in caf["Gene Systematic Name"] else 0, axis=1) 
	cnv["Gene2_caffeine"] = cnv.apply(lambda x: 1 if x.Gene2 in caf["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene1_cuso4"] = cnv.apply(lambda x: 1 if x.Gene1 in cu["Gene Systematic Name"] else 0, axis=1) 
	cnv["Gene2_cuso4"] = cnv.apply(lambda x: 1 if x.Gene2 in cu["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene1_sma"] = cnv.apply(lambda x: 1 if x.Gene1 in sma["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene2_sma"] = cnv.apply(lambda x: 1 if x.Gene2 in sma["Gene Systematic Name"] else 0, axis=1)
	#
	try:
		snp.to_csv(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt", index=False, sep="\t")
	except:
		pass
	pav.to_csv(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt", index=False, sep="\t")
	cnv.to_csv(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt", index=False, sep="\t")
	#
	## number of shap interactions for each literature gene
	try:
		snp_beno1 = snp.where(snp.Gene1.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_beno2 = snp.where(snp.Gene2.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_beno1, snp_beno2], ignore_index=False, axis=1)
		snp_caf1 = snp.where(snp.Gene1.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_caf2 = snp.where(snp.Gene2.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_int_count, snp_caf1, snp_caf2], ignore_index=False, axis=1)
		snp_cu1 = snp.where(snp.Gene1.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_cu2 = snp.where(snp.Gene2.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_int_count, snp_cu1, snp_cu2], ignore_index=False, axis=1)
		snp_sma1 = snp.where(snp.Gene1.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_sma2 = snp.where(snp.Gene2.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_int_count, snp_sma1, snp_sma2], ignore_index=False, axis=1)
		snp_int_count.columns = ["NumInt_Gene1_beno", "NumInt_Gene2_beno", \
			"NumInt_Gene1_caffeine", "NumInt_Gene2_caffeine", "NumInt_Gene1_cuso4", \
			"NumInt_Gene2_cuso4", "NumInt_Gene1_sma", "NumInt_Gene2_sma"]
		snp_int_count.fillna(0, inplace=True)
		snp_int_count.index.name = "Gene"
		snp_int_count.to_csv(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes_int_counts.txt", index=False, sep="\t")
	except:
		pass
	pav_beno1 = pav.where(pav.Gene1.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_beno2 = pav.where(pav.Gene2.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_beno1 = cnv.where(cnv.Gene1.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_beno2 = cnv.where(cnv.Gene2.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_beno1, pav_beno2], ignore_index=False, axis=1)
	cnv_int_count = pd.concat([cnv_beno1, cnv_beno2], ignore_index=False, axis=1)
	#
	pav_caf1 = pav.where(pav.Gene1.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_caf2 = pav.where(pav.Gene2.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_caf1 = cnv.where(cnv.Gene1.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_caf2 = cnv.where(cnv.Gene2.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_int_count, pav_caf1, pav_caf2], ignore_index=False, axis=1)
	cnv_int_count = pd.concat([cnv_int_count, cnv_caf1, cnv_caf2], ignore_index=False, axis=1)
	#
	pav_cu1 = pav.where(pav.Gene1.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_cu2 = pav.where(pav.Gene2.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_cu1 = cnv.where(cnv.Gene1.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_cu2 = cnv.where(cnv.Gene2.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_int_count, pav_cu1, pav_cu2], ignore_index=False, axis=1)
	cnv_int_count = pd.concat([cnv_int_count, cnv_cu1, cnv_cu2], ignore_index=False, axis=1)
	#
	pav_sma1 = pav.where(pav.Gene1.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_sma2 = pav.where(pav.Gene2.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_sma1 = cnv.where(cnv.Gene1.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_sma2 = cnv.where(cnv.Gene2.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_int_count, pav_sma1, pav_sma2], ignore_index=False, axis=1)
	pav_int_count.columns = snp_int_count.columns
	pav_int_count.fillna(0, inplace=True)
	pav_int_count.index.name = "Gene"
	cnv_int_count = pd.concat([cnv_int_count, cnv_sma1, cnv_sma2], ignore_index=False, axis=1)
	cnv_int_count.columns = snp_int_count.columns
	cnv_int_count.fillna(0, inplace=True)
	cnv_int_count.index.name = "Gene"
	#
	pav_int_count.to_csv(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes_int_counts.txt", index=False, sep="\t")
	cnv_int_count.to_csv(f"{d}/PAV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes_int_counts.txt", index=False, sep="\t")
	#
	try:
		del snp, pav, cnv
	except:
		del pav, cnv

## Compare SHAP-based interactions to known GIs from the literature
# BioGRID database genetic interactions
biogrid = dt.fread("Data/BioGRID/yeast_gi_biogrid.txt").to_pandas()
biogrid = biogrid.iloc[:,[5,6,7,8,11,13,14,17,36]]
biogrid.columns = ["Systematic Name Interactor A", "Systematic Name Interactor B",
				   "Standard Name Interactor A", "Standard name Interactor B",
				   "Evidence", "Author", "PMID", "Throughput", "Organism"]
biogrid = biogrid.loc[biogrid.Organism=="Saccharomyces cerevisiae (S288c)",:]
biogrid = biogrid.loc[biogrid.Evidence.isin(["Synthetic Growth Defect",\
			"Synthetic Lethality", "Synthetic Rescue", "Negative Genetic",\
			"Positive Genetic"]),:] # remove overexpression gene pairs
biogrid_gp = biogrid.apply(lambda x: set([x["Systematic Name Interactor A"], \
			x["Systematic Name Interactor B"]]), axis=1).values # gene pairs
# biogrid_gp = {tuple(value): index for index, value in enumerate(list(biogrid_gp))}
biogrid_gp = {frozenset(set) for set in biogrid_gp}
len(biogrid_gp) # 438546

# Costanzo 2021 benomyl genetic interactions
costanzo = pd.read_excel("Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx",\
	engine="openpyxl", sheet_name="Genome-scale_Benomyl")
# costanzo = dt.fread("Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx/Genome-scale_Benomyl") # xlrd no longer supports xlsx
# costanzo = costanzo.to_pandas()
costanzo = costanzo.iloc[:,0:10]
costanzo = costanzo.loc[(costanzo.mean_condition_epsilon.abs() > 0.12) & \
	(costanzo.condition_p_value < 0.05),:] # Stringent criteria filter
costanzo.insert(0, "array_gene", costanzo.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo.insert(0, "query_gene", costanzo.apply(lambda x: x.query_orf.split("_")[0], axis=1))
costanzo_gp = costanzo.apply(lambda x: set([x.query_gene, x.array_gene]), axis=1).values # gene pairs
costanzo_gp = {frozenset(set) for set in costanzo_gp}
len(costanzo_gp) # 3472
len(biogrid_gp.union(costanzo_gp)) # 440463

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction"
num_int = {}
for env in target_envs:
	print(env)
	try:
		snp = dt.fread(glob.glob(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt")[0]).to_pandas()
	except:
		pass # missing benomyl snp summed file
	pav = dt.fread(glob.glob(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt")[0]).to_pandas()
	cnv = dt.fread(glob.glob(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt")[0]).to_pandas()#
	#
	# Is there experimental evidence for these SHAP interactions?
	try:
		snp["Known_biogrid"] = snp.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
		snp["Known_costanzo_2021"] = snp.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
		print('snp biogrid', snp["Known_biogrid"].value_counts())
		print('snp costanzo', snp["Known_costanzo_2021"].value_counts())
		num_int[("biogrid", "snp", env)] = snp["Known_biogrid"].value_counts()
		num_int[("costanzo", "snp", env)] = snp["Known_costanzo_2021"].value_counts()
	except:
		pass
	pav["Known_biogrid"] = pav.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
	pav["Known_costanzo_2021"] = pav.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
	print('pav biogrid', pav["Known_biogrid"].value_counts())
	print('pav costanzo', pav["Known_costanzo_2021"].value_counts())
	num_int[("biogrid", "pav", env)] = pav["Known_biogrid"].value_counts()
	num_int[("costanzo", "pav", env)] = pav["Known_costanzo_2021"].value_counts()
	cnv["Known_biogrid"] = cnv.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
	cnv["Known_costanzo_2021"] = cnv.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
	print('cnv biogrid', cnv["Known_biogrid"].value_counts())
	print('cnv costanzo', cnv["Known_costanzo_2021"].value_counts())
	num_int[("biogrid", "cnv", env)] = cnv["Known_biogrid"].value_counts()
	num_int[("costanzo", "cnv", env)] = cnv["Known_costanzo_2021"].value_counts()
	#
	# make sure it's saved as 0s and 1s not a boolean
	try:
		snp[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
			"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
			snp[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
			"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
		snp.to_csv(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	except:
		pass
	pav[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
		pav[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
	cnv[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
		cnv[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
	pav.to_csv(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	cnv.to_csv(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	# try:
	# 	del snp, pav, cnv
	# except:
	# 	del pav, cnv

num_int = pd.DataFrame(num_int).T
num_int.fillna(0, inplace=True)

biogrid_sub = num_int.loc[num_int.index.get_level_values(0)=="biogrid",:]
sns.barplot(y = biogrid_sub[1].values,\
	x = biogrid_sub.index.get_level_values(2).values,\
	hue = biogrid_sub.index.get_level_values(1).values,
	palette="viridis")
plt.xticks(rotation=45, size=7)
plt.ylabel("Number of known GIs identified")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/Num_known_GIs_biogrid_identified_barplot.pdf")
plt.close()

costanzo_sub = num_int.loc[num_int.index.get_level_values(0)=="costanzo",:]
sns.barplot(y = costanzo_sub[1].values,\
	x = costanzo_sub.index.get_level_values(2).values,\
	hue = costanzo_sub.index.get_level_values(1).values,
	palette="viridis")
plt.xticks(rotation=45, size=7)
plt.ylabel("Number of known GIs identified")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/Num_known_GIs_costanzo_identified_barplot.pdf")
plt.close()