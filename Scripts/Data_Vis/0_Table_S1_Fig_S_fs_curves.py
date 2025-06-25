#!/usr/bin/env python3

import os
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt

################################################################################
### TABLE S1 + FIGURE S? (FS curves)
################################################################################
os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Isolate growth condition labels
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

### RANDOM FOREST (RF) PREDICTION PERFORMANCES (BASELINE USING ALL FEATURES) ###
## PC Models
# path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/PC_yeast_RF_results"
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/baseline"
rf_pc = pd.read_csv(f"{path}/RESULTS_reg.txt", sep="\t")
rf_pc = rf_pc.loc[rf_pc.ID.str.contains("PCs_sklearn")]
rf_pc.insert(0, "new_cond", rf_pc.replace({"Y": mapping})["Y"]) # add full condition names
rf_pc = rf_pc.sort_values(by="r2_test", ascending=False)
cond_order = rf_pc["new_cond"]
rf_pc.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t", index=False)

## SNP Models
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/baseline"
rf_snp = pd.read_csv(f"{path}/RESULTS_reg.txt", sep="\t")
rf_snp = rf_snp.loc[~rf_snp.ID.str.contains("PCs_sklearn")]
rf_snp.insert(0, "new_cond", rf_snp.replace({"Y": mapping})["Y"]) # add full condition names
rf_snp.set_index("new_cond", inplace=True)
rf_snp = rf_snp.loc[cond_order] # sort by PC order
rf_snp.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_baseline.txt", sep="\t")

## PAV Models
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results"
rf_orf = pd.read_csv(f"{path}/baseline/RESULTS_reg.txt", sep="\t")
rf_pav = rf_orf.loc[rf_orf.ID.str.contains("_pav_baseline")]
rf_pav.insert(0, "new_cond", rf_pav.replace({"Y":mapping})["Y"])
rf_pav.set_index("new_cond", inplace=True)
rf_pav = rf_pav.loc[cond_order]
rf_pav.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_baseline.txt", sep="\t")

## CNV Models
rf_cnv = rf_orf[rf_orf.ID.str.contains("_cnv_baseline")]
rf_cnv.insert(0, "new_cond", rf_cnv.replace({"Y":mapping})["Y"])
rf_cnv.set_index("new_cond", inplace=True)
rf_cnv = rf_cnv.loc[cond_order]
rf_cnv.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_baseline.txt", sep="\t")

## Merge them into one table
combined = pd.concat([rf_pc, rf_snp.reset_index(), rf_pav.reset_index(), rf_cnv.reset_index()])
combined.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_ALL_baseline.txt", sep="\t", index=None)

################### XGBOOST BASELINE PREDICTION PERFORMANCES ###################
## PCs and SNPs
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/"
xgb_snp = pd.read_csv(f"{d}SNP_yeast_XGBoost_results/baseline/RESULTS_xgboost.tsv", sep="\t")
# xgb_pc = pd.read_csv(f"{d}PC_yeast_XGBoost_results/RESULTS_xgboost.tsv", sep="\t")
xgb_pc = xgb_snp.loc[xgb_snp.Tag=="PCs_sklearn"]
xgb_pc.insert(0, "Alg", "XGBoost")
xgb_pc.insert(1, "Data", "PCs_sklearn")
xgb_snp = xgb_snp.loc[xgb_snp.Tag=="SNP"]
xgb_snp.insert(0, "Alg", "XGBoost")
xgb_snp.insert(1, "Data", "SNP")
xgb_snp = pd.concat([xgb_pc, xgb_snp], axis=0, ignore_index=True)

# Keep only relevant columns and insert missing information
xgb_snp.rename(columns={"R2_val":"r2_val", "R2_test":"r2_test",
	"R2_val_sd":"r2_val_sd", "R2_test_sd":"r2_test_sd"}, inplace=True)
xgb_snp.drop(columns=["MSE_val", "MSE_val_sd", "RMSE_val", "RMSE_val_sd",
	"EVS_val", "EVS_val_sd", "MSE_test", "MSE_test_sd", "RMSE_test", "RMSE_test_sd",
	"EVS_test", "EVS_test_sd"], inplace=True)
xgb_snp.insert(12, "r2_val_se", xgb_snp.r2_val_sd/np.sqrt(10))
xgb_snp.insert(15, "PCC_val_se", xgb_snp.PCC_val_sd/np.sqrt(10))
xgb_snp.insert(18, "r2_test_se", xgb_snp.r2_test_sd/np.sqrt(10))
xgb_snp.insert(21, "PCC_test_se", xgb_snp.PCC_test_sd/np.sqrt(10))
xgb_snp.rename(columns={"NumFeatures":"FeatureNum", "CV_fold":"CVfold", "Y":"Trait"}, inplace=True)

## PAVs and CNVs
xgb_orf = pd.read_csv(f"{d}ORF_yeast_XGBoost_results/baseline/RESULTS_xgboost.tsv", sep="\t")
xgb_orf.rename(columns={"R2_val":"r2_val", "R2_test":"r2_test",
	"R2_val_sd":"r2_val_sd", "R2_test_sd":"r2_test_sd", "Tag":"Data"}, inplace=True)
xgb_orf.drop(columns=["MSE_val", "MSE_val_sd", "RMSE_val", "RMSE_val_sd",
	"EVS_val", "EVS_val_sd", "MSE_test", "MSE_test_sd", "RMSE_test", "RMSE_test_sd",
	"EVS_test", "EVS_test_sd"], inplace=True)
xgb_orf.insert(10, "r2_val_se", xgb_orf.r2_val_sd/np.sqrt(10))
xgb_orf.insert(13, "PCC_val_se", xgb_orf.PCC_val_sd/np.sqrt(10))
xgb_orf.insert(16, "r2_test_se", xgb_orf.r2_test_sd/np.sqrt(10))
xgb_orf.insert(19, "PCC_test_se", xgb_orf.PCC_test_sd/np.sqrt(10))
xgb_orf.insert(0, "Alg", "XGBoost")
xgb_orf.rename(columns={"NumFeatures":"FeatureNum", "CV_fold":"CVfold", "Y":"Trait"}, inplace=True)

################### rrBLUP BASELINE PREDICTION PERFORMANCES ####################
## PCs
# from sklearn.metrics import r2_score
# pheno = pd.read_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", index_col=0)

# rrblup_pc = []
# for env in mapping.keys():
#     test_preds = pd.read_csv(
#         f"{d}PC_yeast_rrBLUP_results/Predict_value_test_PCs_tassel_rrblup_{env}.csv",
#         index_col=0)
#     cv_preds = pd.read_csv(
#         f"{d}PC_yeast_rrBLUP_results/Predict_value_cv_PCs_tassel_rrblup_{env}.csv",
#         index_col=0)
#     test_r2 = test_preds.apply(lambda x: r2_score(pheno.loc[x.index, env], x))
#     cv_r2 = cv_preds.apply(lambda x: r2_score(pheno.loc[x.index, env], x))
#     test_pcc = test_preds.apply(lambda x: np.corrcoef(pheno.loc[x.index, env], x)[0,1])
#     cv_pcc = cv_preds.apply(lambda x: np.corrcoef(pheno.loc[x.index, env], x)[0,1])
#     rrblup_pc.append(["rrBLUP", "PCs_tassel", env, f"PCs_tassel_{env}_rrblup",
#         625, 5, 5, 20, cv_r2.mean(), cv_r2.std(), cv_r2.std()/np.sqrt(20),
#         cv_pcc.mean(), cv_pcc.std(), cv_pcc.std()/np.sqrt(20), test_r2.mean(),
#         test_r2.std(), test_r2.std()/np.sqrt(20), test_pcc.mean(), test_pcc.std(),
#         test_pcc.std()/np.sqrt(20)])

# rrblup_pc = pd.DataFrame(rrblup_pc, columns=["Alg", "Data", "Trait", "ID",
#     "NumInstances", "FeatureNum", "CVfold", "CV_rep", "r2_val", "r2_val_sd",
#     "r2_val_se", "PCC_val", "PCC_val_sd", "PCC_val_se", "r2_test", "r2_test_sd",
#     "r2_test_se", "PCC_test", "PCC_test_sd", "PCC_test_se"]) # make results file

## SNPs & PCs
rrblup_snp = pd.read_csv(f"{d}SNP_yeast_rrBLUP_results/baseline/RESULTS_rrblup.txt", sep="\t")
rrblup_snp.rename(columns={"NumFeatures":"FeatureNum", "CV-Fold":"CVfold",
	"NumRepetitions":"CV_rep", "R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_val_se":"r2_val_se", "R2_test":"r2_test", "R2_test_sd":"r2_test_sd",
	"R2_test_se":"r2_test_se"}, inplace=True)
rrblup_snp.NumInstances = 625
rrblup_snp.insert(0, "Alg", "rrBLUP")
rrblup_snp.insert(1, "Data", rrblup_snp.apply(lambda x: "SNP" if "snp" in x.ID else "PCs_sklearn", axis=1))
# rrblup_snp = pd.concat([rrblup_pc, rrblup_snp.drop(columns="Date")], axis=0)

## PAVs and CNVs
rrblup_orf = pd.read_csv(f"{d}ORF_yeast_rrBLUP_results/baseline/RESULTS_rrblup.txt", sep="\t")
rrblup_orf.rename(columns={"NumFeatures":"FeatureNum", "CV-Fold":"CVfold",
	"NumRepetitions":"CV_rep", "R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_val_se":"r2_val_se", "R2_test":"r2_test", "R2_test_sd":"r2_test_sd",
	"R2_test_se":"r2_test_se"}, inplace=True)
rrblup_orf.NumInstances = 625
rrblup_orf.insert(0, "Alg", "rrBLUP")
rrblup_orf.insert(1, "Data", rrblup_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))

############### BAYESIAN LASSO BASELINE PREDICTION PERFORMANCES ################
## PCs and SNPs
# bl_pc = pd.read_csv(f"{d}PC_yeast_BL_results/RESULTS_BL.txt", sep="\t")
bl_snp = pd.read_csv(f"{d}SNP_yeast_BL_results/baseline/RESULTS_BL.txt", sep="\t")
# bl_snp = pd.concat([bl_pc, bl_snp], axis=0)
bl_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_snp.insert(3, "Data", bl_snp.apply(lambda x: "SNP" if "snp" in x.ID else "PCs_sklearn", axis=1))

## PAVs and CNVs
bl_orf = pd.read_csv(f"{d}ORF_yeast_BL_results/baseline/RESULTS_BL.txt", sep="\t")
bl_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_orf.insert(3, "Data", bl_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))

################### BAYESC BASELINE PREDICTION PERFORMANCES ####################
## PCs and SNPs
# bayesc_pc = pd.read_csv(f"{d}PC_yeast_BayesC_results/RESULTS_BayesC.txt", sep="\t")
bayesc_snp = pd.read_csv(f"{d}SNP_yeast_BayesC_results/baseline/RESULTS_BayesC.txt", sep="\t")
# bayesc_snp = pd.concat([bayesc_pc, bayesc_snp], axis=0)
bayesc_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_snp.insert(3, "Data", bayesc_snp.apply(lambda x: "SNP" if "snp" in x.ID else "PCs_sklearn", axis=1))

## PAVs and CNVs
bayesc_orf = pd.read_csv(f"{d}ORF_yeast_BayesC_results/baseline/RESULTS_BayesC.txt", sep="\t")
bayesc_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_orf.insert(3, "Data", bayesc_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))

############# MERGE ALL ALGORITHM BASELINE RESULTS INTO ONE TABLE ##############
combined = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_ALL_baseline.txt", sep="\t")
combined.insert(0, "Data", combined.apply(lambda x: "PCs_sklearn" if "PCs" in x.ID else \
								 ("PAV" if "pav" in x.ID else \
								 ("CNV" if "cnv" in x.ID else "SNP")), axis=1))
combined.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se", "Tag", "new_cond", \
						 "EVS_val", "EVS_val_sd", "EVS_val_se", "EVS_test", "EVS_test_sd", \
						 "EVS_test_se"], inplace=True)
combined.rename(columns={"Y":"Trait"}, inplace=True)

all_alg = pd.concat([combined, xgb_snp, rrblup_snp, bl_snp, bayesc_snp, \
					 xgb_orf, rrblup_orf, bl_orf, bayesc_orf])
all_alg.drop(columns=["Date", "RunTime", "ID", "Tag"], inplace=True)
all_alg.insert(1, "new_cond", all_alg.replace({"Trait":mapping})["Trait"])

# Some rows have missing metrics that didn't get printed in their respective results files
all_alg.isna().sum()
pd.set_option("display.max_columns", None)
all_alg.reset_index(inplace=True)
all_alg.drop(columns="index", inplace=True)

# Fix the missing values
from sklearn.metrics import r2_score
pheno = pd.read_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", index_col=0)
test = pd.read_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/Test.txt", header=None)

all_alg.loc[all_alg.r2_val.isna(),:] # 610  PAV  YPD CuSO4 10 mM  YPDCUSO410MM  Bayesian LASSO

bl_cuso4_baseline_cv = pd.read_csv(
    f"{d}ORF_yeast_BL_results/baseline/Predict_value_cv_baseline_pav_BL_YPDCUSO410MM.csv") # some cv reps have NaNs
bl_cuso4_baseline_test = pd.read_csv(
    f"{d}ORF_yeast_BL_results/baseline/Predict_value_test_baseline_pav_BL_YPDCUSO410MM.csv") # some cv reps have NaNs
bl_cuso4_baseline_cv = bl_cuso4_baseline_cv.dropna(axis=1)
bl_cuso4_baseline_test = bl_cuso4_baseline_test.dropna(axis=1)
bl_cuso4_baseline_cv.index = pheno.loc[~pheno.index.isin(test[0])].index
bl_cuso4_baseline_test.index = pheno.loc[pheno.index.isin(test[0])].index
cv_r2 = bl_cuso4_baseline_cv.apply(lambda x: r2_score(pheno.loc[x.index, "YPDCUSO410MM"], x))
test_r2 = bl_cuso4_baseline_test.apply(lambda x: r2_score(pheno.loc[x.index, "YPDCUSO410MM"], x))
cv_pcc = bl_cuso4_baseline_cv.apply(lambda x: np.corrcoef(pheno.loc[x.index, "YPDCUSO410MM"], x)[0,1])
test_pcc = bl_cuso4_baseline_test.apply(lambda x: np.corrcoef(pheno.loc[x.index, "YPDCUSO410MM"], x)[0,1])

all_alg.iloc[610,:] = ["PAV", "YPD CuSO4 10 mM", "YPDCUSO410MM", "Bayesian LASSO",
    625, 7708, 5, 20, cv_r2.mean(), cv_r2.std(), cv_r2.std()/np.sqrt(20),
    cv_pcc.mean(), cv_pcc.std(), cv_pcc.std()/np.sqrt(20), test_r2.mean(),
    test_r2.std(), test_r2.std()/np.sqrt(20), test_pcc.mean(), test_pcc.std(),
    test_pcc.std()/np.sqrt(20)]

all_alg.isna().sum() # 7 rows
missing = all_alg.loc[all_alg.PCC_test.isna(),:]
            # Data                  new_cond          Trait      Alg  \
# 144  PCs_sklearn                  YPD 14ºC          YPD14  XGBoost   
# 165  PCs_sklearn     YPD Nystatin 10 µg/ml    YPDNYSTATIN  XGBoost   
# 174  PCs_sklearn          YPD Formamide 5%  YPDFORMAMIDE5  XGBoost   
# 197          SNP            YPD NaCl 1.5 M     YPDNACL15M  XGBoost   
# 424          PAV  YPD Methylviologen 20 mM          YPDMV  XGBoost   
# 459          PAV            YPD NaCl 1.5 M     YPDNACL15M  XGBoost   
# 485          CNV          YPD Formamide 5%  YPDFORMAMIDE5  XGBoost   
#      NumInstances  FeatureNum  CVfold  CV_rep    r2_val  r2_val_sd  r2_val_se  \
# 144           625           5       5      20 -0.000108   0.000084   0.000026   
# 165           625           5       5      20 -0.000974   0.000287   0.000091   
# 174           625           5       5      20 -0.000843   0.000307   0.000097   
# 197           625      118382       5      20 -0.000181   0.000103   0.000033   
# 424           625        7708       5      20 -0.000407   0.000251   0.000079   
# 459           625        7708       5      20 -0.000360   0.000274   0.000087   
# 485           625        7708       5      20 -0.000659   0.000287   0.000091 

for i in missing.index:
    row = all_alg.iloc[i,:]
    if "PCs" in row.Data:
        cv_preds = pd.read_csv(
            # f"{d}PC_yeast_XGBoost_results/{row.Trait}_PCs_tassel_baseline_cv_preds.csv", index_col=0)
            f"{d}SNP_yeast_XGBoost_results/baseline/{row.Trait}_PCs_sklearn_baseline_cv_preds.csv", index_col=0)
    elif "SNP" in row.Data:
        cv_preds = pd.read_csv(
            f"{d}SNP_yeast_XGBoost_results/baseline/{row.Trait}_snp_baseline_cv_preds.csv", index_col=0)
    elif "PAV" in row.Data:
        cv_preds = pd.read_csv(
            f"{d}ORF_yeast_XGBoost_results/baseline/{row.Trait}_pav_baseline_cv_preds.csv", index_col=0)
    elif "CNV" in row.Data:
        cv_preds = pd.read_csv(
            f"{d}ORF_yeast_XGBoost_results/baseline/{row.Trait}_cnv_baseline_cv_preds.csv", index_col=0)
    test_preds = cv_preds.loc[test[0],[c for c in cv_preds.columns if "test" in c]]
    print(test_preds.head()) # each column has a value that is repeated in each row. the values just change from one column to another
    test_pcc = test_preds.apply(lambda x: np.corrcoef(pheno.loc[x.index, row.Trait], x)[0,1])
    all_alg.iloc[i,:].loc[["PCC_test", "PCC_test_sd", "PCC_test_se"]] = [
        test_pcc.mean(), test_pcc.std(), test_pcc.std()/np.sqrt(20)]

all_alg.isna().sum() # the 8 missing remained NaNs because the test predictions
# have the exact same predicted value for each instance. Thus the pcc returned NaN.
all_alg.fillna(value=0, inplace=True)
all_alg.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t", index=False)

######## RF PREDICTION PERFORMANCES USING RF FEATURE SELECTION FEATURES ########
## SNPs
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project"
rf_snp = pd.read_csv(f"{d}/SNP_yeast_RF_results/fs/RESULTS_reg.txt", sep="\t")
rf_snp.insert(0, "new_cond", rf_snp.replace({"Y": mapping})["Y"]) # add full condition names
rf_fs = pd.DataFrame(columns=rf_snp.columns) # Optimal FS features
fig, ax = plt.subplots(nrows=7, ncols=5, sharex=True, sharey=True, figsize=(8.5,11)) # FS curves
for i,env in enumerate(mapping.keys()):
	try:
		tmp = rf_snp.loc[rf_snp.Y==env,:].sort_values(by="FeatureNum")
		tmp = tmp.loc[tmp.FeatureNum <= 30000,:]
		# Check if any models are missing (should be 40 rows)
		print(env, tmp.shape)
		print(tmp.FeatureNum.values)
		# Now get the optimal number of features
		rf_fs = pd.concat([rf_fs, tmp.loc[tmp.r2_val==max(tmp.r2_val),:]])
		# Generate the feature selection curve
		row_idx, col_idx = divmod(i, 5) # subplot coordinates
		print(row_idx, col_idx)
		ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum, y=tmp.r2_test, yerr=tmp.r2_test_sd, ecolor="black", linewidth="0.5")
		ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum, y=tmp.r2_val, yerr=tmp.r2_val_sd, ecolor="black", linewidth="0.5")
		# plot verticle line to denote optimal number of features
		id = tmp.r2_val.idxmax()
		opt = tmp.loc[tmp.index==id,"FeatureNum"]
		ax[row_idx][col_idx].vlines(x=opt, ymin=0, ymax=1, label=str(opt.values[0]), colors="black", linewidth=0.5)
		ax[row_idx][col_idx].legend()
		ax[row_idx][col_idx].set_title(env, fontsize="small")
	except:
		print(env)
		continue

plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_2/FS_curves_snp_v2.pdf")
plt.close()

rf_fs.set_index("new_cond", inplace=True)
rf_fs.sort_values(by="r2_test", ascending=False, inplace=True)
rf_fs.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t", index=True)

## PAVs and CNVs
rf_orf = pd.read_csv(f"{d}/ORF_yeast_RF_results/fs/RESULTS_reg.txt", sep="\t") # Random Forest results
rf_orf.insert(0, "new_cond", rf_orf.replace({"Y": mapping})["Y"]) # add full condition names
rf_orf = rf_orf.iloc[rf_orf.ID.drop_duplicates(keep="last", ignore_index=False).index,:]
rf_pav = pd.DataFrame(columns=rf_orf.columns) # Optimal FS features
rf_cnv = pd.DataFrame(columns=rf_orf.columns)
for data in ["pav", "cnv"]:
	fig, ax = plt.subplots(nrows=7, ncols=5, sharex=True, sharey=True, figsize=(8.5,11))
	for i,env in enumerate(mapping.keys()):
		try:
			tmp = rf_orf.loc[(rf_orf.Y==env) & (rf_orf.ID.str.contains(data)),:]
			tmp.sort_values(by="FeatureNum", inplace=True)
			# Check if any models are missing (should be 40 rows)
			print(env, tmp.shape)
			print(tmp.FeatureNum.values)
			# Now get the optimal number of features
			if data =="pav":
				rf_pav = pd.concat([rf_pav, tmp.loc[tmp.r2_val==max(tmp.r2_val),:]])
			if data == "cnv":
				rf_cnv = pd.concat([rf_cnv, tmp.loc[tmp.r2_val==max(tmp.r2_val),:]])
			# Generate the feature selection curve
			row_idx, col_idx = divmod(i, 5) # subplot coordinates
			ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum.values, y=tmp.r2_test.values, yerr=tmp.r2_test_sd.values, ecolor="black", linewidth="0.5")
			ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum.values, y=tmp.r2_val.values, yerr=tmp.r2_val_sd.values, ecolor="black", linewidth="0.5")
			# plot verticle line to denote optimal number of features
			id = tmp.r2_val.idxmax()
			opt = tmp.loc[tmp.index==id,"FeatureNum"]
			ax[row_idx][col_idx].vlines(x=opt, ymin=0, ymax=1, label=str(opt.values[0]), colors="black", linewidth=0.5)
			ax[row_idx][col_idx].legend()
			ax[row_idx][col_idx].set_title(env, fontsize="small")
		except:
			print(env)
			continue
	try:
		plt.savefig(f"Scripts/Data_Vis/Section_2/FS_curves_{data}_v2.pdf")
	except TypeError:
		continue
	plt.close()
	del fig, ax

rf_pav.set_index("new_cond", inplace=True)
rf_cnv.set_index("new_cond", inplace=True)
rf_pav.sort_values(by="r2_test", ascending=False, inplace=True)
rf_cnv.sort_values(by="r2_test", ascending=False, inplace=True)
rf_pav.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t", index=True)
rf_cnv.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t", index=True)

###################### XGBOOST FS PREDICTION PERFORMANCES ######################
## SNPs
xgb_snp = pd.read_csv(f"{d}/SNP_yeast_XGBoost_results/fs/RESULTS_xgboost.tsv", sep="\t")
xgb_snp.drop(columns=["MSE_val", "MSE_val_sd", "RMSE_val", "RMSE_val_sd",
	"EVS_val", "EVS_val_sd", "MSE_test", "MSE_test_sd", "RMSE_test", "RMSE_test_sd",
	"EVS_test", "EVS_test_sd", "Date", "RunTime"], inplace=True)
xgb_snp.rename(columns={"R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_test":"r2_test", "R2_test_sd":"r2_test_sd", "Y":"Trait",
	"NumFeatures":"FeatureNum", "CV_fold":"CVfold",
	"NumRepetitions":"CV_rep"}, inplace=True)
xgb_snp.insert(8, "r2_val_se", xgb_snp.r2_val_sd/np.sqrt(10))
xgb_snp.insert(11, "PCC_val_se", xgb_snp.PCC_val_sd/np.sqrt(10))
xgb_snp.insert(14, "r2_test_se", xgb_snp.r2_test_sd/np.sqrt(10))
xgb_snp.insert(17, "PCC_test_se", xgb_snp.PCC_test_sd/np.sqrt(10))
xgb_snp.insert(0, "Alg", "XGBoost")
xgb_snp.insert(1, "Data", "SNP")

## PAVs and CNVs
xgb_orf = pd.read_csv(f"{d}/ORF_yeast_XGBoost_results/fs/RESULTS_xgboost.tsv", sep="\t")
xgb_orf.rename(columns={"R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_test":"r2_test", "R2_test_sd":"r2_test_sd", "Y":"Trait",
	"NumFeatures":"FeatureNum", "CV_fold":"CVfold", "NumRepetitions":"CV_rep"},
	inplace=True)
xgb_orf.drop(columns=["MSE_val", "MSE_val_sd", "RMSE_val", "RMSE_val_sd",
	"EVS_val", "EVS_val_sd", "MSE_test", "MSE_test_sd", "RMSE_test", "RMSE_test_sd",
	"EVS_test", "EVS_test_sd", "Date", "RunTime"], inplace=True)
xgb_orf.insert(8, "r2_val_se", xgb_orf.r2_val_sd/np.sqrt(10))
xgb_orf.insert(11, "PCC_val_se", xgb_orf.PCC_val_sd/np.sqrt(10))
xgb_orf.insert(14, "r2_test_se", xgb_orf.r2_test_sd/np.sqrt(10))
xgb_orf.insert(17, "PCC_test_se", xgb_orf.PCC_test_sd/np.sqrt(10))
xgb_orf.insert(0, "Alg", "XGBoost")
xgb_orf.insert(1, "Data", xgb_orf.apply(lambda x: "PAV" if "pav" in x.Tag else "CNV", axis=1))

###################### rrBLUP FS PREDICTION PERFORMANCES #######################
## SNPs
rrblup_snp = pd.read_csv(f"{d}/SNP_yeast_rrBLUP_results/fs/RESULTS_rrblup.txt", sep="\t")
rrblup_snp.rename(columns={"NumFeatures":"FeatureNum", "CV-Fold":"CVfold",
	"NumRepetitions":"CV_rep", "R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_val_se":"r2_val_se", "R2_test":"r2_test", "R2_test_sd":"r2_test_sd",
	"R2_test_se":"r2_test_se"}, inplace=True)
rrblup_snp.NumInstances = 625
rrblup_snp.insert(0, "Alg", "rrBLUP")
rrblup_snp.insert(1, "Data", "SNP")

## PAVs and CNVs
rrblup_orf = pd.read_csv(f"{d}/ORF_yeast_rrBLUP_results/fs/RESULTS_rrblup.txt", sep="\t")
rrblup_orf.rename(columns={"NumFeatures":"FeatureNum", "CV-Fold":"CVfold",
	"NumRepetitions":"CV_rep", "R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_val_se":"r2_val_se", "R2_test":"r2_test", "R2_test_sd":"r2_test_sd",
	"R2_test_se":"r2_test_se"}, inplace=True)
rrblup_orf.NumInstances = 625
rrblup_orf.insert(0, "Alg", "rrBLUP")
rrblup_orf.insert(1, "Data", rrblup_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))

################## BAYESIAN LASSO FS PREDICTION PERFORMANCES ###################
## SNPs
bl_snp = pd.read_csv(f"{d}/SNP_yeast_BL_results/fs/RESULTS_BL.txt", sep="\t")
bl_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_snp.insert(0, "Data", "SNP")

## PAVs and CNVs
bl_orf = pd.read_csv(f"{d}/ORF_yeast_BL_results/fs/RESULTS_BL.txt", sep="\t")
bl_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_orf.insert(0, "Data", bl_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))

###################### BAYESC FS PREDICTION PERFORMANCES #######################
## SNPs
bayesc_snp = pd.read_csv(f"{d}/SNP_yeast_BayesC_results/fs/RESULTS_BayesC.txt", sep="\t")
bayesc_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_snp.insert(0, "Data", "SNP")

## PAVs and CNVs
bayesc_orf = pd.read_csv(f"{d}/ORF_yeast_BayesC_results/fs/RESULTS_BayesC.txt", sep="\t")
bayesc_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_orf.insert(0, "Data", bayesc_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1))

############# MERGE ALL ALGORITHM BASELINE RESULTS INTO ONE TABLE ##############
rf_fs.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se", "EVS_val", "EVS_val_sd",
				"EVS_val_se", "EVS_test", "EVS_test_sd", "EVS_test_se", "Tag"], inplace=True)
rf_fs.insert(0, "Data", "SNP")
rf_fs.rename(columns={"Y":"Trait"}, inplace=True)
rf_pav.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se", "EVS_val", "EVS_val_sd",
				"EVS_val_se", "EVS_test", "EVS_test_sd", "EVS_test_se", "Tag"], inplace=True)
rf_pav.insert(0, "Data", "PAV")
rf_pav.rename(columns={"Y":"Trait"}, inplace=True)
rf_cnv.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se", "EVS_val", "EVS_val_sd",
				"EVS_val_se", "EVS_test", "EVS_test_sd", "EVS_test_se", "Tag"], inplace=True)
rf_cnv.insert(0, "Data", "CNV")
rf_cnv.rename(columns={"Y":"Trait"}, inplace=True)

all_fs = pd.concat([rf_fs, xgb_snp, rrblup_snp, bl_snp, bayesc_snp, rf_pav, rf_cnv,
			xgb_orf, rrblup_orf, bl_orf, bayesc_orf], axis=0)
all_fs.drop(columns=["Date", "Tag", "ID"], inplace=True)
all_fs.insert(0, "new_cond", all_fs.replace({"Trait": mapping})["Trait"])

# missing values due to xgboost identical predict output as with the baseline models
all_fs.loc[all_fs.PCC_test.isna(),:] # 6 rows
all_fs.fillna(value=0, inplace=True)
all_fs.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_FS.txt", sep="\t", index=None)

all_baseline = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t")
all_single_env = pd.concat([all_baseline, all_fs], axis=0)
# all_single_env.to_csv("Scripts/Data_Vis/Section_2/Table_S1_RESULTS_ALL_SINGLE_ENV.txt", sep="\t", index=None)
all_single_env.to_csv("Scripts/Data_Vis/Section_2/Table_S1_RESULTS_ALL_SINGLE_ENV_CORRECTED.txt", sep="\t", index=None)