#!/usr/bin/env python3

import os
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt

################################################################################
### TABLE S1 + FIGURE S? (FS curves)
################################################################################
os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

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
path = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/baseline"
rf_pc = pd.read_csv(f"{path}/RESULTS_PC_reg.txt", sep="\t")
rf_pc.insert(0, "new_cond", rf_pc.replace({"Y": mapping})["Y"]) # add full condition names
rf_pc = rf_pc.sort_values(by="r2_test", ascending=False)
cond_order = rf_pc["new_cond"]
rf_pc.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t", index=False) #********* done

## SNP Models
rf_snp = pd.read_csv(f"{path}/RESULTS_SNP_reg.txt", sep="\t")
rf_snp.insert(0, "new_cond", rf_snp.replace({"Y": mapping})["Y"]) # add full condition names
rf_snp.set_index("new_cond", inplace=True)
rf_snp = rf_snp.loc[cond_order] # sort by PC order
rf_snp.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_baseline.txt", sep="\t") #********* done

## PAV Models
path = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
rf_orf = pd.read_csv(f"{path}/baseline/RESULTS_reg.txt", sep="\t")
rf_pav = rf_orf.loc[rf_orf.ID.str.contains("_pav_baseline")]
rf_pav.insert(0, "new_cond", rf_pav.replace({"Y":mapping})["Y"])
rf_pav.set_index("new_cond", inplace=True)
rf_pav = rf_pav.loc[cond_order]
rf_pav.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_baseline.txt", sep="\t") #********* done

## CNV Models
rf_cnv = rf_orf[rf_orf.ID.str.contains("_cnv_baseline")]
rf_cnv.insert(0, "new_cond", rf_cnv.replace({"Y":mapping})["Y"])
rf_cnv.set_index("new_cond", inplace=True)
rf_cnv = rf_cnv.loc[cond_order]
rf_cnv.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_baseline.txt", sep="\t") #********* done

## Merge them into one table
combined = pd.concat([rf_pc, rf_snp.reset_index(), rf_pav.reset_index(), rf_cnv.reset_index()])
combined.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_ALL_baseline.txt", sep="\t", index=None) #********* done

################### XGBOOST BASELINE PREDICTION PERFORMANCES ###################
## PCs and SNPs
d = "/mnt/gs21/scratch/seguraab/yeast_project/"
xgb_snp = dt.fread(f"{d}SNP_yeast_XGBoost_results/baseline/RESULTS_xgboost.txt").to_pandas() #********* done

# Keep only relevant columns and insert missing information
xgb_snp.rename(columns={"R2_val":"r2_val", "R2_test":"r2_test",
	"R2_val_sd":"r2_val_sd", "R2_test_sd":"r2_test_sd"}, inplace=True)
xgb_snp.drop(columns=['MSE_val', 'MSE_val_sd', 'RMSE_val', 'RMSE_val_sd',
	'EVS_val', 'EVS_val_sd', 'MSE_test', 'MSE_test_sd', 'RMSE_test', 'RMSE_test_sd',
	'EVS_test', 'EVS_test_sd'], inplace=True)
xgb_snp.insert(10, "r2_val_se", xgb_snp.r2_val_sd/np.sqrt(10))
xgb_snp.insert(13, "PCC_val_se", xgb_snp.PCC_val_sd/np.sqrt(10))
xgb_snp.insert(16, "r2_test_se", xgb_snp.r2_test_sd/np.sqrt(10))
xgb_snp.insert(19, "PCC_test_se", xgb_snp.PCC_test_sd/np.sqrt(10))
xgb_snp.insert(0, "Alg", "XGBoost")

## PAVs and CNVs
xgb_orf = dt.fread(f"{d}ORF_yeast_XGBoost_results/baseline/RESULTS_xgboost.txt").to_pandas() #********* done
xgb_orf.rename(columns={"R2_val":"r2_val", "R2_test":"r2_test",
	"R2_val_sd":"r2_val_sd", "R2_test_sd":"r2_test_sd"}, inplace=True)
xgb_orf.drop(columns=['MSE_val', 'MSE_val_sd', 'RMSE_val', 'RMSE_val_sd',
	'EVS_val', 'EVS_val_sd', 'MSE_test', 'MSE_test_sd', 'RMSE_test', 'RMSE_test_sd',
	'EVS_test', 'EVS_test_sd'], inplace=True)
xgb_orf.insert(10, "r2_val_se", xgb_orf.r2_val_sd/np.sqrt(10))
xgb_orf.insert(13, "PCC_val_se", xgb_orf.PCC_val_sd/np.sqrt(10))
xgb_orf.insert(16, "r2_test_se", xgb_orf.r2_test_sd/np.sqrt(10))
xgb_orf.insert(19, "PCC_test_se", xgb_orf.PCC_test_sd/np.sqrt(10))
xgb_orf.insert(0, "Alg", "XGBoost")

################### rrBLUP BASELINE PREDICTION PERFORMANCES ####################
## PCs and SNPs
rrblup_snp = pd.read_csv(f"{d}SNP_yeast_rrBLUP_results/baseline/RESULTS_rrblup.txt", sep="\t")
rrblup_snp.rename(columns={"NumFeatures":"FeatureNum", "CV-Fold":"CVfold",
	"NumRepetitions":"CV_rep", "R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_val_se":"r2_val_se", "R2_test":"r2_test", "R2_test_sd":"r2_test_sd",
	"R2_test_se":"r2_test_se"}, inplace=True)
rrblup_snp.NumInstances = 625
rrblup_snp.insert(0, "Alg", "rrBLUP")
rrblup_snp.insert(1, "Data", rrblup_snp.apply(lambda x: "SNP" if "snp" in x.ID else "PCs_tassel", axis=1)) #********* done

## PAVs and CNVs
rrblup_orf = pd.read_csv(f"{d}ORF_yeast_rrBLUP_results/baseline/RESULTS_rrblup.txt", sep="\t")
rrblup_orf.rename(columns={"NumFeatures":"FeatureNum", "CV-Fold":"CVfold",
	"NumRepetitions":"CV_rep", "R2_val":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_val_se":"r2_val_se", "R2_test":"r2_test", "R2_test_sd":"r2_test_sd",
	"R2_test_se":"r2_test_se"}, inplace=True)
rrblup_orf.NumInstances = 625
rrblup_orf.insert(0, "Alg", "rrBLUP")
rrblup_orf.insert(1, "Data", rrblup_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1)) #********* done

############### BAYESIAN LASSO BASELINE PREDICTION PERFORMANCES ################
## SNPs
bl_snp = pd.read_csv(f"{d}SNP_yeast_BL_results/baseline/RESULTS_BL.txt", sep="\t")
bl_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_snp.insert(3, "Data", bl_snp.apply(lambda x: "SNP" if "snp" in x.ID else "PCs_tassel", axis=1)) #********* done

## PAVs and CNVs
bl_orf = pd.read_csv(f"{d}ORF_yeast_BL_results/baseline/RESULTS_BL.txt", sep="\t")
bl_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_orf.insert(3, "Data", bl_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1)) #********* done

################### BAYESC BASELINE PREDICTION PERFORMANCES ####################
## SNPs
bayesc_snp = pd.read_csv(f"{d}SNP_yeast_BayesC_results/baseline/RESULTS_BayesC.txt", sep="\t")
bayesc_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_snp.insert(3, "Data", bayesc_snp.apply(lambda x: "SNP" if "snp" in x.ID else "PCs_tassel", axis=1)) #********* done

## PAVs and CNVs
bayesc_orf = pd.read_csv(f"{d}ORF_yeast_BayesC_results/baseline/RESULTS_BayesC.txt", sep="\t")
bayesc_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_orf.insert(3, "Data", bayesc_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1)) #********* done

############# MERGE ALL ALGORITHM BASELINE RESULTS INTO ONE TABLE ##############
combined.insert(0, "Data", combined.apply(lambda x: "PCs_tassel" if "PC" in x.ID else \
								 ("PAV" if "pav" in x.ID else \
								 ("CNV" if "cnv" in x.ID else "SNP")), axis=1))
combined.drop(columns=["DateTime", "RunTime", "ID", "MSE_val", "MSE_val_sd", "MSE_val_se", \
						 "MSE_test", "MSE_test_sd", "MSE_test_se", "Tag", "new_cond", \
						 "EVS_val", "EVS_val_sd", "EVS_val_se", "EVS_test", "EVS_test_sd", \
						 "EVS_test_se"], inplace=True)
combined.rename(columns={"Y":"Trait"}, inplace=True)

all_alg = pd.concat([combined, xgb_snp, rrblup_snp, bl_snp, bayesc_snp, \
					 xgb_orf, rrblup_orf, bl_orf, bayesc_orf])
all_alg.drop(columns=["Date", "RunTime", "ID"], inplace=True)
all_alg.insert(1, "new_cond", all_alg.replace({"Trait":mapping})["Trait"])
all_alg.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t", index=False) #********* done

######## RF PREDICTION PERFORMANCES USING RF FEATURE SELECTION FEATURES ########
## SNPs
d = "/mnt/gs21/scratch/seguraab/yeast_project"
rf_snp = pd.read_csv(f"{d}/SNP_yeast_RF_results/fs/RESULTS_reg.txt", sep="\t")
rf_snp = rf_snp.loc[(rf_snp.DateTime.str.contains("2024-04-1") |
	rf_snp.DateTime.str.contains("2024-04-2")),:]
rf_snp.insert(0, "new_cond", rf_snp.replace({"Y": mapping})["Y"]) # add full condition names
rf_fs = pd.DataFrame(columns=rf_snp.columns) # Optimal FS features
fig, ax = plt.subplots(nrows=7, ncols=5, sharex=True, sharey=True, figsize=(8.5,11)) # FS curves
for i,env in enumerate(mapping.keys()):
	tmp = rf_snp.loc[rf_snp.Y==env,:].sort_values(by="FeatureNum")
	tmp = tmp.loc[tmp.FeatureNum <= 30000,:]
	tmp.drop_duplicates(subset="FeatureNum", keep="last", inplace=True)
	# Check if any models are missing (should be 40 rows)
	print(env, tmp.shape)
	print(tmp.FeatureNum.values)
	# Now get the optimal number of features
	rf_fs = pd.concat([rf_fs, tmp.loc[tmp.r2_val==max(tmp.r2_val),:]])
	# Generate the feature selection curve
	row_idx, col_idx = divmod(i, 5) # subplot coordinates
	ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum, y=tmp.r2_test, yerr=tmp.r2_test_sd, ecolor="black", linewidth="0.7")
	ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum, y=tmp.r2_val, yerr=tmp.r2_val_sd, ecolor="black", linewidth="0.7")
	# plot verticle line to denote optimal number of features
	id = tmp.r2_val.idxmax()
	opt = tmp.loc[tmp.index==id,"FeatureNum"]
	ax[row_idx][col_idx].vlines(x=opt, ymin=0, ymax=1, label=str(opt.values[0]), colors="black", linewidth=0.5)
	ax[row_idx][col_idx].legend()
	ax[row_idx][col_idx].set_title(env, fontsize="small")

plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_2/FS_curves_snp.pdf")
plt.close()

rf_fs.set_index("new_cond", inplace=True)
rf_fs = rf_fs.loc[cond_order,] # set to PC env order
rf_fs.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t", index=True) #********* done

## PAVs and CNVs
rf_orf = pd.read_csv(f"{d}/ORF_yeast_RF_results/fs/RESULTS_reg.txt", sep="\t") # Random Forest results
rf_orf = rf_orf.loc[rf_orf.DateTime.str.contains("2024-04-1"),:]
rf_orf.insert(0, "new_cond", rf_orf.replace({"Y": mapping})["Y"]) # add full condition names
rf_pav = pd.DataFrame(columns=rf_orf.columns) # Optimal FS features
rf_cnv = pd.DataFrame(columns=rf_orf.columns)
for data in ["pav", "cnv"]:
	fig, ax = plt.subplots(nrows=7, ncols=5, sharex=True, sharey=True, figsize=(8.5,11))
	for i,env in enumerate(mapping.keys()):
		tmp = rf_orf.loc[(rf_orf.Y==env) & (rf_orf.ID.str.contains(data)),:]
		tmp.drop_duplicates(subset="FeatureNum", keep="last", inplace=True)
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
		ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum, y=tmp.r2_test, yerr=tmp.r2_test_sd, ecolor="black", linewidth="0.7")
		ax[row_idx][col_idx].errorbar(x=tmp.FeatureNum, y=tmp.r2_val, yerr=tmp.r2_val_sd, ecolor="black", linewidth="0.7")
		# plot verticle line to denote optimal number of features
		id = tmp.r2_val.idxmax()
		opt = tmp.loc[tmp.index==id,"FeatureNum"]
		ax[row_idx][col_idx].vlines(x=opt, ymin=0, ymax=1, label=str(opt.values[0]), colors="black", linewidth=0.5)
		ax[row_idx][col_idx].legend()
		ax[row_idx][col_idx].set_title(env, fontsize="small")
	try:
		plt.savefig(f"Scripts/Data_Vis/Section_2/FS_curves_{data}.pdf")
	except TypeError:
		continue
	plt.close()

rf_pav.set_index("new_cond", inplace=True)
rf_cnv.set_index("new_cond", inplace=True)
rf_pav = rf_pav.loc[cond_order,] # set to PC env order
rf_cnv = rf_cnv.loc[cond_order,]
rf_pav.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t", index=True)
rf_cnv.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t", index=True)

###################### XGBOOST FS PREDICTION PERFORMANCES ######################
## SNPs
xgb_snp = pd.read_csv(f"{d}/SNP_yeast_XGBoost_results/fs/RESULTS_xgboost.txt", sep="\t") #********* done
xgb_snp.rename(columns={"R2_val            ":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_test            ":"r2_test", "R2_test_sd":"r2_test_sd",
	"NumFeatures            ":"FeatureNum", "CV-Fold":"CVfold", "NumRepetitions":"CV_rep"},
	inplace=True)
xgb_snp.drop(columns=['MSE_val', 'MSE_val_sd            ', 'RMSE_val', 'RMSE_val_sd',
	'EVS_val', 'EVS_val_sd', 'MSE_test', 'MSE_test_sd            ', 'RMSE_test', 'RMSE_test_sd',
	'EVS_test', 'EVS_test_sd', "Date", "RunTime"], inplace=True)
xgb_snp.insert(8, "r2_val_se", xgb_snp.r2_val_sd/np.sqrt(10))
xgb_snp.insert(11, "PCC_val_se", xgb_snp.PCC_val_sd/np.sqrt(10))
xgb_snp.insert(14, "r2_test_se", xgb_snp.r2_test_sd/np.sqrt(10))
xgb_snp.insert(17, "PCC_test_se", xgb_snp.PCC_test_sd/np.sqrt(10))
xgb_snp.insert(0, "Alg", "XGBoost")

## PAVs and CNVs
xgb_orf = pd.read_csv(f"{d}/ORF_yeast_XGBoost_results/fs/RESULTS_xgboost.txt", sep="\t") #********* done
xgb_orf.rename(columns={"R2_val            ":"r2_val", "R2_val_sd":"r2_val_sd",
	"R2_test            ":"r2_test", "R2_test_sd":"r2_test_sd",
	"NumFeatures            ":"FeatureNum", "CV-Fold":"CVfold", "NumRepetitions":"CV_rep"},
	inplace=True)
xgb_orf.drop(columns=['MSE_val', 'MSE_val_sd            ', 'RMSE_val', 'RMSE_val_sd',
	'EVS_val', 'EVS_val_sd', 'MSE_test', 'MSE_test_sd            ', 'RMSE_test', 'RMSE_test_sd',
	'EVS_test', 'EVS_test_sd', "Date", "RunTime"], inplace=True)
xgb_orf.insert(8, "r2_val_se", xgb_orf.r2_val_sd/np.sqrt(10))
xgb_orf.insert(11, "PCC_val_se", xgb_orf.PCC_val_sd/np.sqrt(10))
xgb_orf.insert(14, "r2_test_se", xgb_orf.r2_test_sd/np.sqrt(10))
xgb_orf.insert(17, "PCC_test_se", xgb_orf.PCC_test_sd/np.sqrt(10))
xgb_orf.insert(0, "Alg", "XGBoost")

###################### rrBLUP FS PREDICTION PERFORMANCES #######################
## SNPs
rrblup_snp = pd.read_csv(f"{d}/SNP_yeast_rrBLUP_results/fs_based_on_rf/RESULTS_rrblup.txt", sep="\t") #********* done
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
rrblup_orf.insert(1, "Data", rrblup_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1)) #********* done

################## BAYESIAN LASSO FS PREDICTION PERFORMANCES ###################
## SNPs
bl_snp = pd.read_csv(f"{d}/SNP_yeast_BL_results/fs/RESULTS_BL.txt", sep="\t") #********* done
bl_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_snp.insert(0, "Data", "SNP")

## PAVs and CNVs
bl_orf = pd.read_csv(f"{d}/ORF_yeast_BL_results/fs/RESULTS_BL.txt", sep="\t")
bl_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
					 "MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bl_orf.insert(0, "Data", bl_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1)) #********* done

###################### BAYESC FS PREDICTION PERFORMANCES #######################
## SNPs
bayesc_snp = pd.read_csv(f"{d}/SNP_yeast_BayesC_results/fs/RESULTS_BayesC.txt", sep="\t") #********* done
bayesc_snp.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_snp.insert(0, "Data", "SNP")

## PAVs and CNVs
bayesc_orf = pd.read_csv(f"{d}/ORF_yeast_BayesC_results/fs/RESULTS_BayesC.txt", sep="\t")
bayesc_orf.drop(columns=["Date", "RunTime", "MSE_val", "MSE_val_sd", "MSE_val_se", \
				"MSE_test", "MSE_test_sd", "MSE_test_se"], inplace=True)
bayesc_orf.insert(0, "Data", bayesc_orf.apply(lambda x: "PAV" if "pav" in x.ID else "CNV", axis=1)) #********* done

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
all_fs.insert(0, "new_cond", all_fs.replace({"Trait": mapping})["Trait"])
all_fs.drop(columns=["Date", "ID"], inplace=True)
all_fs.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_FS.txt", sep="\t", index=None) #********* done

all_baseline = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t")
all_single_env = pd.concat([all_baseline, all_fs], axis=0)
all_single_env["new_cond"] = all_single_env.replace({"Trait": mapping})["Trait"]
all_single_env.to_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_SINGLE_ENV.txt", sep="\t", index=None) #********* done

