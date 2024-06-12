#!/usr/bin/env python

import os
import re
import glob
import datatable as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GxE")

def se(x):
	"""Calculate the standard error"""
	return np.std(x, ddof=1)/np.sqrt(len(x))


def summarize_metrics(test_file, val_file, data_type):
	df_test = pd.read_csv(test_file)
	df_val = pd.read_csv(val_file)
	ID_test = re.search("[A-Z0-9]+_[A-Z0-9]+_baseline", test_file).group()
	ID_val = re.search("[A-Z0-9]+_[A-Z0-9]+_baseline", val_file).group()
	mean_test = df_test.groupby("trait").mean()
	sd_test = df_test.groupby("trait").std()
	se_test = df_test.groupby("trait").apply(lambda x: se(x))
	mean_val = df_val.groupby("trait").mean()
	sd_val = df_val.groupby("trait").std()
	se_val = df_val.groupby("trait").apply(lambda x: se(x))
	env1 = df_test.trait.unique()[0]
	env2 = df_test.trait.unique()[1]
	to_append1 = {"Data":data_type, "ID":ID_test, "Env":env1,
				  "r2_val":mean_val.loc[mean_val.index==env1, "r2"].values[0],
				  "r2_val_sd":sd_val.loc[sd_val.index==env1, "r2"].values[0],
				  "r2_val_se":se_val.loc[se_val.index==env1, "r2"].values[0],
				  "MSE_val":mean_val.loc[mean_val.index==env1, "mse"].values[0],
				  "MSE_val_sd":sd_val.loc[sd_val.index==env1, "mse"].values[0],
				  "MSE_val_se":se_val.loc[se_val.index==env1, "mse"].values[0],
				  "PCC_val":mean_val.loc[mean_val.index==env1, "pcc"].values[0],
				  "PCC_val_sd":sd_val.loc[sd_val.index==env1, "pcc"].values[0],
				  "PCC_val_se":se_val.loc[se_val.index==env1, "pcc"].values[0],
				  "r2_test":mean_test.loc[mean_test.index==env1, "r2"].values[0],
				  "r2_test_sd":sd_test.loc[sd_test.index==env1, "r2"].values[0],
				  "r2_test_se":se_test.loc[se_test.index==env1, "r2"].values[0],
				  "MSE_test":mean_test.loc[mean_test.index==env1, "mse"].values[0],
				  "MSE_test_sd":sd_test.loc[sd_test.index==env1, "mse"].values[0],
				  "MSE_test_se":se_test.loc[se_test.index==env1, "mse"].values[0],
				  "PCC_test":mean_test.loc[mean_test.index==env1, "pcc"].values[0],
				  "PCC_test_sd":sd_test.loc[sd_test.index==env1, "pcc"].values[0],
				  "PCC_test_se":se_test.loc[se_test.index==env1, "pcc"].values[0]}
	to_append2 = {"Data":data_type, "ID":ID_test, "Env":env2,
				  "r2_val":mean_val.loc[mean_val.index==env2, "r2"].values[0],
				  "r2_val_sd":sd_val.loc[sd_val.index==env2, "r2"].values[0],
				  "r2_val_se":se_val.loc[se_val.index==env2, "r2"].values[0],
				  "MSE_val":mean_val.loc[mean_val.index==env2, "mse"].values[0],
				  "MSE_val_sd":sd_val.loc[sd_val.index==env2, "mse"].values[0],
				  "MSE_val_se":se_val.loc[se_val.index==env2, "mse"].values[0],
				  "PCC_val":mean_val.loc[mean_val.index==env2, "pcc"].values[0],
				  "PCC_val_sd":sd_val.loc[sd_val.index==env2, "pcc"].values[0],
				  "PCC_val_se":se_val.loc[se_val.index==env2, "pcc"].values[0],
				  "r2_test":mean_test.loc[mean_test.index==env2, "r2"].values[0],
				  "r2_test_sd":sd_test.loc[sd_test.index==env2, "r2"].values[0],
				  "r2_test_se":se_test.loc[se_test.index==env2, "r2"].values[0],
				  "MSE_test":mean_test.loc[mean_test.index==env2, "mse"].values[0],
				  "MSE_test_sd":sd_test.loc[sd_test.index==env2, "mse"].values[0],
				  "MSE_test_se":se_test.loc[se_test.index==env2, "mse"].values[0],
				  "PCC_test":mean_test.loc[mean_test.index==env2, "pcc"].values[0],
				  "PCC_test_sd":sd_test.loc[sd_test.index==env2, "pcc"].values[0],
				  "PCC_test_se":se_test.loc[se_test.index==env2, "pcc"].values[0]}
	return(to_append1, to_append2)

# Get files with performance metrics
dir = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_MORF_results/"
test_files = glob.glob(f"{dir}*_results_test.csv")
val_files = glob.glob(f"{dir}*_results_cv.csv")

res = pd.DataFrame(columns=["Data", "ID", "Env", "r2_val", "r2_val_sd",
							"r2_val_se", "MSE_val", "MSE_val_sd",
							"MSE_val_se", "PCC_val", "PCC_val_sd",
							"PCC_val_se", "r2_test", "r2_test_sd",
							"r2_test_se", "MSE_test", "MSE_test_sd",
							"MSE_test_se", "PCC_test", "PCC_test_sd",
							"PCC_test_se"]) # results

# Gather results
for i in range(len(test_files)):
	val_file = re.sub("_test.csv", "_cv.csv", test_files[i])
	print(f"TEST FILE: {test_files[i]}")
	print(f"VAL FILE: {val_file}")
	to_append1, to_append2 = summarize_metrics(test_files[i], val_file, "SNPs")
	res = res.append(to_append1, ignore_index=True)
	res = res.append(to_append2, ignore_index=True)

res.to_csv(f"{dir}RESULTS_MORF.txt", sep="\t", index=False)

# boxplots for each trait
res2 = res.copy(deep=True)
res2 = res2.sort_values(by="r2_test", ascending=False)
res2.boxplot(column="r2_test", by="Env", figsize=(11, 5))
plt.xticks(rotation=50, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.savefig("MORF_R2_test_ALL_envs_summary.pdf")
plt.close()
