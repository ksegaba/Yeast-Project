#!/usr/bin/env python3
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


def calc_metrics(file, cv_type):
	print(file)
	df = dt.fread(file)
	df = df.to_pandas().T
	print(df.shape)
	df.columns = df.iloc[0,:] # set Rep as column label
	df = df.iloc[1:,:] # remove the row used as column labels
	ID = re.search("baseline_[A-Z0-9]+_[A-Z0-9]+", file).group()
	trait = re.search("[A-Z0-9]+.txt", file).group().split(".")[0]
	df = df.loc[Y_test.index,] # make sure row order is the same
	r2 = df.apply(lambda x: r2_score(Y_test[trait], x), axis=0) # coefficient of determination
	mse = df.apply(lambda x: mean_squared_error(Y_test[trait], x), axis=0) # mean squared error
	pcc = df.corrwith(Y_test[trait], method="pearson") # pearson's correlation
	to_append = {"CV_type":cv_type, "Data":"ORFs", "ID":ID, "Env":trait,
				"r2_test":r2.mean(), "r2_test_sd":r2.std(),
				"r2_test_se":se(r2), "MSE_test":mse.mean(),
				"MSE_test_sd":mse.std(), "MSE_test_se":se(mse),
				"PCC_test":pcc.mean(), "PCC_test_sd":pcc.std(),
				"PCC_test_se":se(pcc)}
	return (to_append)


if __name__=="__main__":
	# Actual trait values
	Y = pd.read_csv("../../Data/Peter_2018/pheno.csv", index_col=0)
	test = pd.read_csv("../../Data/Peter_2018/Test.txt", header=None)
	Y_test = Y.loc[test[0],:]
	
	# Get files with predicted values
	dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/markerXenv/"
	cv1_files = glob.glob(f"{dir}cv1/GxE_cv1_orf_baseline_*_preds_test_*.txt")
	cv2_files = glob.glob(f"{dir}cv2/GxE_cv2_orf_baseline_*_preds_test_*.txt")

	res = pd.DataFrame(columns=["CV_type", "Data", "ID", "Env", "r2_test",
								"r2_test_sd", "r2_test_se", "MSE_test",
								"MSE_test_sd", "MSE_test_se", "PCC_test",
								"PCC_test_sd", "PCC_test_se"]) # results
	# Gather CV1 results
	for i in range(1110, len(cv1_files)):
		to_append = calc_metrics(cv1_files[i], "cv1")
		res = res.append(to_append, ignore_index=True)

	# Gather CV2 results
	for i in range(len(cv2_files)):
		to_append = calc_metrics(cv2_files[i], "cv2")
		res = res.append(to_append, ignore_index=True)

	res.to_csv(f"{dir}RESULTS_GxE.txt", sep="\t", index=False)

	# boxplots for each trait
	# res2 = res.copy(deep=True)
	# res2 = res2.sort_values(by="PCC_test", ascending=False)
	# res2["r2_test"] = res2["PCC_test"].values**2
	# res2.boxplot(column="r2_test", by="Env", figsize=(11, 5))
	# plt.xticks(rotation=50, ha="right", rotation_mode="anchor")
	# plt.tight_layout()
	# plt.savefig("CV1_R2_ALL_envs_summary.pdf")
	# plt.close()
	# res2 = res2[res2["r2_test"]=="Series([], dtype: float64)",:]
	# PCC = res2.groupby("Env").mean()
	# PCC.sd = res2.groupby("Env").std()
	# PCC = PCC.sort_values(by="PCC_test", ascending=False)
	# y = PCC["PCC_test"]
	# y.plot.bar(yerr=[PCC.sd["PCC_test"], PCC.sd["PCC_test"]], figsize=(11, 5))
	# plt.xticks(rotation=50, ha="right", rotation_mode="anchor")
	# plt.tight_layout()
	# 