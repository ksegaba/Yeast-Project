#!/usr/bin/env python3
"""
Calculate the performance of the GxE models from GxE_models_multialg.py

Arguments:
	-i, --input: Predicted trait values files from GxE models
	-d, --dir: Working directory
	-a, --alg: Algorithm used
	-c, --cv_type: Cross-validation type

Usage:
	python get_GxE_performance.py -i <input> -d <dir> -a <alg> -c <cv_type>

Example:
	input="across-env/cv1/GxE_cv1_5envs_orf_baseline_*_preds_test_*.txt"
	dir="/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/"
	python get_GxE_performance.py -i $input -d $dir -a across-env -c cv1

Output:
	- RESULTS_GxE_performance.csv: Performance metrics for each GxE model
"""


__author__ = "Kenia Segura Aba"


import os
import argparse, sys
import re
import glob
import datatable as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm
from sklearn.metrics import r2_score, mean_squared_error


os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GxE")


def parse_arguments():
	"""Argument parser"""
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-i", "--input", type=str, help="Predicted trait values files from GxE models",
		required=True)
	parser.add_argument(
		"-d", "--dir", type=str, help="Working directory",
		required="True")
	parser.add_argument(
		"-a", "--alg", type=str, help="Algorithm used",
		required="True")
	parser.add_argument(
		"-c", "--cv_type", type=str, help="Cross-validation type",
		default="cv1")
	return parser.parse_args()


def se(x):
	"""Calculate the standard error"""
	return np.std(x, ddof=1)/np.sqrt(len(x))


def calc_metrics(file, cv_type, alg):
	# print(file)
	df = dt.fread(file)
	df = df.to_pandas()
	# print(df.shape)
	df = df.set_index(["Rep", "Trait"]) # set multiindex with Rep and Trait
	ID = re.search("baseline_([A-Z0-9]+_){5}", file).group()
	df = df.loc[:,Y_test.index] # make sure column order is the same
	to_append = pd.DataFrame()
	for trait in re.search("([A-Z0-9]+_){5}", ID).group().split("_"):
		# print(trait)
		try:
			r2 = df.apply(lambda x: r2_score(Y_test[trait], x), axis=1) # coefficient of determination
			mse = df.apply(lambda x: mean_squared_error(Y_test[trait], x), axis=1) # mean squared error
			pcc = df.corrwith(Y_test[trait], method="pearson", axis=1) # pearson's correlation
		except KeyError:
			pass
		to_append = to_append.append(
			{"CV_type":cv_type, "Alg":alg, "Data":"ORFs", "ID":ID, "Env":trait,
			"r2_test":r2.loc[r2.index.get_level_values("Trait")==trait,].mean(),
			"r2_test_sd":r2.loc[r2.index.get_level_values("Trait")==trait,].std(),
			"r2_test_se":se(r2.loc[r2.index.get_level_values("Trait")==trait,]),
			"MSE_test":mse.loc[mse.index.get_level_values("Trait")==trait,].mean(),
			"MSE_test_sd":mse.loc[mse.index.get_level_values("Trait")==trait,].std(),
			"MSE_test_se":se(mse.loc[mse.index.get_level_values("Trait")==trait,]),
			"PCC_test":pcc.loc[pcc.index.get_level_values("Trait")==trait,].mean(),
			"PCC_test_sd":pcc.loc[pcc.index.get_level_values("Trait")==trait,].std(),
			"PCC_test_se":se(pcc.loc[pcc.index.get_level_values("Trait")==trait,])},
			ignore_index=True)
	return (to_append)


def calc_metrics2(file, cv_type):
	# print(file)
	df = dt.fread(file)
	df = df.to_pandas().T
	# print(df.shape)
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
	# Arguments help
	if len(sys.argv)==1:
		print("Input file and save name required")
		sys.exit()
	
	# Arguments
	args = parse_arguments()
	dir = Path(args.dir)
	input = dir / Path(args.input)
	alg = args.alg
	cv_type = args.cv_type

	# Actual trait values
	Y = pd.read_csv("../../Data/Peter_2018/pheno.csv", index_col=0)
	test = pd.read_csv("../../Data/Peter_2018/Test.txt", header=None)
	Y_test = Y.loc[test[0],:]
	
	# Get files with predicted values
	# dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/markerXenv/"
	# cv1_files = glob.glob(f"{dir}cv1/GxE_cv1_orf_baseline_*_preds_test_*.txt")
	# cv2_files = glob.glob(f"{dir}cv2/GxE_cv2_orf_baseline_*_preds_test_*.txt")
	input_files = glob.glob(str(input))

	# save results to this dataframe
	save_file = dir/"RESULTS_GxE.txt"
	res = pd.DataFrame(columns=["CV_type", "Alg", "Data", "ID", "Env",
								"r2_test", "r2_test_sd", "r2_test_se",
								"MSE_test", "MSE_test_sd", "MSE_test_se",
								"PCC_test", "PCC_test_sd", "PCC_test_se"])
	
	if not os.path.exists(dir/"RESULTS_GxE.txt"): # print header if file doesn't exist
		res.to_csv(save_file, sep="\t", index=False)

	# Gather CV1 results
	# for i in range(1110, len(cv1_files)):
	# 	to_append = calc_metrics2(cv1_files[i], "cv1")
	# 	res = res.append(to_append, ignore_index=True)

	# # Gather CV2 results
	# for i in range(len(cv2_files)):
	# 	to_append = calc_metrics2(cv2_files[i], "cv2")
	# 	res = res.append(to_append, ignore_index=True)

	# Gather trait performances
	for i in tqdm(range(len(input_files))):
		to_append = calc_metrics(input_files[i], cv_type, alg)
		res = res.append(to_append, ignore_index=True)
	
	res.to_csv(save_file, sep="\t", index=False, mode="a", header=False)

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