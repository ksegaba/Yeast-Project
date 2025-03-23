#!/usr/bin/env python3
################################################################################
# TABLE S7: SPEARMAN's RHO FEATURE RANK PERCENTILE CORRELATIONS BETWEEN ENVS
################################################################################
import os
import pandas as pd
import datatable as dt
import numpy as np
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

## SPEARMAN's RHO CORRELATIONS BETWEEN ENVIRONMENTS
for i,mod_type in enumerate(["baseline", "FS"]):
	for j,data_type in enumerate(["snp", "pav", "cnv"]):
		for k,imp_type in enumerate(["imp", "shap"]):
			print(f"Calculating spearman's rho for Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_{data_type}_rank_per.tsv")
			df = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_{data_type}_rank_per.tsv").to_pandas()
			df.set_index(df.columns[0], inplace=True)
			# calculate spearman's rho
			rho = df.corr(method=lambda x, y: spearmanr(x, y).statistic)
			pval = df.corr(method=lambda x, y: spearmanr(x, y).pvalue)
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
			if (i==0 and j==0 and k==0):
				res = out.copy(deep=True)
				res.columns = [f"{imp_type}_{data_type}_{mod_type}_rho", \
							   f"{imp_type}_{data_type}_{mod_type}_pvalue"]
			else:
				out.columns = [f"{imp_type}_{data_type}_{mod_type}_rho", \
							   f"{imp_type}_{data_type}_{mod_type}_pvalue"]
				res = pd.concat([res, out], ignore_index=False, axis=1)
			del df, rho, pval, out

res.to_csv("Scripts/Data_Vis/Section_4/Table_S7_rank_per_corr_btwn_envs.tsv", sep="\t")
