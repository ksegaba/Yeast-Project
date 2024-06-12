#!/usr/bin/env python3

import os
import re
import multiprocessing
import pandas as pd
import datatable as dt
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from tqdm import tqdm
from scipy.stats import ks_2samp
from functools import partial

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

################################################################################
### TABLE S4
################################################################################
## SPEARMAN'S RHO CORRELATIONS BETWEEN IMPORTANCE MEASURES
res = [["Model Type", "Data Type", "Env", "rho", "pval"]]
for data_type in ["snp", "pav", "cnv"]:
    # Feature selection models
    gini = dt.fread(f"Scripts/Data_Vis/Section_4/RF_FS_imp_{data_type}.tsv").to_pandas()
    shap = dt.fread(f"Scripts/Data_Vis/Section_4/RF_FS_shap_{data_type}.tsv").to_pandas()
    #
    # Baseline models
    gini_base = dt.fread(f"Scripts/Data_Vis/Section_4/RF_baseline_imp_{data_type}.tsv").to_pandas()
    shap_base = dt.fread(f"Scripts/Data_Vis/Section_4/RF_baseline_shap_{data_type}.tsv").to_pandas()
    if data_type == "snp":
        gini.set_index("snp", inplace=True)
        shap.set_index("snp", inplace=True)
        gini_base.set_index("snp", inplace=True)
        shap_base.set_index("snp", inplace=True)
    else:
        gini.set_index("orf", inplace=True)
        shap.set_index("orf", inplace=True)
        gini_base.set_index("orf", inplace=True)
        shap_base.set_index("orf", inplace=True)
    #
    ## Convert data to rank percentiles
    gini_env_rank_out = pd.DataFrame() # FS + remaining baseline features
    shap_env_rank_out = pd.DataFrame()
    gini_rank_out = pd.DataFrame() # just FS features
    shap_rank_out = pd.DataFrame()
    gini_base_rank_out = pd.DataFrame() # just baseline features
    shap_base_rank_out = pd.DataFrame()
    for env in mapping.keys():
        #### First rank the feature selection features
        # Concatenate the remaining baseline model features to the end of 
        # the feature selection models (for proper ranking of FS features)
        gini_env = gini.loc[:,env]
        shap_env = shap.loc[:,env]
        gini_env.dropna(inplace=True) # remove the extra features from the feature selection dataset
        shap_env.dropna(inplace=True)
        fs_gini_feat = gini_env[gini_env != 0].index
        fs_shap_feat = shap_env[shap_env != 0].index
        gini_env = pd.concat([gini_env, 
            gini_base.loc[~gini_base.index.isin(gini_env.index),env]], axis=0) # add the remaining features from baseline models
        shap_env = pd.concat([shap_env,
            shap_base.loc[~shap_base.index.isin(shap_env.index),env]], axis=0)
        gini_env = gini_env[gini_env != 0] # Drop features with zero importance
        shap_env = shap_env[shap_env != 0]
        #
        # Calculate the rank percentiles
        gini_env.sort_values(ascending=False, inplace=True)
        shap_env.sort_values(ascending=False, inplace=True)
        gini_env_rank = gini_env.rank(axis=0, method="average", numeric_only=True, pct=True) # ranks as percentiles (1= most important)
        shap_env_rank = shap_env.abs().rank(axis=0, method="average", numeric_only=True, pct=True) # based on absolute shap value
        gini_env_rank_out = pd.concat([gini_env_rank_out, gini_env_rank], axis=1, ignore_index=False)
        shap_env_rank_out = pd.concat([shap_env_rank_out, shap_env_rank], axis=1, ignore_index=False)
        # drop the baseline features (keep only the feature selection features)
        gini_env_rank = gini_env_rank[fs_gini_feat]
        shap_env_rank = shap_env_rank[fs_shap_feat]
        gini_rank_out = pd.concat([gini_rank_out, gini_env_rank], axis=1, ignore_index=False)
        shap_rank_out = pd.concat([shap_rank_out, shap_env_rank], axis=1, ignore_index=False)
        #### Now rank the baseline features
        gini_base_env = gini_base.loc[:,env].dropna()
        shap_base_env = shap_base.loc[:,env].dropna()
        gini_base_env = gini_base_env[gini_base_env != 0]
        shap_base_env = shap_base_env[shap_base_env != 0]
        gini_base_env.sort_values(ascending=False, inplace=True)
        shap_base_env.sort_values(ascending=False, inplace=True)
        gini_base_env_rank = gini_base_env.rank(axis=0, method="average", numeric_only=True, pct=True)
        shap_base_env_rank = shap_base_env.abs().rank(axis=0, method="average", numeric_only=True, pct=True)
        gini_base_rank_out = pd.concat([gini_base_rank_out, gini_base_env_rank], axis=1, ignore_index=False)
        shap_base_rank_out = pd.concat([shap_base_rank_out, shap_base_env_rank], axis=1, ignore_index=False)
    gini_env_rank_out.to_csv(f"Scripts/Data_Vis/Section_4/RF_FS_plus_baseline_imp_{data_type}_rank_per.tsv", sep="\t")
    shap_env_rank_out.to_csv(f"Scripts/Data_Vis/Section_4/RF_FS_plus_baseline_shap_{data_type}_rank_per.tsv", sep="\t")
    gini_rank_out.to_csv(f"Scripts/Data_Vis/Section_4/RF_FS_imp_{data_type}_rank_per.tsv", sep="\t")
    shap_rank_out.to_csv(f"Scripts/Data_Vis/Section_4/RF_FS_shap_{data_type}_rank_per.tsv", sep="\t")
    gini_base_rank_out.to_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_imp_{data_type}_rank_per.tsv", sep="\t")
    shap_base_rank_out.to_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_shap_{data_type}_rank_per.tsv", sep="\t")
    #
    ## Calculate spearman's rho
    for env in mapping.keys():
        # First, for the baseline features
        df = pd.concat([gini_base_rank_out.loc[:,env], shap_base_rank_out.loc[:,env]], ignore_index=False, axis=1).dropna()
        rho = df.corr(method=lambda x, y: spearmanr(x, y, alternative="two-sided").statistic)
        pval = df.corr(method=lambda x, y: spearmanr(x, y, alternative="two-sided").pvalue)
        res.append(["baseline", data_type, env, rho.iloc[0,1], pval.iloc[0,1]])
        # Second, for the feature selection features
        df = pd.concat([gini_rank_out.loc[:,env], shap_rank_out.loc[:,env]], ignore_index=False, axis=1).dropna()
        rho = df.corr(method=lambda x, y: spearmanr(x, y, alternative="two-sided").statistic)
        pval = df.corr(method=lambda x, y: spearmanr(x, y, alternative="two-sided").pvalue)
        res.append(["FS", data_type, env, rho.iloc[0,1], pval.iloc[0,1]])

res = pd.DataFrame(res)
res.columns = res.iloc[0,:]
res = res.iloc[1:,:]
res.sort_values(by='rho', ascending=False, inplace=True)
res.to_csv("Scripts/Data_Vis/Section_4/Table_S4_gini_vs_shap_rank_per_corr.tsv", sep="\t", index=False)


## Hypothesis test to prove ranking distributions are not random
# SH said this is not necessary, we only do this when we have a point estimate and are comparing it to a distribution of means
n = 1000000 # number of repetitions
res = pd.read_csv("Scripts/Data_Vis/Section_4/Table_S4_gini_vs_shap_rank_per_corr.tsv", sep="\t")
res["D_gini_mean"] = 0
res["D_gini_sd"] = 0
res["pval_gini_mean"] = 0
res["pval_gini_sd"] = 0
res["D_shap_mean"] = 0
res["D_shap_sd"] = 0
res["pval_shap_mean"] = 0
res["pval_shap_sd"] = 0

def get_rand_ranks(j, i):
	""" Function to get random ranks and perform a KS test to compare to actual ranks"""
	row = res.iloc[i,:]
	# Read in optimized model (after feature selection) rank percentile datasets
	gini_rank = dt.fread(f"Scripts/Data_Vis/Section_4/RF_FS_imp_{row['Data Type']}_rank_per.tsv").to_pandas()
	shap_rank = dt.fread(f"Scripts/Data_Vis/Section_4/RF_FS_shap_{row['Data Type']}_rank_per.tsv").to_pandas()
	gini_rank.set_index("C0", inplace=True)
	shap_rank.set_index("C0", inplace=True)
	actual_rank_gini = gini_rank.loc[:,row['Env']].dropna() # Drop extra features from other envs
	actual_rank_shap = shap_rank.loc[:,row['Env']].dropna()

	# Generate random ranks
	rand_samp_gini = pd.Series(np.random.random_sample(size=len(actual_rank_gini)))
	rand_samp_shap = pd.Series(np.random.random_sample(size=len(actual_rank_shap)))
	rand_samp_gini.index = actual_rank_gini.index
	rand_samp_shap.index = actual_rank_shap.index

	## KS-test
	D_gini, p_gini = ks_2samp(actual_rank_gini, rand_samp_gini, alternative="greater")
	D_shap, p_shap = ks_2samp(actual_rank_shap, rand_samp_shap, alternative="greater")
	
	return (D_gini, p_gini, D_shap, p_shap)

pool = multiprocessing.Pool(40)
results = []
for i in tqdm(range(len(res))):
	row = res.iloc[i,:]
	## Generate random ranks and perform KS test
	result = pool.map(partial(get_rand_ranks, i=i), range(n), chunksize=1)
	results.append(result)
	
	## Unpack KS test results and add to res table
	gini_D, gini_p, shap_D, shap_p = zip(*results[i])
	res.iloc[i, res.columns.get_loc("D_gini_mean")] = np.mean(gini_D)
	res.iloc[i, res.columns.get_loc("D_gini_sd")] = np.std(gini_D)
	res.iloc[i, res.columns.get_loc("pval_gini_mean")] = np.mean(gini_p)
	res.iloc[i, res.columns.get_loc("pval_gini_sd")] = np.std(gini_p)
	res.iloc[i, res.columns.get_loc("D_shap_mean")] = np.mean(shap_D)
	res.iloc[i, res.columns.get_loc("D_shap_sd")] = np.std(shap_D)
	res.iloc[i, res.columns.get_loc("pval_shap_mean")] = np.mean(shap_p)
	res.iloc[i, res.columns.get_loc("pval_shap_sd")] = np.std(shap_p)

res.to_csv("Scripts/Data_Vis/Section_4/Table_S_gini_vs_shap_rank_per_corr_hypothesis_test_vfinal.tsv", sep="\t", index=False)

pool.close()
pool.join()
