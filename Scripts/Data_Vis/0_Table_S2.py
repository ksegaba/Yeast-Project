#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

################################################################################
### TABLE S2
################################################################################
# Correlations between performance of models trained on different data types
# Model results files:
res_base = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t")
res_fs = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_FS.txt", sep="\t")

############################ BASELINE MODELS FIRST #############################
# Compare between data types
pcc_btwn_data = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_base.loc[res_base.Alg==alg, ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_data = pd.concat([pcc_btwn_data, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_data = pd.concat([pcc_btwn_data, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")

pcc_btwn_data.columns = columns
pcc_btwn_data.to_csv("Scripts/Data_Vis/Section_2/Baseline_performance_correlations_between_data_types.tsv", index=True, sep="\t")

# Compare between algorithms
pcc_btwn_alg = pd.DataFrame()
columns = []
for data_type in ["PCs_tassel", "SNP", "PAV", "CNV"]:
    df = res_base.loc[res_base.Data==data_type, ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_alg = pd.concat([pcc_btwn_alg, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_alg = pd.concat([pcc_btwn_alg, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_alg.columns = columns
pcc_btwn_alg.to_csv("Scripts/Data_Vis/Section_2/Baseline_performance_correlations_between_algorithms.tsv", index=True, sep="\t")

# compare between environments across algorithms
pcc_btwn_env = pd.DataFrame()
columns = []
for data_type in ["PCs_tassel", "SNP", "PAV", "CNV"]:
    df = res_base.loc[(res_base.Data==data_type), ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv("Scripts/Data_Vis/Section_2/Baseline_performance_correlations_between_envs_across_data_types.tsv", index=True, sep="\t")

# compare between environments across data types
pcc_btwn_env = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_base.loc[(res_base.Alg==alg), ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv("Scripts/Data_Vis/Section_2/Baseline_performance_correlations_between_envs_across_algorithms.tsv", index=True, sep="\t")

########################### FEATURE SELECTION MODELS ###########################
# Compare between data types
pcc_btwn_data_fs = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_fs.loc[res_fs.Alg==alg, ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_data_fs = pd.concat([pcc_btwn_data_fs, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_data_fs = pd.concat([pcc_btwn_data_fs, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")

pcc_btwn_data_fs.columns = columns
pcc_btwn_data_fs.to_csv("Scripts/Data_Vis/Section_2/FS_performance_correlations_between_data_types.tsv", index=True, sep="\t")

# Compare between algorithms
pcc_btwn_alg_fs = pd.DataFrame()
columns = []
for data_type in ["SNP", "PAV", "CNV"]:
    df = res_fs.loc[res_fs.Data==data_type, ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_alg_fs = pd.concat([pcc_btwn_alg_fs, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_alg_fs = pd.concat([pcc_btwn_alg_fs, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_alg_fs.columns = columns
pcc_btwn_alg_fs.to_csv("Scripts/Data_Vis/Section_2/FS_performance_correlations_between_algorithms.tsv", index=True, sep="\t")

# compare between environments across algorithms
pcc_btwn_env = pd.DataFrame()
columns = []
for data_type in ["SNP", "PAV", "CNV"]:
    df = res_fs.loc[(res_fs.Data==data_type), ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv("Scripts/Data_Vis/Section_2/FS_performance_correlations_between_envs_across_data_types.tsv", index=True, sep="\t")

# compare between environments across data types
pcc_btwn_env = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_fs.loc[(res_fs.Alg==alg), ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x,y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x,y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv("Scripts/Data_Vis/Section_2/FS_performance_correlations_between_envs_across_algorithms.tsv", index=True, sep="\t")
