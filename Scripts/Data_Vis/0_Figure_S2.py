#!/usr/bin/env python3
############################################################################
# Figure S2
############################################################################

""" Yeast Fitness GxE Manuscript Figures """
import os
import re
import math
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy
from matplotlib import cm
from matplotlib.colors import Normalize, ListedColormap, rgb2hex
from matplotlib.patches import Rectangle
from scipy.interpolate import interpn
from scipy.stats import binned_statistic, spearmanr, linregress, ks_2samp
from tqdm import tqdm

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")


    pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # fitness data
    h2 = pd.read_csv("Results/Heritability_h2_H2_sommer.csv", index_col=0) # narrow-sense heritability
    
    # Single-environment models
    pcs = pd.read_csv("Results/RESULTS_RF_PCs_sorted.txt", sep="\t", index_col=0)
    snp = pd.read_csv("Results/RESULTS_RF_SNPs_FS.txt", sep="\t", index_col=0)
    orf = pd.read_csv("Results/RESULTS_RF_ORFs_FS.txt", sep="\t", index_col=0)
    cnv = pd.read_csv("Results/RESULTS_RF_CNVs_FS.txt", sep="\t", index_col=0)
    single = pd.concat([pcs.r2_test, snp.r2_test, orf.r2_test, cnv.r2_test], ignore_index=False, axis=1)
    single.columns = ["PCs", "SNPs", "PAVs", "CNVs"]
    single.to_csv("Results/RESULTS_RF_ALL_FS.txt", sep="\t")

    # Multi-environment models
    brr_res = pd.read_csv("/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/RESULTS_GxE.txt", sep="\t")
    brr_res["ID"] = brr_res.apply(lambda x: x.ID.split("_")[1], axis=1)
    brr_res = brr_res.loc[brr_res.ID==brr_res.Env,:]
    dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/multi-trait/cv1"
    mo_res =  [f"{dir}/{file}" for file in os.listdir(dir) if file.endswith("_preds_test_multitrait.txt")]
    multi = brr_res[["Alg", "Env", "r2_test"]]
    from sklearn.metrics import r2_score
    for file in mo_res:
        Env =  file.split("_")[-8]
        df = pd.read_csv(file, sep="\t", index_col=[0,1])
        means = df.groupby("Trait").mean().T
        r2 = r2_score(pheno.loc[means.index,Env], means[Env])
        toplot = toplot._append({"Alg":"Multi-output", "Env":Env, "r2_test":r2}, ignore_index=True)
    multi.to_csv("Results/RESULTS_ALL_ALG_GxE_PAVs.txt", sep="\t", index=False)

    #### Plot Single- and Multi-environment models results
    single = pd.read_csv("Results/RESULTS_RF_ALL_FS.txt", sep="\t", index_col=0) # single-environment models
    multi = pd.read_csv("Results/RESULTS_ALL_ALG_GxE_PAVs.txt", sep="\t") # multi-environment models
    meta = pd.concat([pheno.median(), h2["h2"]], ignore_index=False, axis=1) # median fitness and narrow-sense heritability
    media = pd.Series([v[0] for v in pheno.columns.str.split(r"(?<=YP(?!D))|(?<=YPD).+", regex=True)], index=pheno.columns) # growth media
    envs = [v[1] for v in pheno.columns.str.split(r"^YP(?!D)|YPD", regex=True)] # environments
    envs = {pheno.columns[i]:envs[i] for i in range(len(envs))}
    meta = pd.concat([media, meta], ignore_index=False, axis=1)
    meta.columns = ["media", "median", "h2"]
    
    # cluster single-environment models and reorder meta and multi accordingly
    dend = hierarchy.dendrogram(hierarchy.linkage(single, method="average"))
    single = single.iloc[dend["leaves"],:]
    multi = multi.pivot(index="Env", columns="Alg", values="r2_test") # reshape multi-env matrix
    multi = multi.loc[single.index,:]
    meta = meta.loc[single.index,:]

    # plot heatmaps
    single.rename(index=envs, inplace=True)
    multi.rename(index=envs, inplace=True)
    meta.rename(index=envs, inplace=True)
    meta.media.replace({"YPD":0, "YP":1}, inplace=True) # recode meta.media to binary
    fig, ax = plt.subplots(1, 4, figsize=(10, 11), sharey=True)
    sns.heatmap(meta.iloc[:,:2], cmap="plasma", ax=ax[0], cbar=True, cbar_kws={"location":"bottom"}, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    sns.heatmap(meta.iloc[:,2:], cmap="magma", ax=ax[1], cbar=True, cbar_kws={"location":"bottom"}, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    sns.heatmap(single, cmap="viridis", ax=ax[2], cbar=True, cbar_kws={"location":"bottom"}, center=0.3, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    sns.heatmap(multi, cmap="viridis", ax=ax[3], cbar=True, cbar_kws={"location":"bottom"}, center=0.3, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    plt.tight_layout()
    plt.savefig("Scripts/Data_Vis/Figure_2_model_performances.pdf")
    plt.close()
