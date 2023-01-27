"""
Make feature map for Gustavo.
0 means the feature was lost after feature selection.
1 means the feature was kept.
"""
author = "Kenia Segura Aba"

import os
import pandas as pd
import datatable as dt
import glob

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/")
res = pd.read_csv("Results/RESULTS_RF_SNPs_FS.txt", sep="\t") # RF FS Results
geno = dt.fread("Data/Peter_2018/geno.csv") # to get SNP IDs

map = pd.DataFrame(index=geno.names, columns=res.cond.unique()) # binary map

for i in range(res.shape[0]):
    # get feature file
    f = glob.glob(
        f"/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs/feat*_rf_{res.cond[i]}_top_{res.FeatureNum[i]}")
    feat = open(f[0], "r").readlines()
    feat = [snp.strip() for snp in feat]
    map.loc[map.index.isin(feat), res.cond[i]] = 1

map.fillna(0, inplace=True)
map.sum() # confirmation
map.to_csv("Scripts/Genomic_Prediction_RKHS/data/feat_map_for_GDLC.csv", index=True)

del geno, res, map