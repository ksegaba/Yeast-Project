#!/usr/bin/env python3

import os
import re
import multiprocessing
import pandas as pd
import datatable as dt
import numpy as np
from scipy import stats
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
### TABLE S5
################################################################################
## SPEARMAN'S RHO CORRELATIONS BEWTEEN DATA TYPES
map_snps = pd.read_csv("Scripts/Data_Vis/biallelic_snps_diploid_and_S288C_genes.txt",
					   header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_16_removed.txt", sep="\t")

res = [["Model Type", "Importance Type", "Env", "Comparison", "rho", "pval"]]

for mod_type in ["baseline", "FS"]:
	for imp_type in ["imp", "shap"]:
		df_snp = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_snp_rank_per.tsv").to_pandas()
		df_pav = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_pav_rank_per.tsv").to_pandas()
		df_cnv = dt.fread(f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_cnv_rank_per.tsv").to_pandas()
		
		## Add gene information
		df_snp = map_snps[["snp", "gene"]].merge(df_snp, left_on="snp", right_on="C0", how="right")
		df_pav = map_orfs[["orf", "gene"]].merge(df_pav, left_on="orf", right_on="C0", how="right")
		df_cnv = map_orfs[["orf", "gene"]].merge(df_cnv, left_on="orf", right_on="C0", how="right")

		## Calculate spearman's rho
		if imp_type == "imp":
			imp_score_type = "mean gini importance"
		else:
			imp_score_type = "mean absolute SHAP value"
		
		# prep data for SNP to PAV/CNV comparisons
		tmp_snp = df_snp.loc[df_snp.gene != "intergenic",:] # drop intergenic snps
		tmp_pav = df_pav.loc[~df_pav.gene.isna(),:] # drop ORFs with no gene mapping
		tmp_cnv = df_cnv.loc[~df_cnv.gene.isna(),:]
		tmp_snp.drop(columns=["snp", "C0"], inplace=True) # drop snp column
		tmp_pav.drop(columns=["orf", "C0"], inplace=True) # drop orf and organism columns
		tmp_cnv.drop(columns=["orf", "C0"], inplace=True)
		tmp_snp = tmp_snp.groupby("gene").max() # get the max rank percentile per gene
		tmp_pav = tmp_pav.groupby("gene").max()
		tmp_cnv = tmp_cnv.groupby("gene").max()
		
		for env in mapping.keys():
			## SNP vs PAV
			try:
				df = pd.concat([tmp_snp.loc[:,env], tmp_pav.loc[:,env]], ignore_index=False, axis=1).dropna()
				if len(df) > 2:
					rho = df.corr(method=lambda x, y: spearmanr(x,y).statistic)
					pval = df.corr(method=lambda x, y: spearmanr(x,y).pvalue)
					res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs", \
							rho.iloc[0,1], pval.iloc[0,1]])
				else:
					res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs", \
							None, None])
				del df
			except:
				res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs", \
					None, None])
			## SNP vs CNV
			try:
				df = pd.concat([tmp_snp.loc[:,env], tmp_cnv.loc[:,env]], ignore_index=False, axis=1).dropna()
				if len(df) > 2:
					rho = df.corr(method=lambda x, y: spearmanr(x,y).statistic)
					pval = df.corr(method=lambda x, y: spearmanr(x,y).pvalue)
					res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs", \
							rho.iloc[0,1], pval.iloc[0,1]])
				else:
					res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs", \
							None, None])
				del df
			except:
				res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs", \
					None, None])
			## PAV vs CNV (keep all ORFs regardless of gene mapping)
			try:
				pav = df_pav.loc[~df_pav.orf.isna(),["orf", env]].set_index("orf")
				cnv = df_cnv.loc[~df_cnv.orf.isna(),["orf", env]].set_index("orf")
				df = pd.concat([pav, cnv], ignore_index=False, axis=1).dropna()
				if len(df) > 2:
					rho = df.corr(method=lambda x, y: spearmanr(x,y).statistic)
					pval = df.corr(method=lambda x, y: spearmanr(x,y).pvalue)
					res.append([mod_type, imp_score_type, env, "PAVs vs CNVs", \
							rho.iloc[0,1], pval.iloc[0,1]])
				else:
					res.append([mod_type, imp_score_type, env, "PAVs vs CNVs", \
							None, None])
				del df
			except:
				res.append([mod_type, imp_score_type, env, "PAVs vs CNVs", \
					None, None])
		del df_snp, df_pav, df_cnv, tmp_snp, tmp_pav, tmp_cnv
res = pd.DataFrame(res)
res.columns = res.iloc[0,:]
res = res.iloc[1:,:]
res.to_csv("Scripts/Data_Vis/Section_4/Table_S5_rank_per_corr_btwn_data_types.tsv", sep="\t", index=False)