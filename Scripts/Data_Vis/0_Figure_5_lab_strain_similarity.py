#!/usr/bin/env python3
"""
Figure 5: SHAP values of top genes/benchmark genes VS genetic similarity to
S288C or W303 (SACE_GAV).
"""

import os, glob, re, swifter, tqdm, threading
import pandas as pd
import numpy as np
import datatable as dt
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import cupy as cp
import statsmodels.api as sm
from functools import partial
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial.distance import euclidean
from scipy.stats import mannwhitneyu

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

# Read in the benchmark gene and gene map data
ben_orf = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes_orfs.txt")
ben_snp = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes_snps.txt")
caf_orf = pd.read_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes_orfs.txt")
caf_snp = pd.read_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes_snps.txt")
cu_orf = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes_orfs.txt")
cu_snp = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes_snps.txt")
sma_orf = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes_orfs.txt")
sma_snp = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes_snps.txt")

map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                       sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")
map_orfs.loc[map_orfs.organism=="Saccharomyces cerevisiae", "organism"] = "Saccharomyces cerevisiae S288C"

""" Code to generate the genetic distance matrices from SNP and PAV data
# Read in the SNP, PAV, and test set data
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
geno = dt.fread("Data/Peter_2018/geno.csv").to_pandas()
geno.set_index("ID", inplace=True)

# Add S288C SNP genotypes as a row of -1s (homozygous for the reference allele)
geno = geno.T
geno["S288C"] = -1
geno = geno.T
# geno.to_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_with_S288C.csv")
geno = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_with_S288C.csv").to_pandas()
geno.set_index("ID", inplace=True)

# Calculate the SNP-based euclidean distance and the SNP-based jaccard distance
geno_train = geno.loc[~geno.index.isin(test[0]),:] # remove the test set before calculating genetic distances
eu_dist_snp = pdist(geno_train.values, metric="euclidean")
eu_dist_snp = pd.DataFrame(squareform(eu_dist_snp), columns=geno_train.index, index=geno_train.index)
ja_dist_snp = pdist(geno_train.values, metric="jaccard") # proportion of dissimilarity
ja_dist_snp = pd.DataFrame(squareform(ja_dist_snp), columns=geno_train.index, index=geno_train.index)
# eu_dist_snp.to_csv("Scripts/Data_Vis/Section_5/genetic_distance_snp_euclidean_to_S288C.csv")
# ja_dist_snp.to_csv("Scripts/Data_Vis/Section_5/genetic_distance_snp_jaccard_to_S288C.csv")
eu_dist_snp = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_snp_euclidean_to_S288C.csv", index_col=0)
ja_dist_snp = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_snp_jaccard_to_S288C.csv", index_col=0)

# Determine which ORFs are present in the reference genome
os.system("grep '>' Data/S288C_reference_genome_R64-3-1_20210421/orf_coding_all_R64-3-1_20210421.fasta | wc -l") # 6716 reference ORFs
# os.system("grep '>' Data/S288C_reference_genome_R64-3-1_20210421/orf_coding_all_R64-3-1_20210421.fasta | \
#           cut -d ' ' -f1 | cut -d '>' -f2 > \
#           Data/S288C_reference_genome_R64-3-1_20210421/orf_coding_all_R64-3-1_20210421_IDs_only.txt") # extract reference ORF IDs
ref_orf_ids = pd.read_csv("Data/S288C_reference_genome_R64-3-1_20210421/orf_coding_all_R64-3-1_20210421_IDs_only.txt", header=None)

pav = dt.fread("Data/Peter_2018/ORFs_pres_abs.csv", max_nrows=1, header=False) # Peter et al., 2018 ORFs
pav = pav[:,1:].to_pandas().T # exclude 'ID' column
pav[1] = pav.replace("X[0-9]+\.", "", regex=True) # remove the ORF ID prefix
pav[1] = pav[1].str.split("_", expand=True)[0] # get the genes, non-reference ORFs have identifiers that do not match the systematic gene names
pav[1] = pav[1].str.replace(".", "-")

# in_ref_orfs = np.intersect1d(pav[1], ref_orf_ids[0]) # based on pattern matching the ORF IDs in S288C reference genome
# len(in_ref_orfs) # 6059, I'm not sure if I can trust this, since the BLAST resulted in less ORFs, plus some ORFs mapped to a different gene than what is in the ORF identifier

pav["orf"] = pav.apply(lambda x: re.sub("^X", "", x[0]), axis=1) # fix orf IDs
pav["orf"] = pav.apply(lambda x: re.sub("\.", "-", x["orf"]), axis=1)
pav = pav.merge(map_orfs[map_orfs.organism=="Saccharomyces cerevisiae S288C"],
          left_on="orf", right_on="orf", how="left") # ORF to gene map is based on BLAST results
pav["in_ref_blast"] = pav.apply(lambda x: 1 if not pd.isna(x["gene"]) else 0, axis=1)
in_ref_orfs = pav.loc[pav["in_ref_blast"]==1, "gene"].to_list()

# sanity check
# >>> pav["organism"].unique()
# array([nan, 'Saccharomyces cerevisiae S288C'], dtype=object)
# >>> pav.apply(lambda x: 1 if not pd.isna(x["gene"]) else 0, axis=1).value_counts()
# 1    5561
# 0    2147
# Name: count, dtype: int64
# >>> pav.apply(lambda x: 1 if x["organism"]=='Saccharomyces cerevisiae S288C' else 0, axis=1).value_counts()
# 1    5561
# 0    2147
# Name: count, dtype: int64

len(np.intersect1d(in_ref_orfs, pav.gene.to_list())) # 5533 I WILL USE THE BLAST RESULTS TO CALCULATE THE GENETIC DISTANCES FOR PAV

## Calculate the PAV-based euclidean distance and the PAV-based jaccard distance
pav_df = dt.fread("Data/Peter_2018/ORFs_pres_abs.csv").to_pandas() # read in PAV genotype data
pav_df.set_index("ID", inplace=True)
pav_df = pav_df.T.merge(pav.set_index(0)["in_ref_blast"], left_index=True, right_index=True) # add the S288C PAV genotypes
pav_df.rename(columns={"in_ref_blast": "in_S288C"}, inplace=True)
pav_df.to_csv("Data/Peter_2018/ORFs_pres_abs_with_S288C.csv")

pav_train = pav_df.loc[:,~pav_df.columns.isin(test[0])] # remove the test set before calculating genetic distances
pav_train = pav_train.astype(int)
eu_dist_pav = pdist(pav_train.T.values, metric="euclidean")
eu_dist_pav = pd.DataFrame(squareform(eu_dist_pav), columns=pav_train.columns, index=pav_train.columns)
ja_dist_pav = pdist(pav_train.T.values, metric="jaccard") # proportion of dissimilarity
ja_dist_pav = pd.DataFrame(squareform(ja_dist_pav), columns=pav_train.columns, index=pav_train.columns)
# eu_dist_pav.to_csv("Scripts/Data_Vis/Section_5/genetic_distance_pav_euclidean_to_S288C.csv")
# ja_dist_pav.to_csv("Scripts/Data_Vis/Section_5/genetic_distance_pav_jaccard_to_S288C.csv")
"""

# Read in the distance matrices
eu_dist_snp = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_snp_euclidean_to_S288C.csv", index_col=0)
ja_dist_snp = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_snp_jaccard_to_S288C.csv", index_col=0)
eu_dist_pav = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_pav_euclidean_to_S288C.csv", index_col=0)
ja_dist_pav = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_pav_jaccard_to_S288C.csv", index_col=0)

# Ensure isolates are all in the same order
ja_dist_snp = ja_dist_snp.loc[eu_dist_snp.index, eu_dist_snp.columns]
eu_dist_pav.rename(columns={"in_S288C": "S288C"}, index={"in_S288C": "S288C"}, inplace=True)
eu_dist_pav = eu_dist_pav.loc[eu_dist_snp.index, eu_dist_snp.columns]
ja_dist_pav.rename(columns={"in_S288C": "S288C"}, index={"in_S288C": "S288C"}, inplace=True)
ja_dist_pav = ja_dist_pav.loc[eu_dist_snp.index, eu_dist_snp.columns]

#### Plot the benchmark gene shap values and distance value for each isolate ###
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
            "YPDSODIUMMETAARSENITE"]
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP"

def plot_benchmark_shap_vs_distance_scatter(df, x, y, model, save="", **kwargs):
    """ Create a scatter plot. This function is called within lin_reg_stats()
    to plot shap values vs distance values.
    """
    plt.figure(figsize=(10, 4))
    sns.scatterplot(data=df, x=x, y=y, **kwargs)
    
    # draw the regression line
    print(model.params)
    coef = model.params
    print(np.array(df[x]))
    plt.plot(df[x], coef[0]*np.array(df[x]), color="red")
    
    plt.savefig(save, bbox_inches="tight", dpi=300)
    plt.close()
    
    return None


def lin_reg_stats(shap_df, dist_df, target_col, save, **kwargs):
    """ Ordinary Least Squares Regression
    Returns:
        model.params: coefficients
        model.rsquared: R-squared value
        model.pvalues: p-values for each coefficient
    """
    
    # flatten the shap value matrix and map to distance values
    if target_col=="S288C":
        y = np.array(shap_df.iloc[:,1:].astype("float").loc[:,dist_df.index[:625]]).flatten()
        x = np.tile(dist_df["S288C"][:625].to_numpy(), len(shap_df))
    
    if target_col=="W303":
        cols = dist_df.drop(["SACE_GAV", "S288C"], axis=1).columns
        y = np.array(shap_df.drop("SACE_GAV", axis=1).\
            iloc[:,1:].astype("float").loc[:,cols]).flatten()
        x = np.tile(dist_df["SACE_GAV"][cols].to_numpy(), len(shap_df))
        
    print(y.shape, x.shape)
    
    # fit the linear model
    model = sm.OLS(endog=y, exog=x).fit()
    print(model.summary())
    summary = model.summary2()
    summary_table = summary.tables[1].to_dict()
    summary_table["R2"] = {"x1": model.rsquared}
    
    # plot the regression
    df = pd.DataFrame([x,y], index=["distance", "shap"]).T
    plot_benchmark_shap_vs_distance_scatter(df, "distance", "shap", model, save, **kwargs)
    
    return summary_table


lin_res_ja = {} # collect the linear regression results when jaccard distance is considered for S288C
lin_res_eu = {} # when euclidean distance is considered for S288C
lin_res_ja_w303 = {}
lin_res_eu_w303 = {}
for env in target_envs:
    print(env)
    # read in the shap values
    snp_shap_all = dt.fread(f"{d}/SNP/baseline/SHAP_values_sorted_{env}_snp_rf_baseline_training.txt").to_pandas()
    pav_shap_all = dt.fread(f"{d}/PAV/baseline/SHAP_values_sorted_{env}_pav_rf_baseline_training.txt").to_pandas()
    pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix orf IDs
    pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
    
    # map genes to shap data
    snp_shap_all = snp_shap_all.T.reset_index()
    pav_shap_all = pav_shap_all.T.reset_index()
    snp_shap_all.columns = snp_shap_all.iloc[0,:] # set column names
    snp_shap_all = snp_shap_all.iloc[1:,:] # remove first row
    pav_shap_all.columns = pav_shap_all.iloc[0,:]
    pav_shap_all = pav_shap_all.iloc[1:,:]
    snp_shap_all["ID"] = snp_shap_all["ID"].map(map_snps.set_index("snp")["gene"]) # map features to genes
    pav_shap_all["ID"] = pav_shap_all["ID"].map(map_orfs.set_index("orf")["gene"])
    
    # Ensure isolate columns are the same order as the distance matrices
    snp_shap_all = snp_shap_all.loc[:,["ID"] + eu_dist_snp.columns.tolist()[:-1]]
    pav_shap_all = pav_shap_all.loc[:,["ID"] + eu_dist_pav.columns.tolist()[:-1]]
    
    if env == "YPDBENOMYL500":
        # subset the benchmark genes from the shap data
        ben_snp_sub = snp_shap_all.loc[snp_shap_all["ID"].isin(ben_snp["0"]),:]
        ben_snp_sub.iloc[:,1:] = ben_snp_sub.iloc[:,1:].abs() # take the absolute value of the shap values, only magnitude represents importance
        ben_snp_sub_max = ben_snp_sub.groupby("ID").max().reset_index() # set the max shap value per gene
        
        ben_pav_sub = pav_shap_all.loc[pav_shap_all["ID"].isin(ben_orf["0"]),:]
        ben_pav_sub.iloc[:,1:] = ben_pav_sub.iloc[:,1:].abs()
        ben_pav_sub_max = ben_pav_sub.groupby("ID").max().reset_index()
        
        # calculate linear regression statistics for each benchmark gene
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_benomyl_snp_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_benomyl_pav_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        lin_res_ja[env] = {"snp": lin_reg_stats(ben_snp_sub_max, ja_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(ben_pav_sub_max, ja_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_benomyl_snp_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_benomyl_pav_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        lin_res_eu[env] = {"snp": lin_reg_stats(ben_snp_sub_max, eu_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(ben_pav_sub_max, eu_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_benomyl_snp_shap_vs_jaccard_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_benomyl_pav_shap_vs_jaccard_distance_to_W303_v2.pdf"
        lin_res_ja_w303[env] = {"snp": lin_reg_stats(ben_snp_sub_max, ja_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(ben_pav_sub_max, ja_dist_pav, "W303", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_benomyl_snp_shap_vs_euclidean_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_benomyl_pav_shap_vs_euclidean_distance_to_W303_v2.pdf"
        lin_res_eu_w303[env] = {"snp": lin_reg_stats(ben_snp_sub_max, eu_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(ben_pav_sub_max, eu_dist_pav, "W303", save_pav)}
        
    if (env == "YPDCAFEIN40") or (env == "YPDCAFEIN50"):
        caf_snp_sub = snp_shap_all.loc[snp_shap_all["ID"].isin(caf_snp["0"]),:]
        caf_snp_sub.iloc[:,1:] = caf_snp_sub.iloc[:,1:].abs() # take the absolute value of the shap values, only magnitude represents importance
        caf_snp_sub_max = caf_snp_sub.copy(deep=True).groupby("ID").max().reset_index() # set the max shap value per gene
        
        caf_pav_sub = pav_shap_all.loc[pav_shap_all["ID"].isin(caf_orf["0"]),:]
        caf_pav_sub.iloc[:,1:] = caf_pav_sub.iloc[:,1:].abs()
        caf_pav_sub_max = caf_pav_sub.copy(deep=True).groupby("ID").max().reset_index()
        
        # calculate linear regression statistics for each benchmark gene
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_caffeine_snp_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_caffeine_pav_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        lin_res_ja[env] = {"snp": lin_reg_stats(caf_snp_sub_max, ja_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(caf_pav_sub_max, ja_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_caffeine_snp_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_caffeine_pav_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        lin_res_eu[env] = {"snp": lin_reg_stats(caf_snp_sub_max, eu_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(caf_pav_sub_max, eu_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_caffeine_snp_shap_vs_jaccard_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_caffeine_pav_shap_vs_jaccard_distance_to_W303_v2.pdf"
        lin_res_ja_w303[env] = {"snp": lin_reg_stats(caf_snp_sub_max, ja_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(caf_pav_sub_max, ja_dist_pav, "W303", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_caffeine_snp_shap_vs_euclidean_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_caffeine_pav_shap_vs_euclidean_distance_to_W303_v2.pdf"
        lin_res_eu_w303[env] = {"snp": lin_reg_stats(caf_snp_sub_max, eu_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(caf_pav_sub_max, eu_dist_pav, "W303", save_pav)}
        
    if env == "YPDCUSO410MM":
        cu_snp_sub = snp_shap_all.loc[snp_shap_all["ID"].isin(cu_snp["0"]),:]
        cu_snp_sub.iloc[:,1:] = cu_snp_sub.iloc[:,1:].abs() # take the absolute value of the shap values, only magnitude represents importance
        cu_snp_sub_max = cu_snp_sub.groupby("ID").max().reset_index() # set the max shap value per gene
        
        cu_pav_sub = pav_shap_all.loc[pav_shap_all["ID"].isin(cu_orf["0"]),:]
        cu_pav_sub.iloc[:,1:] = cu_pav_sub.iloc[:,1:].abs()
        cu_pav_sub_max = cu_pav_sub.groupby("ID").max().reset_index()
        
        # calculate linear regression statistics for each benchmark gene
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_cuso4_snp_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_cuso4_pav_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        lin_res_ja[env] = {"snp": lin_reg_stats(cu_snp_sub_max, ja_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(cu_pav_sub_max, ja_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_cuso4_snp_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_cuso4_pav_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        lin_res_eu[env] = {"snp": lin_reg_stats(cu_snp_sub_max, eu_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(cu_pav_sub_max, eu_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_cuso4_snp_shap_vs_jaccard_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_cuso4_pav_shap_vs_jaccard_distance_to_W303_v2.pdf"
        lin_res_ja_w303[env] = {"snp": lin_reg_stats(cu_snp_sub_max, ja_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(cu_pav_sub_max, ja_dist_pav, "W303", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_cuso4_snp_shap_vs_euclidean_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_cuso4_pav_shap_vs_euclidean_distance_to_W303_v2.pdf"
        lin_res_eu_w303[env] = {"snp": lin_reg_stats(cu_snp_sub_max, eu_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(cu_pav_sub_max, eu_dist_pav, "W303", save_pav)}
        
    if env == "YPDSODIUMMETAARSENITE":
        sma_snp_sub = snp_shap_all.loc[snp_shap_all["ID"].isin(sma_snp["0"]),:]
        sma_snp_sub.iloc[:,1:] = sma_snp_sub.iloc[:,1:].abs() # take the absolute value of the shap values, only magnitude represents importance
        sma_snp_sub_max = sma_snp_sub.groupby("ID").max().reset_index() # set the max shap value per gene
        
        sma_pav_sub = pav_shap_all.loc[pav_shap_all["ID"].isin(sma_orf["0"]),:]
        sma_pav_sub.iloc[:,1:] = sma_pav_sub.iloc[:,1:].abs()
        sma_pav_sub_max = sma_pav_sub.groupby("ID").max().reset_index()
        
        # calculate linear regression statistics for each benchmark gene
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_snp_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_pav_shap_vs_jaccard_distance_to_S288C_v2.pdf"
        lin_res_ja[env] = {"snp": lin_reg_stats(sma_snp_sub_max, ja_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(sma_pav_sub_max, ja_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_snp_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_pav_shap_vs_euclidean_distance_to_S288C_v2.pdf"
        lin_res_eu[env] = {"snp": lin_reg_stats(sma_snp_sub_max, eu_dist_snp, "S288C", save_snp),
                       "pav": lin_reg_stats(sma_pav_sub_max, eu_dist_pav, "S288C", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_snp_shap_vs_jaccard_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_pav_shap_vs_jaccard_distance_to_W303_v2.pdf"
        lin_res_ja_w303[env] = {"snp": lin_reg_stats(sma_snp_sub_max, ja_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(sma_pav_sub_max, ja_dist_pav, "W303", save_pav)}
        
        save_snp = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_snp_shap_vs_euclidean_distance_to_W303_v2.pdf"
        save_pav = f"Scripts/Data_Vis/Section_5/{env}_sodium_meta-arsenite_pav_shap_vs_euclidean_distance_to_W303_v2.pdf"
        lin_res_eu_w303[env] = {"snp": lin_reg_stats(sma_snp_sub_max, eu_dist_snp, "W303", save_snp),
                       "pav": lin_reg_stats(sma_pav_sub_max, eu_dist_pav, "W303", save_pav)}


# reshape lin_reg_res to have the environment as a column and statistics as rows
df = pd.json_normalize(lin_res_ja, sep="_")
df = df.transpose()
df.index = pd.MultiIndex.from_tuples(df.index.str.split("_").map(tuple))
df.index.names = ["Env", "DataType", "Statistics", "Vals"]
df = df.pivot_table(index=["DataType", "Statistics"], columns="Env", values=0)
df.to_csv("Scripts/Data_Vis/Section_5/Table_S_benchmark_gene_shap_vs_jaccard_distance_to_S288C_lm_stats_v2.csv")

df = pd.json_normalize(lin_res_eu, sep="_")
df = df.transpose()
df.index = pd.MultiIndex.from_tuples(df.index.str.split("_").map(tuple))
df.index.names = ["Env", "DataType", "Statistics", "Vals"]
df = df.pivot_table(index=["DataType", "Statistics"], columns="Env", values=0)
df.to_csv("Scripts/Data_Vis/Section_5/Table_S_benchmark_gene_shap_vs_euclidean_distance_to_S288C_lm_stats_v2.csv")

df = pd.json_normalize(lin_res_ja_w303, sep="_")
df = df.transpose()
df.index = pd.MultiIndex.from_tuples(df.index.str.split("_").map(tuple))
df.index.names = ["Env", "DataType", "Statistics", "Vals"]
df = df.pivot_table(index=["DataType", "Statistics"], columns="Env", values=0)
df.to_csv("Scripts/Data_Vis/Section_5/Table_S_benchmark_gene_shap_vs_jaccard_distance_to_W303_lm_stats_v2.csv")

df = pd.json_normalize(lin_res_eu_w303, sep="_")
df = df.transpose()
df.index = pd.MultiIndex.from_tuples(df.index.str.split("_").map(tuple))
df.index.names = ["Env", "DataType", "Statistics", "Vals"]
df = df.pivot_table(index=["DataType", "Statistics"], columns="Env", values=0)
df.to_csv("Scripts/Data_Vis/Section_5/Table_S_benchmark_gene_shap_vs_euclidean_distance_to_W303_lm_stats_v2.csv")

################################################################################
"""Compare the fitness distributions and shap values between the cluster
containing the lab strains and the most distinct cluster to it."""

# Genetic distance matrices
eu_dist_snp = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_snp_euclidean_to_S288C.csv",
                          index_col=0)
eu_dist_pav = pd.read_csv("Scripts/Data_Vis/Section_5/genetic_distance_pav_euclidean_to_S288C.csv",
                          index_col=0)
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None) # to get training instances

# Clade & fitness information
clades = pd.read_excel("Data/Peter_2018/0_raw_data/Peter_2018_Supplementary_Tables.xls",
                       sheet_name="Table S1", skiprows=3, nrows=1011) # isolate clades
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # isolate fitness

# Map clades to isolates and create a colormap
clades = clades[["Standardized name", "Clades"]] # subset relevant columns
clades = clades.loc[clades["Standardized name"].isin(pheno.index)] # diploid isolates
clades.set_index("Standardized name", inplace=True)
clades.loc[clades.Clades.isna(),"Clades"] = "Unknown" # replace NaN with Unknown
clades.loc["S288C"] = "Reference" # insert a row for "S288C"
clades.loc[["S288C", "SACE_GAV"], "Clades"]
# S288C                 Reference
# SACE_GAV    M3. Mosaic region 3

########### Principal Component Analysis on genetic distance matrices
snp_train = eu_dist_snp.loc[~eu_dist_snp.index.isin(test[0]), ~eu_dist_snp.index.isin(test[0])]
pav_train = eu_dist_pav.loc[~eu_dist_pav.index.isin(test[0]), ~eu_dist_pav.index.isin(test[0])]

pca = PCA(n_components=5)
pca_snp = pca.fit(snp_train)
pca_snp_df = pca_snp.transform(snp_train)
vexp_pca_snp = pca_snp.explained_variance_ratio_ # variance explained by each component

pca_p = PCA(n_components=5)
pca_pav = pca_p.fit(pav_train)
pca_pav_df = pca_p.transform(pav_train)
vexp_pca_pav = pca_pav.explained_variance_ratio_

pca_snp_df = pd.DataFrame(pca_snp_df, index=snp_train.index)
pca_pav_df = pd.DataFrame(pca_pav_df, index=pav_train.index)
pca_snp_df.columns = [f"PC{i+1}" for i in range(5)]
pca_pav_df.columns = [f"PC{i+1}" for i in range(5)]
pca_pav_df.rename(index={"in_S288C": "S288C"}, inplace=True)

# Apply K-means clustering to the PCA results
inertia_snp = []
inertia_pav = []
for k in range(2, 11):
    kmeans_snp = KMeans(n_clusters=k, random_state=42).fit(pca_snp_df)
    kmeans_pav = KMeans(n_clusters=k, random_state=42).fit(pca_pav_df)
    inertia_snp.append(kmeans_snp.inertia_)
    inertia_pav.append(kmeans_pav.inertia_)

# Plot the elbow plot
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
ax[0].plot([k for k in range(2, 11)], inertia_snp, marker="o")
ax[0].set_title("Elbow Plot for PCA of SNP genetic distance")
ax[0].set_xlabel("Number of clusters (k)")
ax[0].set_ylabel("Inertia")
ax[1].plot([k for k in range(2, 11)], inertia_pav, marker="o")
ax[1].set_title("Elbow Plot for PCA of PAV genetic distance")
ax[1].set_xlabel("Number of clusters (k)")
ax[1].set_ylabel("Inertia")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/PCA_snp_or_pav_genetic_distance_elbow_plot.pdf",
            bbox_inches="tight", dpi=300)
plt.close()

# Fit K-means with the optimal number of clusters
kmeans_snp = KMeans(n_clusters=6, random_state=42).fit(pca_snp_df)
kmeans_pav = KMeans(n_clusters=4, random_state=42).fit(pca_pav_df)

# Plot the PCA of eu_dist_snp and eu_dist_pav for the training data only (shap values were calculated from the training data)
fig, ax = plt.subplots(2, 2, figsize=(12, 12))
sns.scatterplot(x=pca_snp_df["PC1"], y=pca_snp_df["PC2"], hue=kmeans_snp.labels_,
                ax=ax[0][0], palette="tab10", alpha=0.7, edgecolor="none")
sns.scatterplot(x=pca_pav_df["PC1"], y=pca_pav_df["PC2"], hue=kmeans_pav.labels_,
                ax=ax[0][1], palette="tab10", alpha=0.7, edgecolor="none")
sns.scatterplot(x=pca_snp_df["PC1"], y=pca_snp_df["PC2"],
                hue=clades.loc[pca_snp_df.index, "Clades"],
                ax=ax[1][0], palette="tab20", alpha=0.7, edgecolor="none")
sns.scatterplot(x=pca_pav_df["PC1"], y=pca_pav_df["PC2"],
                hue=clades.loc[pca_pav_df.index, "Clades"],
                ax=ax[1][1], palette="tab20", alpha=0.7, edgecolor="none")
for i in range(2):
    ax[i][0].annotate("S288C", (pca_snp_df.loc["S288C","PC1"],
                                pca_snp_df.loc["S288C","PC2"]), color="red")
    ax[i][0].annotate("W303", (pca_snp_df.loc["SACE_GAV","PC1"],
                               pca_snp_df.loc["SACE_GAV","PC2"]), color="red")
    ax[i][1].annotate("S288C", (pca_pav_df.loc["S288C","PC1"],
                                pca_pav_df.loc["S288C","PC2"]), color="red")
    ax[i][1].annotate("W303", (pca_pav_df.loc["SACE_GAV","PC1"],
                               pca_pav_df.loc["SACE_GAV","PC2"]), color="red")

for axis in ax.flat:
    axis.set_box_aspect(1)

plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/PCA_snp_or_pav_genetic_distance.pdf", bbox_inches="tight", dpi=300)
plt.close()

# Calculate distance between the cluster centers
def cluster_distance(centroids):
    distances = np.zeros((len(centroids), len(centroids)))
    for i in range(len(centroids)):
        for j in range(i + 1, len(centroids)):
            distances[i, j] = euclidean(centroids[i], centroids[j])
            distances[j, i] = distances[i, j]
    return pd.DataFrame(distances)


centroids_snp = kmeans_snp.cluster_centers_ # Cluster centers
centroids_pav = kmeans_pav.cluster_centers_
dist_snp = cluster_distance(centroids_snp)
dist_pav = cluster_distance(centroids_pav)

# Identify the clusters most distinct to the cluster in which S288C and W303 are in
cluster_assignments_snp = pd.DataFrame(kmeans_snp.labels_, index=pca_snp_df.index,
                                       columns=["Cluster"]) # S288C is in cluster 0
cluster_assignments_pav = pd.DataFrame(kmeans_pav.labels_, index=pca_pav_df.index,
                                       columns=["Cluster"]) # S288C is in cluster 1
cluster_assignments_snp.loc["SACE_GAV",:] # W303 is in cluster 5
cluster_assignments_pav.loc["SACE_GAV",:] # W303 is in cluster 0
s288c_distinct_clusters_snp = [0, dist_snp.loc[0,:].idxmax()]
s288c_distinct_clusters_pav = [1, dist_pav.loc[1,:].idxmax()]
w303_distinct_clusters_snp = [5, dist_snp.loc[5,:].idxmax()]
w303_distinct_clusters_pav = [0, dist_pav.loc[0,:].idxmax()]

######### Calculate the mann-whitney u test statistic of shap values between clusters
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # fitness data

# Benchmark genes validated in S288C or W303
ben_meta = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf_meta = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt", sep="\t").to_pandas()
cu_meta = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma_meta = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")

ben_s288c = ben_meta.loc[ben_meta["Strain Background"] == "S288C", "Gene Systematic Name"].unique()
ben_w303 = ben_meta.loc[ben_meta["Strain Background"] == "W303", "Gene Systematic Name"].unique()
caf_s288c = caf_meta.loc[caf_meta["Strain Background"] == "S288C", "Gene Systematic Name"].unique()
caf_w303 = caf_meta.loc[caf_meta["Strain Background"] == "W303", "Gene Systematic Name"].unique()
cu_s288c = cu_meta.loc[cu_meta["Strain Background"] == "S288C", "Gene Systematic Name"].unique()
cu_w303 = cu_meta.loc[cu_meta["Strain Background"] == "W303", "Gene Systematic Name"].unique()
sma_s288c = sma_meta.loc[sma_meta["Strain Background"] == "S288C", "Gene Systematic Name"].unique()
sma_w303 = sma_meta.loc[sma_meta["Strain Background"] == "W303", "Gene Systematic Name"].unique()

# Benchmark genes present in SNP and ORF data
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")

ben_snp = map_snps.loc[map_snps["Benomyl"] == 1, "gene"].unique()
ben_orf = map_orfs.loc[map_orfs["Benomyl"] == 1, "gene"].unique()
caf_snp = map_snps.loc[map_snps["Caffeine"] == 1, "gene"].unique()
caf_orf = map_orfs.loc[map_orfs["Caffeine"] == 1, "gene"].unique()
cu_snp = map_snps.loc[map_snps["CuSO4"] == 1, "gene"].unique()
cu_orf = map_orfs.loc[map_orfs["CuSO4"] == 1, "gene"].unique()
sma_snp = map_snps.loc[map_snps["Sodium_meta-arsenite"] == 1, "gene"].unique()
sma_orf = map_orfs.loc[map_orfs["Sodium_meta-arsenite"] == 1, "gene"].unique()

def compared_rand_clusters(k, env, rand_res, snp_shap_bench_s288c, snp_shap_bench_w303,
    snp_s288c_rclust1, snp_s288c_rclust2, snp_w303_rclust1, snp_w303_rclust2,
    pav_shap_bench_s288c, pav_shap_bench_w303, pav_s288c_rclust1,
    pav_s288c_rclust2, pav_w303_rclust1, pav_w303_rclust2):
    
    # Compare the shap values between two clusters with random sets of individuals
    snp_s288c_clust1 = snp_shap_bench_s288c.loc[snp_s288c_rclust1[k]].astype("float").to_numpy().flatten()
    snp_s288c_clust2 = snp_shap_bench_s288c.loc[snp_s288c_rclust2[k]].astype("float").to_numpy().flatten()
    pav_s288c_clust1 = pav_shap_bench_s288c.loc[pav_s288c_rclust1[k]].astype("float").to_numpy().flatten()
    pav_s288c_clust2 = pav_shap_bench_s288c.loc[pav_s288c_rclust2[k]].astype("float").to_numpy().flatten()
    if env != "YPDCUSO410MM":
        snp_w303_clust1 = snp_shap_bench_w303.loc[snp_w303_rclust1[k]].astype("float").to_numpy().flatten()
        snp_w303_clust2 = snp_shap_bench_w303.loc[snp_w303_rclust2[k]].astype("float").to_numpy().flatten()
        pav_w303_clust1 = pav_shap_bench_w303.loc[pav_w303_rclust1[k]].astype("float").to_numpy().flatten()
        pav_w303_clust2 = pav_shap_bench_w303.loc[pav_w303_rclust2[k]].astype("float").to_numpy().flatten()
        s2, p2 = mannwhitneyu(snp_w303_clust1, snp_w303_clust2, alternative="greater")
        s4, p4 = mannwhitneyu(pav_w303_clust1, pav_w303_clust2, alternative="greater")
        del snp_w303_clust1, snp_w303_clust2, pav_w303_clust1, pav_w303_clust2
    else:
        s2 = np.nan ; p2 = np.nan
        s4 = np.nan ; p4 = np.nan
        
    s1, p1 = mannwhitneyu(snp_s288c_clust1, snp_s288c_clust2, alternative="greater")
    s3, p3 = mannwhitneyu(pav_s288c_clust1, pav_s288c_clust2, alternative="greater")
    
    rand_res[k] = [["snp", "s288c", "greater", s1, p1],
                    ["snp", "w303", "greater", s2, p2],
                    ["pav", "s288c", "greater", s3, p3],
                    ["pav", "w303", "greater", s4, p4]]
    
    del snp_s288c_clust1, snp_s288c_clust2
    del pav_s288c_clust1, pav_s288c_clust2
    del s1, p1, s2, p2, s3, p3, s4, p4
    
    return rand_res[k]


# A wrapper function to handle multiprocessing and collect results
def multiprocessing_wrapper(k, env, rand_res, snp_shap_bench_s288c,
    snp_shap_bench_w303, snp_s288c_rclust1, snp_s288c_rclust2, snp_w303_rclust1,
    snp_w303_rclust2, pav_shap_bench_s288c, pav_shap_bench_w303,
    pav_s288c_rclust1, pav_s288c_rclust2, pav_w303_rclust1, pav_w303_rclust2):
    result = compared_rand_clusters(k, env, rand_res, snp_shap_bench_s288c, snp_shap_bench_w303,
        snp_s288c_rclust1, snp_s288c_rclust2, snp_w303_rclust1,
        snp_w303_rclust2, pav_shap_bench_s288c, pav_shap_bench_w303,
        pav_s288c_rclust1, pav_s288c_rclust2, pav_w303_rclust1, pav_w303_rclust2)
    rand_res[k] = result


def comp_clust_main(env, randomized=True):
    '''Compare the shap values between the two clusters containing the lab
    strains and the most distinct cluster to it using Mann-Whitney U (alternative
    hypothesis is "greater").
    
    If randomized=True, will generate 500 random cluster pairs to calculate
    statistical significance.
    
    If randomized=False, will calculate statistical significance for the two
    actual clusters (cluster 1 containing the lab strains; the other cluster
    being the most distinct to cluster 1).'''
    
    # read in shap values
    snp_shap = dt.fread(f"{d}/SNP/baseline/SHAP_values_sorted_{env}_snp_rf_baseline_training.txt").to_pandas()
    pav_shap = dt.fread(f"{d}/PAV/baseline/SHAP_values_sorted_{env}_pav_rf_baseline_training.txt").to_pandas()
    snp_shap.set_index("ID", inplace=True)
    pav_shap.set_index("ID", inplace=True)
    
    # map features to genes
    print("Mapping features to genes...")
    snp_shap.loc["Gene"] = snp_shap.T.index.map(map_snps.set_index("snp")["gene"]) # map features to genes
    pav_shap.columns = pav_shap.columns.map(lambda x: re.sub("^X", "", x))
    pav_shap.columns = pav_shap.columns.map(lambda x: re.sub("\.", "-", x))
    pav_shap.loc["Gene"] = pav_shap.T.index.map(map_orfs.set_index("orf")["gene"])
    
    # ensure isolate columns are the same order as the distance matrices
    snp_shap = snp_shap.loc[eu_dist_snp.columns.tolist()[:-1] + ["Gene"],:]
    pav_shap = pav_shap.loc[eu_dist_pav.columns.tolist()[:-1] + ["Gene"],:]
    
    # subset the benchmark genes
    print("Subsetting benchmark genes...")
    if env == "YPDBENOMYL500":
        snp_shap_bench = snp_shap.loc[:, snp_shap.loc["Gene"].isin(ben_snp)] # all benchmark genes validated in benomyl
        pav_shap_bench = pav_shap.loc[:, pav_shap.loc["Gene"].isin(ben_orf)]
        snp_shap_bench_s288c = snp_shap_bench.loc[:, snp_shap_bench.loc["Gene"].isin(ben_s288c)] # only benchmark genes validated in s288c and benomyl
        snp_shap_bench_w303 = snp_shap_bench.loc[:, snp_shap_bench.loc["Gene"].isin(ben_w303)] # only benchmark genes validated in w303 and benomyl
        pav_shap_bench_s288c = pav_shap_bench.loc[:, pav_shap_bench.loc["Gene"].isin(ben_s288c)]
        pav_shap_bench_w303 = pav_shap_bench.loc[:, pav_shap_bench.loc["Gene"].isin(ben_w303)]
        
    if (env == "YPDCAFEIN40") or (env == "YPDCAFEIN50"):
        snp_shap_bench = snp_shap.loc[:, snp_shap.loc["Gene"].isin(caf_snp)]
        pav_shap_bench = pav_shap.loc[:, pav_shap.loc["Gene"].isin(caf_orf)]
        snp_shap_bench_s288c = snp_shap_bench.loc[:, snp_shap_bench.loc["Gene"].isin(caf_s288c)]
        snp_shap_bench_w303 = snp_shap_bench.loc[:, snp_shap_bench.loc["Gene"].isin(caf_w303)]
        pav_shap_bench_s288c = pav_shap_bench.loc[:, pav_shap_bench.loc["Gene"].isin(caf_s288c)]
        pav_shap_bench_w303 = pav_shap_bench.loc[:, pav_shap_bench.loc["Gene"].isin(caf_w303)]
        
    if env == "YPDCUSO410MM":
        snp_shap_bench = snp_shap.loc[:, snp_shap.loc["Gene"].isin(cu_snp)]
        pav_shap_bench = pav_shap.loc[:, pav_shap.loc["Gene"].isin(cu_orf)]
        snp_shap_bench_s288c = snp_shap_bench.loc[:, snp_shap_bench.loc["Gene"].isin(cu_s288c)]
        pav_shap_bench_s288c = pav_shap_bench.loc[:, pav_shap_bench.loc["Gene"].isin(cu_s288c)]
        # I did not include w303 bc there is only one benchmark gene validated in w303 for CuSO4
        snp_shap_bench_w303 = pd.Series()
        pav_shap_bench_w303 = pd.Series()
        
    if env == "YPDSODIUMMETAARSENITE":
        snp_shap_bench = snp_shap.loc[:, snp_shap.loc["Gene"].isin(sma_snp)]
        pav_shap_bench = pav_shap.loc[:, pav_shap.loc["Gene"].isin(sma_orf)]
        snp_shap_bench_s288c = snp_shap_bench.loc[:, snp_shap_bench.loc["Gene"].isin(sma_s288c)]
        snp_shap_bench_w303 = snp_shap_bench.loc[:, snp_shap_bench.loc["Gene"].isin(sma_w303)]
        pav_shap_bench_s288c = pav_shap_bench.loc[:, pav_shap_bench.loc["Gene"].isin(sma_s288c)]
        pav_shap_bench_w303 = pav_shap_bench.loc[:, pav_shap_bench.loc["Gene"].isin(sma_w303)]
    
    # subset the individuals in the two clusters
    snp_shap_bench.loc[:,"Cluster"] = cluster_assignments_snp.loc[snp_shap_bench.index[:-1], "Cluster"]
    pav_shap_bench.loc[:,"Cluster"] = cluster_assignments_pav.loc[pav_shap_bench.index[:-1], "Cluster"]
    snp_s288c_clust1 = snp_shap_bench_s288c.loc[snp_shap_bench.Cluster == s288c_distinct_clusters_snp[0],:]
    snp_s288c_clust2 = snp_shap_bench_s288c.loc[snp_shap_bench.Cluster == s288c_distinct_clusters_snp[1],:]
    pav_s288c_clust1 = pav_shap_bench_s288c.loc[pav_shap_bench.Cluster == s288c_distinct_clusters_pav[0],:]
    pav_s288c_clust2 = pav_shap_bench_s288c.loc[pav_shap_bench.Cluster == s288c_distinct_clusters_pav[1],:]
    
    if env != "YPDCUSO410MM": # only 1 CuSO4 gene was verified in W303, can't calculate stats on this
        snp_w303_clust1 = snp_shap_bench_w303.loc[snp_shap_bench.Cluster == w303_distinct_clusters_snp[0],:]
        snp_w303_clust2 = snp_shap_bench_w303.loc[snp_shap_bench.Cluster == w303_distinct_clusters_snp[1],:]
        pav_w303_clust1 = pav_shap_bench_w303.loc[pav_shap_bench.Cluster == w303_distinct_clusters_pav[0],:]
        pav_w303_clust2 = pav_shap_bench_w303.loc[pav_shap_bench.Cluster == w303_distinct_clusters_pav[1],:]
    else: # necessary since it's needed in the partial function
        snp_w303_clust1 = pd.Series()
        snp_w303_clust2 = pd.Series()
        pav_w303_clust1 = pd.Series()
        pav_w303_clust2 = pd.Series()
    
    if randomized == True:
        # randomly subset individuals to create two clusters
        print("Creating random clusters to test S288C or W303 benchmark genes...")
        snp_s288c_indices = [np.random.choice(snp_shap_bench_s288c.index[:-1], len(snp_s288c_clust1), replace=False) for _ in range(5000)]
        snp_s288c_indices2 = [np.random.choice(snp_shap_bench_s288c.index[:-1], len(snp_s288c_clust2), replace=False) for _ in range(5000)]
        pav_s288c_indices = [np.random.choice(pav_shap_bench_s288c.index[:-1], len(pav_s288c_clust1), replace=False) for _ in range(5000)]
        pav_s288c_indices2 = [np.random.choice(pav_shap_bench_s288c.index[:-1], len(pav_s288c_clust2), replace=False) for _ in range(5000)]
        snp_s288c_rclust1 = {k: snp_s288c_indices[k] for k in range(5000)}
        snp_s288c_rclust2 = {k: snp_s288c_indices2[k] for k in range(5000)}
        pav_s288c_rclust1 = {k: pav_s288c_indices[k] for k in range(5000)}
        pav_s288c_rclust2 = {k: pav_s288c_indices2[k] for k in range(5000)}
        
        if env != "YPDCUSO410MM":
            snp_w303_indices = [np.random.choice(snp_shap_bench_w303.index[:-1], len(snp_w303_clust1), replace=False) for _ in range(5000)]
            snp_w303_indices2 = [np.random.choice(snp_shap_bench_w303.index[:-1], len(snp_w303_clust2), replace=False) for _ in range(5000)]
            pav_w303_indices = [np.random.choice(pav_shap_bench_w303.index[:-1], len(pav_w303_clust1), replace=False) for _ in range(5000)]
            pav_w303_indices2 = [np.random.choice(pav_shap_bench_w303.index[:-1], len(pav_w303_clust2), replace=False) for _ in range(5000)]
            snp_w303_rclust1 = {k: snp_w303_indices[k] for k in range(5000)}
            snp_w303_rclust2 = {k: snp_w303_indices2[k] for k in range(5000)}
            pav_w303_rclust1 = {k: pav_w303_indices[k] for k in range(5000)}
            pav_w303_rclust2 = {k: pav_w303_indices2[k] for k in range(5000)}
        else: # necessary since it's passed to the partial function
            snp_w303_rclust1 = {k: pd.Series() for k in range(5000)}
            snp_w303_rclust2 = {k: pd.Series() for k in range(5000)}
            pav_w303_rclust1 = {k: pd.Series() for k in range(5000)}
            pav_w303_rclust2 = {k: pd.Series() for k in range(5000)}
        
        # Create a Manager dictionary to store the results
        k_values = list(range(5000))
        with mp.Manager() as manager:
            rand_res = manager.dict()  # shared dictionary
            
            # Create a partial function to pass rand_res
            func = partial(multiprocessing_wrapper, env=env, rand_res=rand_res,
                        snp_shap_bench_s288c=snp_shap_bench_s288c,
                        snp_shap_bench_w303=snp_shap_bench_w303,
                        snp_s288c_rclust1=snp_s288c_rclust1,
                        snp_s288c_rclust2=snp_s288c_rclust2,
                        snp_w303_rclust1=snp_w303_rclust1,
                        snp_w303_rclust2=snp_w303_rclust2,
                        pav_shap_bench_s288c=pav_shap_bench_s288c,
                        pav_shap_bench_w303=pav_shap_bench_w303,
                        pav_s288c_rclust1=pav_s288c_rclust1,
                        pav_s288c_rclust2=pav_s288c_rclust2,
                        pav_w303_rclust1=pav_w303_rclust1,
                        pav_w303_rclust2=pav_w303_rclust2)
            
            # Initialize a multiprocessing pool
            with mp.Pool(processes=mp.cpu_count()) as pool:
                # Map the k values to the function
                pool.map(func, k_values)
            
            # Convert the manager dictionary to a regular dictionary
            results = dict(rand_res)  # Extract the data while the manager is active
        
        return results
    
    if randomized == False:
        mwu_res = {"S288C": {"greater":{"shap":{}, "fitness":{}},
                     "two-sided":{"shap":{}, "fitness":{}}},
           "W303": {"greater":{"shap":{}, "fitness":{}},
                    "two-sided":{"shap":{}, "fitness":{}}}}
        boxplot_data = {"S288C": {"shap":{}, "fitness":{}},
                "W303": {"shap":{}, "fitness":{}}}
        
        # calculate the mann-whitney u test statistic on shap values of clusters
        print("Calculating Mann-Whitney U test statistic for S288C benchmark genes...")
        boxplot_data["S288C"]["shap"] = {"snp": {"clust1": snp_s288c_clust1,
                                                    "clust2": snp_s288c_clust2}}
        s, p = mannwhitneyu(np.array(snp_s288c_clust1).astype("float").flatten(),
                    np.array(snp_s288c_clust2).astype("float").flatten(),
                    alternative="greater")
        mwu_res["S288C"]["greater"]["shap"] = {"snp": {"statistic": s, "p-value": p}} # Ha: greater
        
        s, p = mannwhitneyu(np.array(snp_s288c_clust1).astype("float").flatten(),
                    np.array(snp_s288c_clust2).astype("float").flatten(),
                    alternative="two-sided")
        mwu_res["S288C"]["two-sided"]["shap"] = {"snp": {"statistic": s, "p-value": p}} # Ha: not equal
        
        boxplot_data["S288C"]["shap"]["pav"] = {"clust1": pav_s288c_clust1,
                                                    "clust2": pav_s288c_clust2}
        s, p = mannwhitneyu(np.array(pav_s288c_clust1).astype("float").flatten(),
                    np.array(pav_s288c_clust2).astype("float").flatten(),
                    alternative="greater")
        mwu_res["S288C"]["greater"]["shap"]["pav"] = {"statistic": s, "p-value": p} # Ha: greater
        
        s, p = mannwhitneyu(np.array(pav_s288c_clust1).astype("float").flatten(),
                    np.array(pav_s288c_clust2).astype("float").flatten(),
                    alternative="two-sided")
        mwu_res["S288C"]["two-sided"]["shap"]["pav"] = {"statistic": s, "p-value": p} # Ha: not equal
        
        if env != "YPDCUSO410MM":
            print("Calculating Mann-Whitney U test statistic for W303 benchmark genes...")
            boxplot_data["W303"]["shap"] = {"snp": {"clust1": snp_w303_clust1,
                                                        "clust2": snp_w303_clust2}}
            s, p = mannwhitneyu(np.array(snp_w303_clust1).astype("float").flatten(),
                        np.array(snp_w303_clust2).astype("float").flatten(),
                        alternative="greater")
            mwu_res["W303"]["greater"]["shap"] = {"snp": {"statistic": s, "p-value": p}}
            
            s, p = mannwhitneyu(np.array(snp_w303_clust1).astype("float").flatten(),
                        np.array(snp_w303_clust2).astype("float").flatten(),
                        alternative="two-sided")
            mwu_res["W303"]["two-sided"]["shap"] = {"snp": {"statistic": s, "p-value": p}}
            
            boxplot_data["W303"]["shap"]["pav"] = {"clust1": pav_w303_clust1,
                                                        "clust2": pav_w303_clust2}
            s, p = mannwhitneyu(np.array(pav_w303_clust1).astype("float").flatten(),
                        np.array(pav_w303_clust2).astype("float").flatten(),
                        alternative="greater")
            mwu_res["W303"]["greater"]["shap"]["pav"] = {"statistic": s, "p-value": p}
            
            s, p = mannwhitneyu(np.array(pav_w303_clust1).astype("float").flatten(),
                        np.array(pav_w303_clust2).astype("float").flatten(),
                        alternative="two-sided")
            mwu_res["W303"]["two-sided"]["shap"]["pav"] = {"statistic": s, "p-value": p}
        
        # calculate the mann-whitney u test statistic on fitness values of clusters
        print("Calculating Mann-Whitney U test statistic for fitness of S288C clusters...")
        boxplot_data["S288C"]["fitness"] = {"snp": {"clust1": pheno.loc[snp_s288c_clust1.index, env],
                                                "clust2": pheno.loc[snp_s288c_clust2.index, env]}}
        s, p = mannwhitneyu(pheno.loc[snp_s288c_clust1.index, env],
                        pheno.loc[snp_s288c_clust2.index, env], alternative="greater")
        mwu_res["S288C"]["greater"]["fitness"] = {"snp": {"statistic": s, "p-value": p}}
        
        s, p = mannwhitneyu(pheno.loc[snp_s288c_clust1.index, env],
                        pheno.loc[snp_s288c_clust2.index, env], alternative="two-sided")
        mwu_res["S288C"]["two-sided"]["fitness"] = {"snp": {"statistic": s, "p-value": p}}
        
        boxplot_data["S288C"]["fitness"]["pav"] = {"clust1": pheno.loc[pav_s288c_clust1.index, env],
                                                    "clust2": pheno.loc[pav_s288c_clust2.index, env]}
        s, p = mannwhitneyu(pheno.loc[pav_s288c_clust1.index, env],
                        pheno.loc[pav_s288c_clust2.index, env], alternative="greater")
        mwu_res["S288C"]["greater"]["fitness"]["pav"] = {"statistic": s, "p-value": p}
        
        s, p = mannwhitneyu(pheno.loc[pav_s288c_clust1.index, env],
                        pheno.loc[pav_s288c_clust2.index, env], alternative="two-sided")
        mwu_res["S288C"]["two-sided"]["fitness"]["pav"] = {"statistic": s, "p-value": p}
        
        if env != "YPDCUSO410MM":
            print("Calculating Mann-Whitney U test statistic for fitness of W303 clusters...")
            boxplot_data["W303"]["fitness"] = {"snp": {"clust1": pheno.loc[snp_w303_clust1.index, env],
                                                    "clust2": pheno.loc[snp_w303_clust2.index, env]}}
            s, p = mannwhitneyu(pheno.loc[snp_w303_clust1.index, env],
                        pheno.loc[snp_w303_clust2.index, env], alternative="greater")
            mwu_res["W303"]["greater"]["fitness"] = {"snp": {"statistic": s, "p-value": p}}
            
            s, p = mannwhitneyu(pheno.loc[snp_w303_clust1.index, env],
                        pheno.loc[snp_w303_clust2.index, env], alternative="two-sided")
            mwu_res["W303"]["two-sided"]["fitness"] = {"snp": {"statistic": s, "p-value": p}}
            
            boxplot_data["W303"]["fitness"]["pav"] = {"clust1": pheno.loc[pav_w303_clust1.index, env],
                                                        "clust2": pheno.loc[pav_w303_clust2.index, env]}
            s, p = mannwhitneyu(pheno.loc[pav_w303_clust1.index, env],
                        pheno.loc[pav_w303_clust2.index, env], alternative="greater")
            mwu_res["W303"]["greater"]["fitness"]["pav"] = {"statistic": s, "p-value": p}
            
            s, p = mannwhitneyu(pheno.loc[pav_w303_clust1.index, env],
                        pheno.loc[pav_w303_clust2.index, env], alternative="two-sided")
            mwu_res["W303"]["two-sided"]["fitness"]["pav"] = {"statistic": s, "p-value": p}
        
        del snp_shap, pav_shap, snp_shap_bench, pav_shap_bench
        del snp_shap_bench_s288c, pav_shap_bench_s288c
        del snp_s288c_clust1, snp_s288c_clust2, pav_s288c_clust1, pav_s288c_clust2
        
        # if env != "YPDCUSO410MM":
        del snp_w303_clust1, snp_w303_clust2, pav_w303_clust1, pav_w303_clust2
        del snp_shap_bench_w303, pav_shap_bench_w303
        
        return mwu_res, boxplot_data


target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP"

mwu_res = {}
boxplot_data = {}
for env in target_envs:
    rand_res = comp_clust_main(env, randomized=True)
    mwu_res[env], boxplot_data[env] = comp_clust_main(env, randomized=False)
    
    rand_res_df = pd.DataFrame.from_dict(rand_res, orient='index')
    rand_res_df = rand_res_df.applymap(lambda x: pd.Series(x)).stack().apply(pd.Series)
    rand_res_df.columns = ["Data", "Strain", "Alternative", "Statistic", "P-value"]
    rand_res_df.to_csv(f"Scripts/Data_Vis/Section_5/PCA_snp_or_pav_genetic_distance_randomized_mwu_results_{env}.csv")


mwu_res = pd.DataFrame.from_dict({(i, j, k, h, l): mwu_res[i][j][k][h][l]
                        for i in mwu_res.keys()
                        for j in mwu_res[i].keys()
                        for k in mwu_res[i][j].keys()
                        for h in mwu_res[i][j][k].keys()
                        for l in mwu_res[i][j][k][h].keys()}, orient='index')

mwu_res.sort_values(by="p-value", inplace=True)
mwu_res.loc[mwu_res.index.get_level_values(1)=='greater',:]
mwu_res.index.names = ["Strain", "Alternative", "Factor", "Environment",  "Data"]
mwu_res.to_csv("Scripts/Data_Vis/Section_5/PCA_snp_or_pav_genetic_distance_mwu_results.csv")


#### Determine how many of the optimized features are unknown ORFs or intergenic SNPs.
# Feature to gene maps
map_snps = map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                       sep="\t", header=None, names=["snp", "chr", "pos", "gene"], index_col=0)
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv",
                       sep="\t",index_col=0)

# Optimized feature sets for all single environment models
snp_fs = pd.read_csv("Scripts/Data_Vis/Section_4/RF_FS_imp_snp_rank_per.tsv", sep="\t", index_col=0)
pav_fs = pd.read_csv("Scripts/Data_Vis/Section_4/RF_FS_imp_pav_rank_per.tsv", sep="\t", index_col=0)
cnv_fs = pd.read_csv("Scripts/Data_Vis/Section_4/RF_FS_imp_cnv_rank_per.tsv", sep="\t", index_col=0)

# For each env, count how many of the optimized features are intergenic/unknown ORFs
res = {}
for env in snp_fs.columns[1:]:
    env_snp_fs = snp_fs[[env]].dropna()
    env_pav_fs = pav_fs[[env]].dropna()
    env_cnv_fs = cnv_fs[[env]].dropna()
    
    res[env] = {"SNP": [sum(map_snps.loc[env_snp_fs.index, "gene"] == "intergenic"), # number of intergenic snps
                    len(env_snp_fs)], # total number of snps
                "PAV": [sum(pd.merge(map_orfs["gene"], env_pav_fs, how="right", # number of unknown orfs
                    left_index=True, right_index=True).gene.isna()), 
                    len(env_pav_fs)], # total number of orfs
                "CNV": [sum(pd.merge(map_orfs["gene"], env_cnv_fs, how="right", # number of unknown orfs
                    left_index=True, right_index=True).gene.isna()),
                    len(env_cnv_fs)]} # total number of orfs

res = pd.concat(
    {env: pd.DataFrame(data, index=["Unknown", "Total"]).T for env, data in res.items()}, 
    axis=1
)
res.to_csv("Scripts/Data_Vis/Section_5/Number_unknown_or_intergenic_features_RF_FS_models.csv")

res = pd.read_csv("Scripts/Data_Vis/Section_5/Number_unknown_or_intergenic_features_RF_FS_models.csv", index_col=0, header=[0,1])
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]

res_targ = res.loc[:,res.columns.get_level_values(0).isin(target_envs)]
res_targ_percent = []
for env in target_envs:
    res_targ_percent.append(res_targ[(env, "Unknown")] / res_targ[(env, "Total")])

pd.concat(res_targ_percent, axis=1).T.describe()
#             SNP       PAV       CNV
# count  5.000000  5.000000  5.000000
# mean   0.367927  0.639062  0.550000
# std    0.016731  0.045353  0.177259
# min    0.353000  0.570312  0.250000
# 25%    0.355469  0.632812  0.546875
# 50%    0.366000  0.632812  0.593750
# 75%    0.370167  0.671875  0.671875
# max    0.395000  0.687500  0.687500


#### How many of the unknown ORFs are missing in S288C and W303 (SACE_GAV)?
pav = pd.read_csv("Data/Peter_2018/ORFs_pres_abs_with_S288C.csv", index_col=0)
pav.index = pav.index.str.replace("^X", "", regex=True)
pav.index = pav.index.str.replace("\.", "-", regex=True)

pav.in_S288C.value_counts() # overall, including all 7708 ORFs; 5561 in S288C
pav.SACE_GAV.value_counts() # overall, 6043 in W303

# Optimized feature sets for all single environment PAV models
pav_fs = pd.read_csv("Scripts/Data_Vis/Section_4/RF_FS_imp_pav_rank_per.tsv", sep="\t", index_col=0)

res = []
for env in pav_fs.columns:
    env_pav_fs = pav_fs[[env]].dropna()
    env_counts = pav.loc[env_pav_fs.index, "in_S288C"].value_counts()
    env_counts.name = env
    res.append(env_counts)

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]
res_target = pd.concat(res, axis=1)[target_envs]
np.mean(res_target.T[0] / (res_target.T[0] + res_target.T[1])) # 0.6765625 (67.7% ORFs missing in S288C on average among the target env RF FS PAV features
np.std(res_target.T[0] / (res_target.T[0] + res_target.T[1])) # 0.04324485157218139 (4.3%)

res_w303 = []
for env in pav_fs.columns:
    env_pav_fs = pav_fs[[env]].dropna()
    env_counts = pav.loc[env_pav_fs.index, "SACE_GAV"].value_counts()
    env_counts.name = env
    res_w303.append(env_counts)

res_w303_target = pd.concat(res_w303, axis=1)[target_envs]
np.mean(res_w303_target.T[False] / (res_w303_target.T[True] + res_w303_target.T[False])) # 0.534375 (53.4 %) missing in W303 on average
np.std(res_w303_target.T[False] / (res_w303_target.T[False] + res_w303_target.T[True])) # 0.07680128579652817 (7.7 %)

