#!/usr/bin/env python3
"""From 0_Manuscript_Tables.py code TABLE S8
Need to submit as slurm job for it to completely run."""

import os
import re
import tqdm
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

################################################################################
### TABLE S8
################################################################################
## Compare ranking distribution of literature genes to random chance
target_envs = {"YPDCAFEIN40":"Caffeine", "YPDCAFEIN50":"Caffeine",
               "YPDBENOMYL500":"Benomyl", "YPDCUSO410MM":"CuSO4",
               "YPDSODIUMMETAARSENITE":"Sodium (meta)arsenite"} # most predictive envs

# Remove literature genes that are not found within our datasets
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
                       header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)

num_genes_snps = map_snps['gene'].nunique() - 1  # don't count intergenic
num_genes_orfs = map_orfs['gene'].nunique()

lit_genes = pd.read_csv("Data/Peter_2018/literature_genes.txt", sep="\t") # target genes from the literature
lit_genes_snps = lit_genes.loc[lit_genes.Gene.isin(map_snps.gene),:].reset_index().iloc[:,1:]
lit_genes_orfs = lit_genes.loc[lit_genes.Gene.isin(map_orfs.gene),:].reset_index().iloc[:,1:]

# List to collect results
permut_res = [["Model Type", "DataType", "ImpType", "Env", "NumGenes", \
               "NumLitGenesFound", "AvgRankPerSampDist", "AvgRankPerLitGenes", \
               "Z-score", "p-value", "MedRankPerSampDist", "MedRankPerSampDist",
               "MedRank Z-score", "MedRank p-value"]] # collect results from permutation test

# Compare max rank percentile of literature genes identified by the model to random chance
n = 1000000 # number of repetitions
seeds = 3*2*5*n # random seed counter
for data_type in ["snp", "pav", "cnv"]:
    for imp_type in ["imp", "shap"]:
        # Read and prepare dataset
        df = dt.fread(f"Scripts/Data_Vis/RF_baseline_{imp_type}_{data_type}.tsv").to_pandas()
        if data_type == "snp":
            df.set_index(["snp", "gene"], inplace=True)
            df = df.iloc[:,1:] # remove "gene_with_intergenic" column; ranking will not includer intergenic snps
            df = df.loc[df.index.get_level_values("gene") != "intergenic",:] # drop intergenic snps
        else:
            df.set_index(["orf", "gene"], inplace=True)
        if imp_type == "shap":
            df = df.abs()
        ## Compare ranking distribution of literature genes to random chance
        for env in target_envs.keys():
            print(data_type, imp_type, env)
            # Rank the genes based on gini importance or shap value
            env_imp = df.loc[:,env].dropna().sort_values(ascending=False) # all the genic snps for this env, excluding intergenic snps
            max_env_ranks = env_imp.groupby("gene").max() # take the max importance (either gini or shap) per gene
            max_env_ranks.sort_values(ascending=False, inplace=True) # sort from highest importance to lowest
            max_env_ranks = pd.concat([max_env_ranks, max_env_ranks.rank(method="average", pct=True)], axis=1) # assign ordinal rank
            max_env_ranks.columns = ["importance", "rank_per"]
            max_lit_ranks = max_env_ranks.loc[max_env_ranks.index.get_level_values("gene")\
                                .isin(lit_genes.loc[lit_genes.Env==target_envs[env]].Gene),:]
            ## Begin random sampling
            if max_lit_ranks.shape[0] > 1: # needed for calculating std of t-test and z-test
                rand_ranks = [] # distribution of mean rank percentiles of the literature genes
                rand_ranks_med = [] # medians
                for i in tqdm.tqdm(range(n)):
                    rand_ranks.append(max_env_ranks.sample(n=max_lit_ranks.shape[0],
                                    random_state=seeds, ignore_index=True).rank_per.mean()) # average max rank percentile of random sample
                    rand_ranks_med.append(max_env_ranks.sample(n=max_lit_ranks.shape[0],
                                    random_state=seeds, ignore_index=True).rank_per.median()) # average max rank percentile of random sample
                    seeds -= 1
                # Figure to check for normality
                fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(9, 3))
                ax[0].hist(max_lit_ranks["rank_per"], alpha=0.5, label="original")
                ax[1].hist(rand_ranks, bins=10000, alpha=0.5, label="samp dist of means")
                ax[1].axvlines(x=np.mean(rand_ranks), color="black", linestyles="dashed", label="mean of samp dist of means")
                ax[1].axvlines(x=max_lit_ranks.rank_per.mean(), color="red", linestyles="dashed", label="mean of lit genes")
                ax[2].hist(rand_ranks_med, bins=10000, alpha=0.5, label="samp dist of medians")
                ax[2].axvlines(x=np.median(rand_ranks_med), color="black", linestyles="dashed", label="med of samp dist of meds")
                ax[2].axvlines(x=max_lit_ranks.rank_per.median(), color="red", linestyle="dashed", label="med of lit genes")
                ax[0].legend(loc="upper right")
                ax[1].legend(loc="upper right")
                ax[2].legend(loc="upper right")
                ax[0].set_xlabel("Rank percentile")
                ax[1].set_xlabel("Mean rank percentile")
                ax[2].set_xlabel("Median rank percentile")
                ax[0].set_ylabel("Frequency")
                ax[1].set_ylabel("Frequency")
                ax[2].set_ylabel("Frequency")
                plt.savefig(f"Scripts/Data_Vis/Samp_dist_rank_per_means_{data_type}_{imp_type}_{env}.pdf")
                plt.close()
                ## Right-tailed Z-score
                # note: if z_score >= z_critical, reject the null hypothesis
                z_score = (max_lit_ranks.rank_per.mean() - np.mean(rand_ranks)) / \
                          np.std(rand_ranks)
                z_critical = stats.norm.ppf(1-0.05) # for alpha = 0.05, the z_critical = 1.6448536269514722
                z_p_val = 1 - stats.norm.cdf(np.abs(z_score))
                ## Right-tailed modified Z-score
                deviation_from_med = np.array(max_lit_ranks.rank_per) - np.median(rand_ranks)
                mad = np.median(np.abs(deviation_from_med))
                mod_z_score = deviation_from_med/(1.4826*mad) # the consistency correction k for the normal distribution is 1.4826
                mod_z_p_val = 1 - stats.norm.cdf(np.abs(mod_z_score))
                # Save results
                permut_res.append(["RF baseline", data_type, imp_type, env, \
                    len(max_env_ranks), len(max_lit_ranks), np.mean(rand_ranks), \
                    max_lit_ranks.rank_per.mean(), z_score, z_p_val, \
                    np.median(rand_ranks_med), max_lit_ranks.rank_per.median(), \
                    mod_z_score, mod_z_p_val])
            else:
                permut_res.append(["RF baseline", data_type, imp_type, env, \
                    len(max_env_ranks), len(max_lit_ranks), len(max_lit_ranks), \
                    np.mean(rand_ranks), None, None, np.median(rand_ranks_med), \
                    max_lit_ranks.rank_per.median(), None, None])

pd.DataFrame(permut_res).to_csv("Scripts/Data_Vis/Table_S8_mean_lit_gene_rank_per_vs_random_v2.tsv",
    sep="\t", header=False, index=False)
