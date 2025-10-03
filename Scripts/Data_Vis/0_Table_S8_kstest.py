#!/usr/bin/python3
import os
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import ks_2samp

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

################################################################################
### TABLE S8
################################################################################
## Compare ranking distribution of literature genes to random chance
target_envs = {"YPDCAFEIN40":"Caffeine", "YPDCAFEIN50":"Caffeine",
               "YPDBENOMYL500":"Benomyl", "YPDCUSO410MM":"CuSO4",
               "YPDSODIUMMETAARSENITE":"Sodium (meta)arsenite"} # most predictive envs

# Genes lists from SGD
ben_genes = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf_genes = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt").to_pandas()
cu_genes = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma_genes = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")

# List to collect results
permut_res = [["Model Type", "DataType", "ImpType", "Env", "NumGenes",
               "NumGenesNonzero", "NumLitGenesFound", "NumLitGenesFoundNonzero",
               "AvgRankPerSampDist", "AvgRankPerLitGenesNonzero",
               "D", "p-value", "MedRankPerSampDist", "MedRankPerSampDist",
               "MedRank D", "MedRank p-value"]] # collect results from permutation test

# Compare max rank percentile of literature genes identified by the model to random chance
n = 1000000 # number of repetitions
seeds = 3*2*5*n # random seed counter
for data_type in ["snp", "pav", "cnv"]:
    for imp_type in ["imp", "shap"]:
        # Read and prepare dataset
        df = dt.fread(f"Scripts/Data_Vis/RF_baseline_{imp_type}_{data_type}.tsv").to_pandas()
        
        if data_type == "snp":
            df.set_index(["snp", "gene"], inplace=True)
            df = df.iloc[:,1:] # remove "gene_with_intergenic" column; ranking will not include intergenic snps
            df = df.loc[df.index.get_level_values("gene") != "intergenic",:] # drop intergenic snps
        else:
            df.gene = df.apply(lambda x: x.orf if x.gene=="" else x.gene, axis=1) # orfs with no gene matches
            # df = df.loc[df.gene!="",:] # drop orfs with no gene matches???
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
            max_env_ranks_nonzero = max_env_ranks.loc[max_env_ranks > 0.0] # remove values with gini importance or shap value of zero
            max_env_ranks_nonzero = pd.concat([max_env_ranks_nonzero, 
                max_env_ranks_nonzero.rank(method="average", pct=True)], axis=1) # assign rank percentile
            max_env_ranks_nonzero.columns = ["importance", "rank_per"]
            max_env_ranks_nonzero.to_csv(f"Scripts/Data_Vis/RF_baseline_{imp_type}_{data_type}_{env}_max_gene_rank_per.tsv", sep="\t")
            
            # Get the ranks of the literature genes
            if env == "YPDBENOMYL500":
                max_lit_ranks = max_env_ranks.loc[max_env_ranks.index.isin(ben_genes["Gene Systematic Name"])]
                max_lit_ranks_nonzero = max_env_ranks_nonzero.loc[max_env_ranks_nonzero.index.isin(ben_genes["Gene Systematic Name"]),:]
            if env == "YPDCAFEIN40":
                max_lit_ranks = max_env_ranks.loc[max_env_ranks.index.isin(caf_genes["Gene Systematic Name"])]
                max_lit_ranks_nonzero = max_env_ranks_nonzero.loc[max_env_ranks_nonzero.index.isin(caf_genes["Gene Systematic Name"]),:]
            if env == "YPDCAFEIN50":
                max_lit_ranks = max_env_ranks.loc[max_env_ranks.index.isin(caf_genes["Gene Systematic Name"])]
                max_lit_ranks_nonzero = max_env_ranks_nonzero.loc[max_env_ranks_nonzero.index.isin(caf_genes["Gene Systematic Name"]),:]
            if env == "YPDCUSO410MM":
                max_lit_ranks = max_env_ranks.loc[max_env_ranks.index.isin(cu_genes["Gene Systematic Name"])]
                max_lit_ranks_nonzero = max_env_ranks_nonzero.loc[max_env_ranks_nonzero.index.isin(cu_genes["Gene Systematic Name"]),:]
            if env == "YPDSODIUMMETAARSENITE":
                max_lit_ranks = max_env_ranks.loc[max_env_ranks.index.isin(sma_genes["Gene Systematic Name"])]
                max_lit_ranks_nonzero = max_env_ranks_nonzero.loc[max_env_ranks_nonzero.index.isin(cu_genes["Gene Systematic Name"]),:]
            max_lit_ranks_nonzero.to_csv(f"Scripts/Data_Vis/RF_baseline_{imp_type}_{data_type}_{env}_max_lit_gene_rank_per.tsv", sep="\t")

            ## Begin random sampling
            rand_ranks = [] # sampling distribution of mean rank percentiles of the literature genes
            rand_ranks_med = [] # medians
            for i in tqdm(range(n)):
                rand_ranks.append(max_env_ranks_nonzero.sample(n=max_lit_ranks_nonzero.shape[0],
                                random_state=seeds, ignore_index=True).rank_per.mean()) # average max rank percentile of random sample
                rand_ranks_med.append(max_env_ranks_nonzero.sample(n=max_lit_ranks_nonzero.shape[0],
                                random_state=seeds, ignore_index=True).rank_per.median()) # average max rank percentile of random sample
                seeds -= 1
            
            # Figure to check for normality
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(9, 3))
            ax[0].hist(max_lit_ranks_nonzero["rank_per"], alpha=0.5, label="original")
            ax[1].hist(rand_ranks, bins=10000, alpha=0.5, label="samp dist of means")
            ax[1].axvline(x=np.mean(rand_ranks), color="black", label="mean of samp dist of means")
            ax[1].axvline(x=max_lit_ranks_nonzero.rank_per.mean(), color="red", label="mean of lit genes")
            ax[2].hist(rand_ranks_med, bins=10000, alpha=0.5, label="samp dist of medians")
            ax[2].axvline(x=np.median(rand_ranks_med), color="black", label="med of samp dist of meds")
            ax[2].axvline(x=max_lit_ranks_nonzero.rank_per.median(), color="red", label="med of lit genes")
            ax[0].legend(loc="upper right")
            ax[1].legend(loc="upper right")
            ax[2].legend(loc="upper right")
            ax[0].set_xlabel("Rank percentile")
            ax[1].set_xlabel("Mean rank percentile")
            ax[2].set_xlabel("Median rank percentile")
            ax[0].set_ylabel("Frequency")
            ax[1].set_ylabel("Frequency")
            ax[2].set_ylabel("Frequency")
            plt.savefig(f"Scripts/Data_Vis/Samp_dist_rank_per_means_and_medians_{data_type}_{imp_type}_{env}.pdf")
            plt.close()
            
            ## KS-test
            D_mean, p_mean = ks_2samp(max_lit_ranks_nonzero.rank_per, rand_ranks, alternative="greater")
            D_med, p_med = ks_2samp(max_lit_ranks_nonzero.rank_per, rand_ranks_med, alternative="greater")
            
            # Save results
            permut_res.append(["RF baseline", data_type, imp_type, env, \
                len(max_env_ranks), len(max_env_ranks_nonzero), len(max_lit_ranks),\
                len(max_lit_ranks_nonzero), np.mean(rand_ranks), \
                max_lit_ranks_nonzero.rank_per.mean(), D_mean, p_mean, \
                np.median(rand_ranks_med), max_lit_ranks_nonzero.rank_per.median(), \
                D_med, p_med])

pd.DataFrame(permut_res).to_csv("Scripts/Data_Vis/Table_S8_mean_lit_gene_rank_per_vs_random_v3.tsv",
    sep="\t", header=False, index=False)