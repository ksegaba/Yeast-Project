################################################################################
### Figure 4
################################################################################

import os, glob, re
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")


## Map the optimized RF model feature importance files to the benchmark genes and draw a heatmap ##
# Read in the feature to gene maps with benchmark gene information
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")

d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project"

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]

## Make a heatmap including only the benchmark genes within the optimized features
combined_opt = {"gini": {}, "shap":{}}

for data_type in ["snp", "pav", "cnv"]:
	for imp_type in ["imp", "shap"]:
		imp = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}.tsv", sep="\t", index_col=0)
		imp = imp[["gene"] + target_envs] # keep only the target envs
		imp = imp.loc[~imp.index.isna()]
		
		imp = imp.loc[imp.select_dtypes("number").dropna(how="all").index] # drop irrelevant rows
		
		# keep only the benchmark genes
		if data_type == "snp":
			map_snps_sub = map_snps.set_index("snp").loc[imp.index]
			map_snps_sub = map_snps_sub.select_dtypes("number").\
				loc[map_snps_sub.select_dtypes("number").\
					apply(lambda row: (row==0).sum() != 4, axis=1)] # keep only benchmark genes
			
			imp = imp.loc[map_snps_sub.index]
			print(sum(imp.gene == "intergenic")) # sanity check, should be 0
			
			imp = imp.groupby("gene").max() # max feature importance per gene
		
		else:
			map_orfs_sub = map_orfs.loc[map_orfs.orf.isin(imp.index)]
			map_orfs_sub.set_index("orf", inplace=True)
			map_orfs_sub = map_orfs_sub.select_dtypes("number").\
				loc[map_orfs_sub.select_dtypes("number").\
					apply(lambda row: (row==0).sum() != 4, axis=1)]
			
			imp = imp.loc[map_orfs_sub.index]
			imp = imp.groupby("gene").max()
		
		if imp_type == "imp":
			combined_opt["gini"][data_type] = imp
		else:
			combined_opt["shap"][data_type] = imp
		
		del imp


def combine_2level_dict_dfs(combined_opt):
	out_df = {}
	for outer_key, inner_dict in combined_opt.items():
		out_df[outer_key] = pd.concat(inner_dict.values(), axis=1, keys=inner_dict.keys())
	return pd.concat(out_df.values(), axis=1, keys=out_df.keys())


combined_opt_df = combine_2level_dict_dfs(combined_opt)
combined_opt_df.to_csv("Scripts/Data_Vis/Section_4/Figure_4_benchmark_genes_in_opt_models_data.tsv", sep="\t")

# Draw the heatmap
combined_opt_df.fillna(0, inplace=True) # set NaN values to 0 importance
combined_opt_df_rank = combined_opt_df.rank(method="average", ascending=False) # assign importance scores to rankings
# smallest imp scores will be further down the ranking. The most important gene will have a rank of 1.

sns.clustermap(combined_opt_df_rank, method="average", metric="euclidean",
			   col_cluster=False, yticklabels=False, figsize=(8.5, 11),
			   cbar_kws={"orientation":"horizontal"})
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Figure_4_benchmark_genes_in_opt_models.pdf")
plt.close()

## Make a heatmap for the top 20 genes
# Feature to gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
					   sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")

combined_top_20 = {"gini": {}, "shap":{}} # top 20 gene importances

for data_type in ["snp", "pav", "cnv"]:
	for imp_type in ["imp", "shap"]:
		rank_per = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}_rank_per.tsv", sep="\t", index_col=0)
		# in rank_per: 1.00 is the most important feature. Importance decreases as rank_per reaches 0.
		
		# map the features to genes
		if data_type == "snp":
			map_snps_sub = map_snps.set_index("snp").loc[rank_per.index]
			rank_per.insert(0, "gene", map_snps_sub.gene)
			rank_per = rank_per.groupby("gene").max()
			rank_per.drop(index="intergenic", inplace=True)
		else:
			rank_per["gene"] = rank_per.apply(lambda x: x.name \
				if x.name not in map_orfs.orf.values \
				else map_orfs.loc[map_orfs.orf==x.name, "gene"].values[0], axis=1)
			rank_per = rank_per.groupby("gene").max()
		
		# keep only the top 20 features in the target environments
		top_20 = []
		for env in target_envs:
			top_20.append(rank_per.sort_values(by=env, ascending=False)[env].iloc[:20])
			# tmp = rank_per.loc[rank_per.index != "intergenic"].\
			# 	sort_values(by=env, ascending=False)[env]
		
		top_20 = pd.concat(top_20, axis=1)
		
		if imp_type == "imp":
			combined_top_20["gini"][data_type] = top_20
		else:
			combined_top_20[imp_type][data_type] = top_20

# Combine innermost dataframe values in a nested dictionary with 2 levels
combined_top_20_df = combine_2level_dict_dfs(combined_top_20)

# map the features to genes and benchmark genes
# genes = []
# for feature in combined_top_20_df.index:
#     if feature in map_snps["snp"].values:
#         genes.append(map_snps.loc[map_snps["snp"] == feature, "gene"].values[0])
#     elif feature in map_orfs["orf"].values:
#         genes.append(map_orfs.loc[map_orfs["orf"] == feature, "gene"].values[0])
#     else:
#         genes.append(feature)

# combined_top_20_df.insert(0, "genes", genes)
# combined_top_20_df.reset_index(inplace=True)
# combined_top_20_df = combined_top_20_df.drop(columns="index").set_index("genes")

curated = pd.read_csv("Data/SGD_Experiment_Genes/manually_curated_genes.txt", sep="\t")
benomyl = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caffeine = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt").to_pandas()
cu = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")

combined_top_20_df["Benomyl"] = combined_top_20_df.index.isin(benomyl["Gene Systematic Name"].values).astype(int)
combined_top_20_df["Caffeine"] = combined_top_20_df.index.isin(caffeine["Gene Systematic Name"].values).astype(int)
combined_top_20_df["CuSO4"] = combined_top_20_df.index.isin(cu["Gene Systematic Name"].values).astype(int)
combined_top_20_df["Sodium_meta-arsenite"] = combined_top_20_df.index.isin(sma["Gene Systematic Name"].values).astype(int)
combined_top_20_df["Curated"] = combined_top_20_df.index.isin(curated["gene"].values).astype(int)

combined_top_20_df[['Benomyl', 'Caffeine', 'CuSO4', 'Sodium_meta-arsenite', 'Curated']].sum()
# Benomyl                      9
# Caffeine                    19
# CuSO4                        7
# Sodium_meta-arsenite         9
# Curated                      4
# dtype: int64

combined_top_20_df.shape # (238, 35)
combined_top_20_df.insert(30, "DIVIDER", 0)

# Heatmaps of the top 20 features with a discrete color bar
bounds = [i for i in range(1,21)]
for env in target_envs:
	env_df = combined_top_20_df.loc[:, combined_top_20_df.columns.get_level_values(2) == env].dropna(how="all")
	env_df.sort_values(by=("gini", "snp", env), ascending=False, inplace=True)
	env_df = env_df.rank(ascending=False)
	
	env_df["Benomyl"] = env_df.index.isin(benomyl["Gene Systematic Name"].values).astype(int)
	env_df["Caffeine"] = env_df.index.isin(caffeine["Gene Systematic Name"].values).astype(int)
	env_df["CuSO4"] = env_df.index.isin(cu["Gene Systematic Name"].values).astype(int)
	env_df["Sodium_meta-arsenite"] = env_df.index.isin(sma["Gene Systematic Name"].values).astype(int)
	env_df["Curated"] = env_df.index.isin(curated["gene"].values).astype(int)
	
	fig, ax = plt.subplots(1, 2, figsize=(8.5, 15))
	sns.heatmap(env_df.iloc[:,:6], cmap="viridis", cbar_kws={"label": "Rank",
		"ticks":bounds}, yticklabels=True, ax=ax[0], vmin=0, vmax=20)
	sns.heatmap(env_df.iloc[:,6:], cmap="Reds", cbar_kws={"label": "Benchmark gene"},
		yticklabels=False, ax=ax[1])
	plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_4_top_20_features_{env}_rank_percentiles.pdf")
	plt.close()
	
	# Save env_df
	env_df.to_csv(f"Scripts/Data_Vis/Section_4/Figure_4_top_20_features_{env}_rank_percentiles.tsv", sep="\t")
	del env_df
