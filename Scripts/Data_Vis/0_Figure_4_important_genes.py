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

# Feature to gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
					   sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_16_removed.txt", sep="\t")

## Map the optimized RF model feature importance files to the benchmark genes and draw a heatmap ##

# Read in the feature to gene maps with benchmark gene information
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_16_removed_expanded_benchmark.tsv", sep="\t")

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
					apply(lambda row: (row==0).sum() != 4, axis=1)]
			
			imp = imp.loc[map_snps_sub.index]
			print(sum(imp.gene == "intergenic")) # sanity check, should be 0
			
			imp = imp.groupby("gene").max()
		
		else:
			# imp.merge(map_orfs.set_index("orf")[["Benomyl", "Caffeine",
			# 	"CuSO4", "Sodium_meta-arsenite"]], left_index=True,
			# 	right_index=True, how="left") # there is an issue with map_orfs, so I won't drop the "unknown ORFs", I'll have to manually see where in the blast results they're at.
			# For now, let's keep all the ORFs
			imp.gene = imp.apply(lambda row: row["gene"] if pd.notna(row["gene"]) else row.name, axis=1)
			imp = imp.reset_index().drop(columns="orf")
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
combined_opt_df.to_csv("Scripts/Data_Vis/Section_4/Figure_4_data_NEED_TO_VALIDATE_ORF_GENES.tsv", sep="\t")
# March 14, 2025: get a list of ORFs that should be in map_orfs.
# Consider writing a script to get the weird ncbi gene IDs for all the orfs in pav/cnv datasets from the blastx and tblastx results.
# Then, generate a systematic gene ID to weird nbi gene IDs map.
# Re-generate map_orfs.
# This is a lot, should I leave this to the review phase?? I really need to give SH my manuscript.


## Make a heatmap for the top 20 genes
combined_top_20 = {"gini": {}, "shap":{}}
combined_top_per75 = {"gini": {}, "shap":{}}

for data_type in ["snp", "pav", "cnv"]:
	for imp_type in ["imp", "shap"]:
		adjusted_rank_per = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}_adjusted_rank_per.tsv", sep="\t", index_col=0)
		
		# keep only the rank percentile columns
		adjusted_rank_per = adjusted_rank_per[[c for c in adjusted_rank_per.columns if "rank_per" in c]]
		adjusted_rank_per.columns = [re.sub("_rank_per", "", c) for c in adjusted_rank_per.columns]
		
		# keep only the top 20 features in the target environments
		top_20 = []
		top_per75 = []
		for env in target_envs:
			top_20.append(adjusted_rank_per.loc[
				adjusted_rank_per.index != "intergenic"].\
				sort_values(by=env, ascending=False)[env].iloc[:20])
			
			tmp = adjusted_rank_per.loc[adjusted_rank_per.index != "intergenic"].\
				sort_values(by=env, ascending=False)[env]
			top_per75.append(tmp[tmp >= 0.75])
		
		top_20 = pd.concat(top_20, axis=1)
		top_per75 = pd.concat(top_per75, axis=1)
		
		if imp_type == "imp":
			combined_top_20["gini"][data_type] = top_20
			combined_top_per75["gini"][data_type] = top_per75
		else:
			combined_top_20[imp_type][data_type] = top_20
			combined_top_per75[imp_type][data_type] = top_per75

# Combine innermost dataframe values in a nested dictionary with 2 levels
def combine_2level_dict_dfs(combined_top_20):
	out_df = {}
	for outer_key, inner_dict in combined_top_20.items():
		out_df[outer_key] = pd.concat(inner_dict.values(), axis=1, keys=inner_dict.keys())
	return pd.concat(out_df.values(), axis=1, keys=out_df.keys())


combined_top_20_df = combine_2level_dict_dfs(combined_top_20)
combined_top_per75_df = combine_2level_dict_dfs(combined_top_per75)
combined_top_per75_df.shape # (6174, 30) too big to make heatmap perhaps?

# map the features to genes
curated = pd.read_csv("Data/SGD_Experiment_Genes/manually_curated_genes.txt", sep="\t")
combined_top_20_df["Benomyl"] = combined_top_20_df.index.isin(benomyl["Gene Systematic Name"].values).astype(int)
combined_top_20_df["Caffeine"] = combined_top_20_df.index.isin(caffeine["Gene Systematic Name"].values).astype(int)
combined_top_20_df["CuSO4"] = combined_top_20_df.index.isin(cuso4["Gene Systematic Name"].values).astype(int)
combined_top_20_df["Sodium_meta-arsenite"] = combined_top_20_df.index.isin(sma["Gene Systematic Name"].values).astype(int)
combined_top_20_df["Curated"] = combined_top_20_df.index.isin(curated["gene"].values).astype(int)

combined_top_20_df[['Benomyl', 'Caffeine', 'CuSO4', 'Sodium_meta-arsenite', 'Curated']].sum()
# Benomyl                     10
# Caffeine                    21
# CuSO4                        7
# Sodium_meta-arsenite        11
# Curated                      3
# dtype: int64

combined_top_20_df.shape # (237, 35)
combined_top_20_df.insert(30, "DIVIDER", 0)

# Heatmaps of the top 20 features with a discrete color bar
for env in target_envs:
	env_df = combined_top_20_df.loc[:, combined_top_20_df.columns.get_level_values(2) == env].dropna(how="all")
	env_df.sort_values(by=("gini", "snp", env), ascending=False, inplace=True)
	env_df = env_df.rank()
	
	env_df["Benomyl"] = env_df.index.isin(benomyl["Gene Systematic Name"].values).astype(int)
	env_df["Caffeine"] = env_df.index.isin(caffeine["Gene Systematic Name"].values).astype(int)
	env_df["CuSO4"] = env_df.index.isin(cuso4["Gene Systematic Name"].values).astype(int)
	env_df["Sodium_meta-arsenite"] = env_df.index.isin(sma["Gene Systematic Name"].values).astype(int)
	env_df["Curated"] = env_df.index.isin(curated["gene"].values).astype(int)
	
	fig, ax = plt.subplots(1, 2, figsize=(8.5, 15))
	sns.heatmap(env_df.iloc[:,:6], cmap="viridis", cbar_kws={"label": "Rank",
		"ticks": bounds}, yticklabels=True, ax=ax[0])
	sns.heatmap(env_df.iloc[:,6:], cmap="Reds", cbar_kws={"label": "Benchmark gene"},
		yticklabels=False, ax=ax[1])
	plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_4_top_20_features_{env}_rank_percentiles.pdf")
	plt.close()
	
	# Save env_df
	env_df.to_csv(f"Scripts/Data_Vis/Section_4/Figure_4_top_20_features_{env}_rank_percentiles.tsv", sep="\t")
	del env_df

# Heatmaps of the 75th percentile most important features with a discrete color bar
for env in target_envs:
	env_df = combined_top_per75_df.loc[:, combined_top_per75_df.columns.get_level_values(2) == env].dropna(how="all")
	env_df.sort_values(by=("gini", "snp", env), ascending=False, inplace=True)
	
	env_df["Benomyl"] = env_df.index.isin(benomyl["Gene Systematic Name"].values).astype(int)
	env_df["Caffeine"] = env_df.index.isin(caffeine["Gene Systematic Name"].values).astype(int)
	env_df["CuSO4"] = env_df.index.isin(cuso4["Gene Systematic Name"].values).astype(int)
	env_df["Sodium_meta-arsenite"] = env_df.index.isin(sma["Gene Systematic Name"].values).astype(int)
	env_df["Curated"] = env_df.index.isin(curated["gene"].values).astype(int)
	print(env_df.shape)
	
	# Discrete colorbar
	bounds = [0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 1.0]
	colors = sns.color_palette("viridis", len(bounds) - 1)
	cmap = ListedColormap(colors)
	norm = BoundaryNorm(bounds, cmap.N)
	
	fig, ax = plt.subplots(1, 2, figsize=(8.5, 15))
	sns.heatmap(env_df.iloc[:,:6], cmap=cmap, norm=norm, cbar_kws={"label": "Rank percentile",
		"ticks": bounds}, yticklabels=False, ax=ax[0])
	sns.heatmap(env_df.iloc[:,6:], cmap="Reds", cbar_kws={"label": "Benchmark gene"},
		yticklabels=False, ax=ax[1])
	plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_4_top_per75_features_{env}_rank_percentiles.pdf")
	plt.close()
	
	# Save env_df
	env_df.to_csv(f"Scripts/Data_Vis/Section_4/Figure_4_top_per75_features_{env}_rank_percentiles.tsv", sep="\t")
	del env_df
