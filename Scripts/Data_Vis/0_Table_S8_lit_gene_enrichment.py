################################################################################
### TABLE S
################################################################################

import os
import datatable as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

## ENRICHMENT OF BENCHMARK STRESS-RESPONSE GENES IN RF FS MODELS
# Target phenotype and mutant information annotations
phenotypes = ["resistance to chemicals: decreased ", "viability: decreased ",
			  "resistance to chemicals: decreased", "metal resistance: decreased",
			  "metal resistance: decreased ", "oxidative stress resistance: decreased",
			  "respiratory growth: decreased", "stress resistance: decreased"] # sensitive mutants
mutants = ["null Allele", "reduction of function", "reduction of function Allele"]

## Genes lists from SGD
benomyl_sgd = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations.txt", skiprows=8, sep="\t")
benomyl_sgd[["Mutant Type", "Mutant Description"]] = benomyl_sgd["Mutant Information"].str.split(":", expand=True)
# benomyl_sgd.groupby("Mutant Type").count() # see mutant types
ben_genes = benomyl_sgd.loc[(benomyl_sgd.Phenotype.isin(phenotypes)) & \
							(benomyl_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
ben_genes.to_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

caffeine_sgd = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations.txt", skip_to_line=9, sep="\t").to_pandas()
caffeine_sgd[["Mutant Type", "Mutant Description"]] = caffeine_sgd["Mutant Information"].str.split(":", expand=True)
caf_genes = caffeine_sgd.loc[(caffeine_sgd.Phenotype.isin(phenotypes)) & \
							(caffeine_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
caf_genes.to_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

cuso4_sgd = dt.fread("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations.txt", skip_to_line=9, sep="\t").to_pandas()
cuso4_sgd[["Mutant Type", "Mutant Description"]] = cuso4_sgd["Mutant Information"].str.split(":", expand=True)
cu_genes  = cuso4_sgd.loc[(cuso4_sgd.Phenotype.isin(phenotypes)) & \
							(cuso4_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
cu_genes.to_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

sodmetars_sgd = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations.txt", skiprows=8, sep="\t")
sodmetars_sgd[["Mutant Type", "Mutant Description"]] = sodmetars_sgd["Mutant Information"].str.split(":", expand=True)
sma_genes = sodmetars_sgd.loc[(sodmetars_sgd.Phenotype.isin(phenotypes)) & \
							(sodmetars_sgd["Mutant Type"].isin(mutants)),:].drop_duplicates()
sma_genes.to_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t", index=False)

## Enrichment of literature genes within the feature selection models' features
# Remove literature genes that are not found within our datasets
map_snps = pd.read_csv("~/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
					   sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("~/Shiu_Lab/Project/Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")

ben_genes = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf_genes = dt.fread("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt").to_pandas()
cu_genes = pd.read_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma_genes = pd.read_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")
curated = pd.read_csv("Data/SGD_Experiment_Genes/manually_curated_genes.txt", sep="\t")

ben_snp = ben_genes.loc[ben_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
caf_snp = caf_genes.loc[caf_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
cu_snp = cu_genes.loc[cu_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
sma_snp = sma_genes.loc[sma_genes["Gene Systematic Name"].isin(map_snps.gene),"Gene Systematic Name"].unique()
curated_snp = curated.loc[curated["gene"].isin(map_snps.gene), "gene"].unique()
pd.Series(ben_snp).to_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(caf_snp).to_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(cu_snp).to_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(sma_snp).to_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes_snps.txt", sep="\t", index=False)
pd.Series(curated_snp).to_csv("Data/SGD_Experiment_Genes/manually_curated_genes_snps.txt", sep="\t", index=False)

ben_orf = ben_genes.loc[ben_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
caf_orf = caf_genes.loc[caf_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
cu_orf = cu_genes.loc[cu_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
sma_orf = sma_genes.loc[sma_genes["Gene Systematic Name"].isin(map_orfs.gene),"Gene Systematic Name"].unique()
curated_orf = curated.loc[curated.gene.isin(map_orfs.gene), "gene"].unique()
pd.Series(ben_orf).to_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(caf_orf).to_csv("Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(cu_orf).to_csv("Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(sma_orf).to_csv("Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes_orfs.txt", sep="\t", index=False)
pd.Series(curated_orf).to_csv("Data/SGD_Experiment_Genes/manually_curated_genes_orfs.txt", sep="\t", index=False)

# Add columns to map_snps and map_orfs indicating if the gene is a benchmark gene
map_snps["Benomyl"] = map_snps["gene"].isin(ben_genes["Gene Systematic Name"].values).astype(int)
map_snps["Caffeine"] = map_snps["gene"].isin(caf_genes["Gene Systematic Name"].values).astype(int)
map_snps["CuSO4"] = map_snps["gene"].isin(cu_genes["Gene Systematic Name"].values).astype(int)
map_snps["Sodium_meta-arsenite"] = map_snps["gene"].isin(sma_genes["Gene Systematic Name"].values).astype(int)
map_orfs["Benomyl"] = map_orfs["gene"].isin(ben_genes["Gene Systematic Name"].values).astype(int)
map_orfs["Caffeine"] = map_orfs["gene"].isin(caf_genes["Gene Systematic Name"].values).astype(int)
map_orfs["CuSO4"] = map_orfs["gene"].isin(cu_genes["Gene Systematic Name"].values).astype(int)
map_orfs["Sodium_meta-arsenite"] = map_orfs["gene"].isin(sma_genes["Gene Systematic Name"].values).astype(int)

map_snps.to_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t", index=False)
map_orfs.to_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t", index=False)

# Feature to gene maps with benchmark gene information (from 0_Figure_4_important_genes.py)
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")

ben_snp = map_snps.loc[map_snps.Benomyl == 1, "gene"].unique()
ben_orf = map_orfs.loc[map_orfs.Benomyl == 1, "gene"].unique()
caf_snp = map_snps.loc[map_snps.Caffeine == 1, "gene"].unique()
caf_orf = map_orfs.loc[map_orfs.Caffeine == 1, "gene"].unique()
cu_snp = map_snps.loc[map_snps.CuSO4 == 1, "gene"].unique()
cu_orf = map_orfs.loc[map_orfs.CuSO4 == 1, "gene"].unique()
sma_snp = map_snps.loc[map_snps["Sodium_meta-arsenite"] == 1, "gene"].unique()
sma_orf = map_orfs.loc[map_orfs["Sodium_meta-arsenite"] == 1, "gene"].unique()
curated_snp = pd.read_csv("Data/SGD_Experiment_Genes/manually_curated_genes_snps.txt", sep="\t", header=None)[0].values[1:]
curated_orf = pd.read_csv("Data/SGD_Experiment_Genes/manually_curated_genes_orfs.txt", sep="\t", header=None)[0].values[1:]

# Target environments to test enrichment in
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"] # most predictive envs
lit_gene_envs = ["benomyl", "caffeine", "copper(II) sulfate", "sodium meta-arsenite", "manually curated"]

# Create contingency table
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
def enrichment_direction(k, n, C, G):
    """determine direction of enrichment
    if >= 1: + (overrepresented)
    if < 1: - (underrepresented)
    k: number of genes in top 5% and are known literature genes
    n: total number of genes in top 5%
    C: total number of genes (in top 5% + background) and are known literature genes
    G: total number of genes (in top 5% + background)"""
    return((k/C)/(n/G))


# Check for enrichment at different ranking percentiles
for percentile in [0.99, 0.95, 0.90, 0.85, 0.80, 0.75]:
	# Contingency table header
	enrich_res = [] # collect results from enrichment test
	for data_type in ["snp", "pav", "cnv"]:
		for imp_type in ["imp", "shap"]:
			print(percentile, data_type, imp_type)
			# Read and prepare dataset
			if data_type != "snp":
				df_fs = dt.fread(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}.tsv").to_pandas()
				df_base = dt.fread(f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_{data_type}.tsv").to_pandas()
			else:
				df_fs = dt.fread(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}.tsv").to_pandas()
				df_base = dt.fread(f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_{data_type}.tsv").to_pandas()
			#
			if imp_type == "shap": # take the absolute value of the shap values
				if data_type == "snp":
					df_fs.iloc[:,3:] = df_fs.select_dtypes("float").abs()
					df_base.iloc[:,3:] = df_base.select_dtypes("float").abs()
				else:
					df_fs.iloc[:,2:] = df_fs.select_dtypes("float").abs()
					df_base.iloc[:,2:] = df_base.select_dtypes("float").abs()
			#
			# Count the number of literature genes identified by each model
			env_df_out = []
			for env in target_envs:
				# Add ORF information for those with no gene information
				if data_type != "snp":
					df_fs["gene"] = df_fs.apply(lambda x: x["orf"] if x["gene"] == "" else x["gene"], axis=1)
					df_base["gene"] = df_base.apply(lambda x: x["orf"] if x["gene"] == "" else x["gene"], axis=1)
				
				env_df_fs = df_fs[["gene", env]].dropna() # drop genes/ORFs with no feature importance
				env_df_base = df_base[["gene", env]].dropna()
				
				# Get the feature with the maximum importance value per gene
				env_df_fs = env_df_fs.groupby("gene").max()
				env_df_base = env_df_base.groupby("gene").max()
				if data_type == "snp":
					env_df_fs = env_df_fs.loc[env_df_fs.index != "intergenic",:] # drop intergenic variants
					env_df_base = env_df_base.loc[env_df_base.index != "intergenic",:]
				
				# Add the baseline model genes to the optimized feature selection
				# model genes in order to rank them properly
				env_df = pd.concat([env_df_fs,
					env_df_base.loc[~env_df_base.index.isin(env_df_fs.index),:]], axis=0)
				
				# drop genes with zero feature importance
				env_df = env_df.loc[env_df[env] != 0.0,:]
				
				# Rank the genes by feature importance
				env_df["rank_per"] = env_df.rank(method="average", numeric_only=True, pct=True)
				env_df.sort_values(by="rank_per", ascending=False, inplace=True)
				env_df_out.append(env_df.rename(columns={"rank_per": f"{env}_rank_per"}))
				
			# pd.concat(env_df_out, ignore_index=False, axis=1).to_csv(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}_adjusted_rank_per.tsv", sep="\t") # run lines 119-160 to get these files
				
				env_df_genes = env_df.index.unique()
				#
				# Variable 1: How many genes are ranked in the top #%?
				env_df_top = env_df.loc[env_df.rank_per >= percentile,:] # target gene list
				#
				# Variable 2: How many genes are NOT ranked in the top #%?
				env_df_bg = env_df.loc[env_df.rank_per < percentile,:] # background gene list
				len(env_df) == len(env_df_top) + len(env_df_bg) # sanity check
				print(len(env_df_bg), len(env_df_top), len(env_df))
				#
				# Variable 3: How many genes are known stress-response genes?
				if data_type=="snp":
					num_ben = np.intersect1d(env_df_genes, ben_snp)
					num_caf = np.intersect1d(env_df_genes, caf_snp)
					num_cu = np.intersect1d(env_df_genes, cu_snp)
					num_sma = np.intersect1d(env_df_genes, sma_snp)
					num_man = np.intersect1d(env_df_genes, curated_snp)
				else:
					num_ben = np.intersect1d(env_df_genes, ben_orf)
					num_caf = np.intersect1d(env_df_genes, caf_orf)
					num_cu = np.intersect1d(env_df_genes, cu_orf)
					num_sma = np.intersect1d(env_df_genes, sma_orf)
					num_man = np.intersect1d(env_df_genes, curated_orf)
				#
				# Variable 4: How many genes are NOT known stress-response genes?
				# if data_type=="snp":
				# 	num_not_ben = np.setdiff1d(env_df_genes, ben_snp, assume_unique=True)
				# 	num_not_caf = np.setdiff1d(env_df_genes, caf_snp, assume_unique=True)
				# 	num_not_cu = np.setdiff1d(env_df_genes, cu_snp, assume_unique=True)
				# 	num_not_sma = np.setdiff1d(env_df_genes, sma_snp, assume_unique=True)
				# 	num_not_man = np.setdiff1d(env_df_genes, curated_snp, assume_unique=True)
				# else:
				# 	num_not_ben = np.setdiff1d(env_df_genes, ben_orf, assume_unique=True)
				# 	num_not_caf = np.setdiff1d(env_df_genes, caf_orf, assume_unique=True)
				# 	num_not_cu = np.setdiff1d(env_df_genes, cu_orf, assume_unique=True)
				# 	num_not_sma = np.setdiff1d(env_df_genes, sma_orf, assume_unique=True)
				# 	num_not_man = np.setdiff1d(env_df_genes, curated_orf, assume_unique=True)
				#
				## Create contingency table
				for i,var in enumerate([num_ben, num_caf, num_cu, num_sma, num_man]):
					a = len(env_df_top.loc[env_df_top.index.isin(var),:]) # gene is a known literature gene AND is in top #%
					b = len(env_df_bg.loc[env_df_bg.index.isin(var),:]) # gene is a known literature gene AND is not in top #%
					c = len(env_df_top.loc[~env_df_top.index.isin(var),:]) # gene is not known AND is in top #%
					d = len(env_df_bg.loc[~env_df_bg.index.isin(var),:]) # gene is not known AND is not in top #%
					odds, pval = fisher_exact([[a,b],[c,d]]) # Run fisher's exact test
					#
					# Direction of enrichment
					if enrichment_direction(k=a, n=a+b, C=a+c, G=a+b+c+d) >= 1:
						direction = "+"
					else:
						direction = "-"
					#
					# Save results
					enrich_res.append(["RF_single-env_baseline", data_type, imp_type, env, lit_gene_envs[i], a,\
					b, c, d, odds, pval, np.log2(odds), np.log10(pval), direction])
	#
	enrich_res = pd.DataFrame(enrich_res)
	enrich_res.columns = ["Model Type", "DataType", "ImpType", "Env", "LitGeneList",
					f"Known & above {int(percentile*100)}%",
					f"Not known & above {int(percentile*100)}%",
					f"Known & not above {int(percentile*100)}%",
					f"Not known & not above {int(percentile*100)}%", "Odds Ratio",
					"p-value", "log2(odds ratio)", "log10(p-value)", "direction"]
	enrich_res_sub = enrich_res.loc[enrich_res["p-value"] != 1.0,:]
	q_values = false_discovery_control(enrich_res_sub["p-value"], method="bh") # Adjust the p-values
	enrich_res_sub["q-values"] = q_values
	enrich_res_sub["log10(q-values)"] = enrich_res_sub.apply(lambda x: np.log10(x["q-values"]) if x.direction=="-" else -np.log10(x["q-values"]), axis=1)
	enrich_res_sub.to_csv(f"Scripts/Data_Vis/Section_5/Enrichment_of_literature_genes_in_above_{int(percentile*100)}%_RF.csv", index=False)
	#
	# Visualize the significant enrichment values
	significant = enrich_res_sub.loc[enrich_res_sub["q-values"] <= 0.05, ["DataType", "ImpType", "Env", "LitGeneList", "log2(odds ratio)"]]
	if significant.shape[0] > 0:
		for i_type in significant.ImpType.unique():
			significant_sub = significant.loc[significant.ImpType==i_type,:]
			print(significant)
			print(significant_sub)
			for dat_type in significant_sub.DataType.unique():
				significant_sub2 = significant_sub.loc[significant_sub.DataType==dat_type, \
					["Env", "LitGeneList", "log2(odds ratio)"]].\
					pivot(index="Env", columns="LitGeneList", values="log2(odds ratio)") # pivot wider
				significant_sub2.replace([np.inf, -np.inf], np.nan, inplace=True) # set -Inf to NaN
				plt.figure(figsize=(7, 11))
				sns.heatmap(data=significant_sub2, cmap="RdBu_r", center=0, mask=significant_sub2.isna())
				plt.tight_layout()
				plt.savefig(f"Scripts/Data_Vis/Section_5/Enrichment_of_literature_genes_in_above_{int(percentile*100)}_RF_{dat_type}_{i_type}.pdf")
				plt.close()
			del significant_sub2, significant_sub
