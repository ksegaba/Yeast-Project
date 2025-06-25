#!/usr/bin/env python3

import pickle, re, swifter, os
import pandas as pd
import datatable as dt
import numpy as np
from tqdm import tqdm
from sklearn.metrics import r2_score

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
### TABLE S10
################################################################################
"""Generate the maximum gini importance feature subsets from SNP, PAV, CNV
RF baseline models that will be used to build a new set of RF models, which
only contains one feature per gene, for SHAP interaction score analysis or all
variants per gene (only for benomyl 500 ug/ml)."""

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"

ben_genes = pd.read_csv("Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")

for data_type in ["snp", "pav", "cnv"]:
	# Read and prepare dataset
	df = dt.fread(f"Scripts/Data_Vis/Section_4/RF_baseline_imp_{data_type}.tsv").to_pandas()
	#
	if data_type == "snp":
		df.set_index(["snp", "gene"], inplace=True)
		df = df.iloc[:,1:] # remove "gene_with_intergenic" column; ranking will not include intergenic snps
		df = df.loc[df.index.get_level_values("gene") != "intergenic",:] # drop intergenic snps
	else:
		df.gene = df.apply(lambda x: x.orf if x.gene=="" else x.gene, axis=1) # add orfs with no gene matches to gene column
		df.set_index(["orf", "gene"], inplace=True)
	#
	for env in target_envs:
		if env == "YPDBENOMYL500":
			# For benomyl, we will use all variants per gene
			if data_type == "snp":
				ben_feat = df.loc[df.index.get_level_values("gene").\
						isin(ben_genes['Gene Systematic Name']),:].\
						index.get_level_values("snp")
				ben_feat = pd.Series(ben_feat)
			else:
				ben_feat = df.loc[df.index.get_level_values("gene").\
						isin(ben_genes['Gene Systematic Name']),:].\
						index.get_level_values("orf")
				ben_feat = ben_feat.str.replace("-", ".")
				ben_feat = pd.Series(ben_feat).apply(lambda x: "X" + x if not x.startswith("X") else x)
			
			ben_feat.to_csv(
				f"{d}/Features_all_variants_per_gene_benomyl_500ugml_{data_type}.txt",
				index=False, header=False)
			del ben_feat
		
		# Take the max importance (either gini or shap) per gene (all envs)
		env_imp = df.loc[:,env].dropna().sort_values(ascending=False)
		env_imp = env_imp.loc[env_imp != 0.0,:] # remove features with 0 gini importance
		max_features = pd.concat([env_imp.groupby("gene").idxmax().apply(lambda x: x[0]),
								  env_imp.groupby("gene").max()],
								ignore_index=False, axis=1)
		#
		if data_type=="snp":
			max_features.columns = ["snp", "max_imp"]
		else:
			max_features.columns = ["orf", "max_imp"]
			max_features.orf = max_features.apply(lambda x: "X" + x.orf, axis=1)
			max_features.orf = max_features.apply(lambda x: re.sub("-", ".", x.orf), axis=1)
		#
		max_features.to_csv(
			f"{d}/Feature_map_max_gini_from_RF_baseline_imp_{data_type}_{env}.csv")
		max_features.iloc[:,0].to_csv(
			f"{d}/Features_max_gini_from_RF_baseline_imp_{data_type}_{env}.txt",
			index=False, header=False)

################################################################################
### Feature lists for the "only benchmark genes", "only important non-benchmark
### genes", and "benchmark plus plus important non-benchmark genes" models
### The real Table S11?
################################################################################
# feature to gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t", index_col=0)
map_orfs.index = map_orfs.index.str.replace("-", ".", regex=False) # replace . with - in the index names
map_orfs.index = map_orfs.apply(lambda x: "X" + x.name, axis=1) # remove X from the beginning of the index names
map_orfs.reset_index(inplace=True)
map_orfs.rename(columns={"index": "orf"}, inplace=True) # rename index column to orf

# benchmark genes to features maps
ben_snp = map_snps.loc[map_snps.Benomyl == 1, ["snp", "gene"]]
ben_orf = map_orfs.loc[map_orfs.Benomyl == 1, ["orf", "gene"]]
caf_snp = map_snps.loc[map_snps.Caffeine == 1, ["snp", "gene"]]
caf_orf = map_orfs.loc[map_orfs.Caffeine == 1, ["orf", "gene"]]
cu_snp = map_snps.loc[map_snps.CuSO4 == 1, ["snp", "gene"]]
cu_orf = map_orfs.loc[map_orfs.CuSO4 == 1, ["orf", "gene"]]
sma_snp = map_snps.loc[map_snps["Sodium_meta-arsenite"] == 1, ["snp", "gene"]]
sma_orf = map_orfs.loc[map_orfs["Sodium_meta-arsenite"] == 1, ["orf", "gene"]]

# feature importance scores from RF models built with the complete feature sets
snp_gini = dt.fread("Scripts/Data_Vis/Section_4/RF_baseline_imp_snp.tsv").to_pandas()
pav_gini = dt.fread("Scripts/Data_Vis/Section_4/RF_baseline_imp_pav.tsv").to_pandas()
cnv_gini = dt.fread("Scripts/Data_Vis/Section_4/RF_baseline_imp_cnv.tsv").to_pandas()

snp_gini.set_index("snp", inplace=True)
pav_gini.set_index("orf", inplace=True)
cnv_gini.set_index("orf", inplace=True)

pav_gini.index = pav_gini.index.str.replace("-", ".", regex=False) # replace . with - in the index names
cnv_gini.index = cnv_gini.index.str.replace("-", ".", regex=False)
pav_gini.index = "X" + pav_gini.index.astype(str) # add X to the beginning of the index names
cnv_gini.index = "X" + cnv_gini.index.astype(str) # add X to the beginning of the index names

# feature importance scores from RF models built with the optimized feature sets
snp_gini_fs = dt.fread("Scripts/Data_Vis/Section_4/RF_FS_imp_snp.tsv").to_pandas()
pav_gini_fs = dt.fread("Scripts/Data_Vis/Section_4/RF_FS_imp_pav.tsv").to_pandas()
cnv_gini_fs = dt.fread("Scripts/Data_Vis/Section_4/RF_FS_imp_cnv.tsv").to_pandas()

snp_gini_fs.set_index("snp", inplace=True)
pav_gini_fs.set_index("orf", inplace=True)
cnv_gini_fs.set_index("orf", inplace=True)

pav_gini_fs.index = pav_gini_fs.index.str.replace("-", ".", regex=False) # replace . with - in the index names
cnv_gini_fs.index = cnv_gini_fs.index.str.replace("-", ".", regex=False)
pav_gini_fs.index = "X" + pav_gini_fs.index.astype(str) # add X to the beginning of the index names
cnv_gini_fs.index = "X" + cnv_gini_fs.index.astype(str) # add X to the beginning of the index names


# Create feature tables that only contain benchmark genes (one variant per gene)
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
snp_bench_list = [caf_snp, caf_snp, ben_snp, cu_snp, sma_snp]
orf_bench_list = [caf_orf, caf_orf, ben_orf, cu_orf, sma_orf]
for i,env in enumerate(target_envs):
	for data_type in ["snp", "pav", "cnv"]:
		# subset the benchmark genes from the feature importance tables
		if data_type == "snp":
			gini_env = snp_gini.loc[
				snp_gini.index.isin(snp_bench_list[i].snp), ["gene", env]]
			gini_fs_env = snp_gini_fs.loc[
				~snp_gini_fs.index.isin(snp_bench_list[i].snp), ["gene", env]].dropna() # this still has intergenic snps
			
		elif data_type == "pav":
			gini_env = pav_gini.loc[
				pav_gini.index.isin(orf_bench_list[i].orf), ["gene", env]]
			gini_fs_env = pav_gini_fs.loc[
				~pav_gini_fs.index.isin(orf_bench_list[i].orf), ["gene", env]].dropna()
		
		elif data_type == "cnv":
			gini_env = cnv_gini.loc[cnv_gini.index.isin(orf_bench_list[i].orf), ["gene", env]]
			gini_fs_env = cnv_gini_fs.loc[
				~cnv_gini_fs.index.isin(orf_bench_list[i].orf), ["gene", env]].dropna()
			
		# Are their missing gini importance values?
		assert gini_env[env].isna().sum() == 0, f"Missing {data_type} gini importance values for benchmark genes in {env}!"
		assert gini_fs_env[env].isna().sum() == 0, f"Missing {data_type} gini importance values for important non-benchmark genes in {env}!"
		
		# Keep one feature per gene
		best_feat = gini_env.groupby("gene")[env].idxmax()
		best_fs_feat = gini_fs_env.groupby("gene")[env].idxmax()
		
		# Ensure there are no duplicate features (those that mapped to multiple genes)
		assert best_feat.nunique() == len(best_feat), f"Duplicate {data_type} features found for {env}!"
		assert best_fs_feat.nunique() == len(best_fs_feat), f"Duplicate {data_type} features found in FS for {env}!"
		
		# Write the feature lists to files
		best_feat.to_csv(f"{path}/Features_baseline_only_benchmark_genes_{env}_{data_type}.txt", index=False, header=False)
		best_fs_feat.to_csv(f"{path}/Features_fs_only_important_non_bench_genes_{env}_{data_type}.txt", index=False, header=False)
		
		combined = pd.concat([best_feat, best_fs_feat], axis=0)
		assert combined.nunique() == len(combined), f"Duplicate {data_type} features found in combined for {env}!"
		combined.to_csv(f"{path}/Features_important_non_bench_plus_bench_genes_{env}_{data_type}.txt", index=False, header=False)
		
		del best_feat, best_fs_feat, gini_env, gini_fs_env, combined

# No assertion errors were triggered.

################################################################################
### TABLE S11
################################################################################
"""After running features selection on the RF models used for SHAP interaction
analysis, concatenate the literature genes to the top features and run a new set
of RF models. Also run the models with the label randomized. The purpose of
these is to obtain SHAP gene-gene interactions between top genes and literature
genes and also the feature SHAP values. 
First kind: top FS features plus unimportant benomyl, caffeine, cuso4, or sma genes
Second kind: top FS features plus unimportant non-lit genes
Third kind: only top benchmark genes
Fourth kind: only top non-benchmark genes
"""

# SNP and ORF gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t") # some models may include intergenic snps and snps that mapped to multiple genes in B_nl 
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")

# Read in the RF FS performance results
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs"
res = pd.read_csv(f"{d}/RESULTS_reg.txt", sep="\t") # XGB FS performance results
res = res.loc[res.DateTime.str.contains("03-26"),]

ben_snp = map_snps.loc[map_snps.Benomyl == 1, ["snp", "gene"]]
ben_orf = map_orfs.loc[map_orfs.Benomyl == 1, ["orf", "gene"]]
caf_snp = map_snps.loc[map_snps.Caffeine == 1, ["snp", "gene"]]
caf_orf = map_orfs.loc[map_orfs.Caffeine == 1, ["orf", "gene"]]
cu_snp = map_snps.loc[map_snps.CuSO4 == 1, ["snp", "gene"]]
cu_orf = map_orfs.loc[map_orfs.CuSO4 == 1, ["orf", "gene"]]
sma_snp = map_snps.loc[map_snps["Sodium_meta-arsenite"] == 1, ["snp", "gene"]]
sma_orf = map_orfs.loc[map_orfs["Sodium_meta-arsenite"] == 1, ["orf", "gene"]]

# Now create the feature lists of top features + individual literature gene lists
d2 = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models"
# models_to_run_a = open(f"{d}/fs_plus_lit_genes_models_to_run.txt", "a") # list of models to run
# models_to_run_b = open(f"{d}/fs_plus_unimportant_non_lit_genes_models_to_run.txt", "a")
models_to_run_c = open(f"{d}/only_lit_genes_models_to_run.txt", "a")
models_to_run_d = open(f"{d}/only_top_non_lit_genes_models_to_run.txt", "a")

best_fs_models = open(f"{d}/best_rf_fs_models.txt", "a")
test_set = pd.read_csv('~/Shiu_Lab/Project/Data/Peter_2018/Test.txt', header=None)

for data_type in ["snp", "pav", "cnv"]:
	for env in target_envs:
		if env in ["YPDCAFEIN40", "YPDCAFEIN50"]:
			lit_env = "caffeine"
		elif env == "YPDBENOMYL500":
			lit_env = "benomyl"
		elif env == "YPDCUSO410MM":
			lit_env = "cuso4"
		elif env == "YPDSODIUMMETAARSENITE":
			lit_env = "sodium_meta-arsenite"
		print(data_type, env, lit_env)
		
		# Determine which model has the top validation R2 score
		res_env = res.loc[(res.Tag.str.contains(f"{env}") & \
			res.Tag.str.contains(f"{data_type}")),:]
		top = res_env.loc[res_env["r2_val"] == res_env["r2_val"].max(),"FeatureNum"].values[0]
		
		# Write the best training repetition model to a file (need for shap interaction calculation)
		tag = res_env.loc[res_env['r2_val'] == res_env['r2_val'].max(), 'Tag'].values[0]
		scores = pd.read_csv(f"{d}/{tag}_scores.txt", sep="\t", index_col=0) # predicted trait values for each training repetition
		scores_val = scores.loc[~scores.index.isin(test_set[0]),:].drop(columns=["Mean", "stdev"]) # validation set
		rep = int(scores_val.iloc[:,1:].apply(lambda x: r2_score(scores_val["Y"], x), axis=0).idxmax().split("_")[1]) # best training repetition
		best_fs_models.write(f"{d}/{tag}_models_rep_{rep-1}.pkl\n")
		
		# Read in the top model feature importances
		imp = pd.read_csv(f"{d}/{data_type}_{env}_max_gini_from_RF_baseline_imp_top_{top}_imp", index_col=0, sep="\t")
		imp = imp.loc[imp.mean_imp != 0.0,:] # drop unimportant features
		
		# Get literature genes that map to baseline model
		imp_base = pd.read_csv(f"{d2}/{data_type}_{env}_max_gini_from_RF_baseline_imp_imp", index_col=0, sep="\t")
		
		if data_type != "snp":
			imp.index = imp.apply(lambda x: re.sub("^X", "", x.name), axis=1)
			imp.index = imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
			imp_base.index = imp_base.apply(lambda x: re.sub("^X", "", x.name), axis=1)
			imp_base.index = imp_base.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
		
		if lit_env == "benomyl":
			if data_type=="snp":
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(ben_snp.iloc[:,0]),:] # snp features
				T_l = T_l["mean_imp"].to_frame().merge(ben_snp, how="left", left_index=True, right_on="snp")
				T_nl = imp.loc[~imp.index.isin(ben_snp.iloc[:,0]),:] # snp features  not in benomyl genes
				T_nl = T_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				
				# Identify the unimportant model features found in ben_snp and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(ben_snp.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				# features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp") # fs + unimportant lit genes
				features_c = pd.concat([T_l, B_l]).set_index("snp") # lit genes only
				features_d = T_nl.copy(deep=True).set_index("snp") # important non-lit genes only
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(ben_snp.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				# features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp") # fs + unimportant non-lit genes
					
			else:
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(ben_orf.iloc[:,0]),:] # orf features
				T_l = T_l["mean_imp"].to_frame().merge(
					ben_orf, how="left", left_index=True, right_on="orf")
				T_nl = imp.loc[~imp.index.isin(ben_orf.iloc[:,0]),:] # orf features not in benomyl genes
				T_nl = T_nl["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				
				# Identify the unimportant model features found in ben_orf and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(ben_orf.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_a = pd.concat([T_l, T_nl, B_l])
				features_c = pd.concat([T_l, B_l])
				features_d = T_nl.copy(deep=True)
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(ben_orf.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_b = pd.concat([T_l, T_nl, B_nl])
			
		elif lit_env == "caffeine":
			if data_type=="snp":
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(caf_snp.iloc[:,0]),:] # snp features
				T_l = T_l["mean_imp"].to_frame().merge(
					caf_snp, how="left", left_index=True, right_on="snp")
				T_nl = imp.loc[~imp.index.isin(caf_snp.iloc[:,0]),:] # snp features not in caf genes
				T_nl = T_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				
				# Identify the unimportant model features found in caf_snp and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(caf_snp.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				# features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp")
				features_c = pd.concat([T_l, B_l]).set_index("snp")
				features_d = T_nl.copy(deep=True).set_index("snp")
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(caf_snp.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				# features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp")
					
			else:
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(caf_orf.iloc[:,0]),:] # orf features
				T_l = T_l["mean_imp"].to_frame().merge(caf_orf, how="left", left_index=True, right_on="orf")
				T_nl = imp.loc[~imp.index.isin(caf_orf.iloc[:,0]),:] # orf features  not in caf genes
				T_nl = T_nl["mean_imp"].to_frame().merge(map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				
				# Identify the unimportant model features found in caf_orf and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(caf_orf.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_a = pd.concat([T_l, T_nl, B_l])
				features_c = pd.concat([T_l, B_l])
				features_d = T_nl.copy(deep=True)
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(caf_orf.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_b = pd.concat([T_l, T_nl, B_nl])
		
		elif lit_env == "cuso4":
			if data_type=="snp":
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(cu_snp.iloc[:,0]),:] # snp features
				T_l = T_l["mean_imp"].to_frame().merge(
					cu_snp, how="left", left_index=True, right_on="snp")
				T_nl = imp.loc[~imp.index.isin(cu_snp.iloc[:,0]),:] # snp features not in cuso4 genes
				T_nl = T_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				
				# Identify the unimportant model features found in cu_snp and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(cu_snp.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				# features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp")
				features_c = pd.concat([T_l, B_l]).set_index("snp")
				features_d = T_nl.copy(deep=True).set_index("snp")
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(cu_snp.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				# features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp")
					
			else:
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(cu_orf.iloc[:,0]),:] # orf features
				T_l = T_l["mean_imp"].to_frame().merge(
					cu_orf, how="left", left_index=True, right_on="orf")
				T_nl = imp.loc[~imp.index.isin(cu_orf.iloc[:,0]),:] # orf features not in cuso4 genes
				T_nl = T_nl["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				
				# Identify the unimportant model features found in cu_orf and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(cu_orf.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_a = pd.concat([T_l, T_nl, B_l])
				features_c = pd.concat([T_l, B_l])
				features_d = T_nl.copy(deep=True)
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(cu_orf.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_b = pd.concat([T_l, T_nl, B_nl])
				
		elif lit_env == "sodium_meta-arsenite":
			if data_type=="snp":
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(sma_snp.iloc[:,0]),:] # snp features in sma genes
				T_l = T_l["mean_imp"].to_frame().merge(
					sma_snp, how="left", left_index=True, right_on="snp")
				T_nl = imp.loc[~imp.index.isin(sma_snp.iloc[:,0]),:] # snp features not in sma genes
				T_nl = T_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
				
				# Identify the unimportant model features found in sma_snp and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(sma_snp.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp") # snp features not top and in sma genes
				# features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp")
				features_c = pd.concat([T_l, B_l]).set_index("snp")
				features_d = T_nl.copy(deep=True).set_index("snp")
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(sma_snp.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp") # snp features not top and not in sma genes
				# features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp")
					
			else:
				# Determine literature and non-literature genes in the top features
				T_l = imp.loc[imp.index.isin(sma_orf.iloc[:,0]),:] # orf features
				T_l = T_l["mean_imp"].to_frame().merge(
					sma_orf, how="left", left_index=True, right_on="orf")
				T_nl = imp.loc[~imp.index.isin(sma_orf.iloc[:,0]),:] # orf features not in sma genes
				T_nl = T_nl["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				
				# Identify the unimportant model features found in sma_orf and combine with imp
				B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(imp_base.index.isin(sma_orf.iloc[:,0])),:]
				B_l = B_l["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_a = pd.concat([T_l, T_nl, B_l])
				features_c = pd.concat([T_l, B_l])
				features_d = T_nl.copy(deep=True)
				B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) & \
					(~imp_base.index.isin(sma_orf.iloc[:,0])),:]
				B_nl = B_nl["mean_imp"].to_frame().merge(
					map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
				# features_b = pd.concat([T_l, T_nl, B_nl])
			
		# Ensure that the genes are unique
		if data_type != "snp":
			# Replace missing gene values with the orf ID
			# features_a.gene = features_a.gene.fillna(features_a.orf)
			# features_b.gene = features_b.gene.fillna(features_b.orf)
			features_c.gene = features_c.gene.fillna(features_c.orf)
			features_d.gene = features_d.gene.fillna(features_d.orf)
			# features_a.orf = "X" + features_a.orf.str.replace("-", ".")
			# features_b.orf = "X" + features_b.orf.str.replace("-", ".")
			features_c.orf = "X" + features_c.orf.str.replace("-", ".")
			features_d.orf = "X" + features_d.orf.str.replace("-", ".")
			
			features_c = features_c.set_index("orf")
			features_d = features_d.set_index("orf")
		
		# features_a = features_a.loc[features_a.index.isin(features_a.groupby("gene")["mean_imp"].idxmax()),:] # one feature per gene
		# features_b = features_b.loc[features_b.index.isin(features_b.groupby("gene")["mean_imp"].idxmax()),:]
		features_c = features_c.loc[features_c.index.isin(features_c.groupby("gene")["mean_imp"].idxmax()),:]
		features_d = features_d.loc[features_d.index.isin(features_d.groupby("gene")["mean_imp"].idxmax()),:]
		features_d = features_d.loc[~features_d.gene.str.contains(","),:] # EXCLUDE SNPs THAT MAPPED TO MORE THAN ONE GENE, SINCE ONE COULD BE A BENCHMARK GENE
		# if (features_a.gene.nunique() == len(features_a)) & \
			# (features_b.gene.nunique() == len(features_b)) & \
		if (features_c.gene.nunique() == len(features_c)) & \
			(features_d.gene.nunique() == len(features_d)):
			# np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_{lit_env}_lit_genes_info", features_a, fmt="%s")
			# np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_{lit_env}_lit_genes", features_a.iloc[:,1], fmt="%s")
			# models_to_run_a.write(f"{d}/{data_type}_{env}_top_{top}_plus_{lit_env}_lit_genes\n")
			# np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_{lit_env}_genes_info", features_b, fmt="%s")
			# np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_{lit_env}_genes", features_b.iloc[:,1], fmt="%s")
			# models_to_run_b.write(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_{lit_env}_genes\n")
			
			# Before saving, ensure feature_c and features_d have the same number of features
			np.savetxt(f"{d}/{data_type}_{env}_{lit_env}_lit_genes_info", features_c.reset_index(), fmt="%s")
			np.savetxt(f"{d}/{data_type}_{env}_top_non_{lit_env}_genes_info", features_d.reset_index(), fmt="%s")
			
			if len(features_c) > len(features_d): # ensure the same number of features in both models
				features_c = features_c.sort_values("mean_imp", ascending=False).iloc[:len(features_d),:]
			elif len(features_c) < len(features_d):
				features_d = features_d.sort_values("mean_imp", ascending=False).iloc[:len(features_c),:]
			print(len(features_c), len(features_d))
			
			np.savetxt(f"{d}/{data_type}_{env}_{lit_env}_lit_genes", features_c.index, fmt="%s")
			# models_to_run_c.write(f"{d}/{data_type}_{env}_{lit_env}_lit_genes\n")
			np.savetxt(f"{d}/{data_type}_{env}_top_non_{lit_env}_genes", features_d.index, fmt="%s")
			# models_to_run_d.write(f"{d}/{data_type}_{env}_top_non_{lit_env}_genes\n")
		
		del res_env, top, imp, imp_base, T_l, T_nl, B_l, B_nl, features_c, features_d

# models_to_run_a.close()
# models_to_run_b.close()
models_to_run_c.close()
models_to_run_d.close()
best_fs_models.close()

# Randomize the label and add it to the top features + combined literature gene list datasets
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
pheno_train = pheno.loc[~pheno.index.isin(test[0]),:]
dr = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized"
for i in tqdm(range(100)):
	randomized = pd.DataFrame(index=pheno_train.index)
	for env in target_envs:
		randomized[env] = pheno_train[env].sample(frac=1, random_state=i).reset_index(drop=True).values
	randomized = pd.concat([randomized, pheno.loc[pheno.index.isin(test[0]),target_envs]]) # add the test set back
	randomized.to_csv(f"{dr}/pheno_randomized_{i}.csv")
	del randomized
