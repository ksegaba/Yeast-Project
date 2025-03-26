#!/usr/bin/env python3
'''Figure 6: Benchmark gene model performances
Model set 1: 15 models with the FS features plus respective lit genes
Model set 2: 15 models with the FS features plus unimportant non-lit genes combined
'''

import glob, os, re
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, chisquare, ks_2samp

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

# %%
# Model performance figure
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"
res_fs = pd.read_csv(f"{d}/models_fs/RESULTS_reg.txt", sep="\t") # feature selection on RF models with max_gini_from_RF_baseline_imp features
res_lit = pd.read_csv(f"{d}/only_lit_genes_models/RESULTS_reg.txt", sep="\t") # top lit fs features + unimportant lit features
res_nlit = pd.read_csv(f"{d}/only_top_non_lit_genes_models/RESULTS_reg.txt", sep="\t") # top non-lit fs features

res_fs = res_fs.loc[res_fs.DateTime.str.contains("03-26")].reset_index()
res_lit = res_lit.loc[res_lit.DateTime.str.contains("03-26")].reset_index()
res_nlit = res_nlit.loc[res_nlit.DateTime.str.contains("03-26")].reset_index()

res_fs.insert(2, "Data", res_fs.Tag.str.split("_", expand=True)[0]) # add data column
res_lit.insert(2, "Data", res_lit.Tag.str.split("_", expand=True)[0])
res_nlit.insert(2, "Data", res_nlit.Tag.str.split("_", expand=True)[0])

# res_fs feature selection curves
best_res_fs = []
fig, ax = plt.subplots(nrows=5, ncols=3, sharex=True, sharey=True, figsize=(7.5,11))
for i,data_type in enumerate(["snp", "pav", "cnv"]):
	for j,env in enumerate(["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM", "YPDSODIUMMETAARSENITE"]):
		best_res_fs.append(res_fs.iloc[res_fs.loc[(res_fs.Data==data_type) &\
						  (res_fs.Y==env), "r2_val"].idxmax(),:])
		sns.lineplot(res_fs[(res_fs.Data==data_type) & (res_fs.Y==env)],
					 x="FeatureNum", y="r2_test", color="black", ax=ax[j][i])
		sns.lineplot(res_fs[(res_fs.Data==data_type) & (res_fs.Y==env)],
					 x="FeatureNum", y="r2_val", color="red", ax=ax[j][i])
		ax[j][i].set_title(f"{data_type} {env}")
		ax[j][i].set_xlabel("Number of features")
		ax[j][i].set_ylabel("Performance (R2)")

plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/FS_curves_for_max_gini_from_RF_baseline_imp_models.pdf")
plt.close()

best_res_fs = pd.DataFrame(best_res_fs)

# Plot performances of all models
fig, ax = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(7,8))
for i,data_type in enumerate(["snp", "pav", "cnv"]):
	# Figures for best_res_fs
	sns.barplot(best_res_fs[best_res_fs.Data==data_type], x="Y", y="r2_test", hue="Y", ax=ax[0][i])
	ax[0][i].errorbar(x=best_res_fs.loc[best_res_fs.Data==data_type, "Y"],
					  y=best_res_fs.loc[best_res_fs.Data==data_type, "r2_test"],
					  yerr=best_res_fs.loc[best_res_fs.Data==data_type, "r2_test_sd"],
					  fmt='o', color='black', capsize=5)
	ax[0][i].set_title(f"{data_type} optimized RF models")
	ax[0][i].set_xticklabels(best_res_fs[best_res_fs.Data==data_type].Y, size=7, rotation=50, ha="right")
	
	# Figures for FS features + Lit genes
	sns.barplot(res_lit[res_lit.Data==data_type], x="Y", y="r2_test", hue="Y", ax=ax[1][i])
	ax[1][i].errorbar(x=res_lit.loc[res_lit.Data==data_type, "Y"], 
					  y=res_lit.loc[res_lit.Data==data_type, "r2_test"],
					  yerr=res_lit.loc[res_lit.Data==data_type, "r2_test_sd"],
					  fmt='o', color='black', capsize=5)
	ax[1][i].set_title(f"{data_type}: top lit genes + unimportant lit genes")
	ax[1][i].set_xticklabels(res_lit[res_lit.Data==data_type].Y, size=7, rotation=50, ha="right")
	
	# Figures for FS features + Unimportant genes
	sns.barplot(res_nlit[res_nlit.Data==data_type], x="Y", y="r2_test", hue="Y", ax=ax[2][i])
	ax[2][i].errorbar(x=res_nlit.loc[res_nlit.Data==data_type, "Y"],
					  y=res_nlit.loc[res_nlit.Data==data_type, "r2_test"],
					  yerr=res_nlit.loc[res_nlit.Data==data_type, "r2_test_sd"],
					  fmt='o', color='black', capsize=5)
	ax[2][i].set_title(f"{data_type} RF FS non-lit genes only")
	ax[2][i].set_xticklabels(res_nlit[res_nlit.Data==data_type].Y, size=7, rotation=50, ha="right")

plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/Figure_RF_benchmark_gene_model_performances.pdf")
plt.close()

'''Compare the PAV and CNV distributions of benchmark genes to the important
non-benchmark genes-- to explain difference in model performances'''
# Read in the benchmark gene and gene map data
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv",
                       sep="\t", index_col=0)
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv",
                       sep="\t",index_col=0)

# Read in the feature values and the testing set isolates
snp = dt.fread("Data/Peter_2018/geno.csv").to_pandas()
pav = dt.fread("Data/Peter_2018/ORFs_pres_abs.csv").to_pandas()
cnv = dt.fread("Data/Peter_2018/ORFs_no_NA.csv").to_pandas()
snp.set_index("ID", inplace=True)
pav.set_index("ID", inplace=True)
cnv.set_index("ID", inplace=True)
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)

snp_train = snp.loc[~snp.index.isin(test[0]),:] # keep only the training set
pav_train = pav.loc[~pav.index.isin(test[0]),:]
cnv_train = cnv.loc[~cnv.index.isin(test[0]),:]
del snp, pav, cnv

pav_train.columns = pav_train.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix ORF IDs to match map_orfs
pav_train.columns = pav_train.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
cnv_train.columns = cnv_train.apply(lambda x: re.sub("^X", "", x.name), axis=0)
cnv_train.columns = cnv_train.apply(lambda x: re.sub("\.", "-", x.name), axis=0)

# Plot the distributions of feature importance to feature values
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM", "YPDSODIUMMETAARSENITE"]
results = {"Benchmark": {}, "Non-benchmark": {}, "Comparison": {}}
for data_type in ["snp", "pav", "cnv"]:
	for env in target_envs:
		print(data_type, env)
		# Load the feature importances
		only_lit_imp = pd.read_csv(
			glob.glob(f"{d}/only_lit_genes_models/{data_type}_{env}_*_imp")[0],
			sep="\t", index_col=0)
		only_nlit_imp = pd.read_csv(glob.glob(
			f"{d}/only_top_non_lit_genes_models/{data_type}_{env}_*_imp")[0],
			sep="\t", index_col=0)
		
		# separate the benchmark genes from the important non-benchmark genes
		if data_type == "snp":
			x_bench = pd.Series(snp_train.loc[:,only_lit_imp.index].to_numpy().flatten()).value_counts() # feature values
			x_nonbench = pd.Series(snp_train.loc[:,only_nlit_imp.index].to_numpy().flatten()).value_counts()
			
			# What is the relationship between feature values and feature importance?
			rb, pb = pearsonr(snp_train.loc[:,only_lit_imp.index].to_numpy().flatten(),
							  np.tile(only_lit_imp.mean_imp, 625))
			rnb, pnb = pearsonr(snp_train.loc[:,only_nlit_imp.index].to_numpy().flatten(),
								np.tile(only_nlit_imp.mean_imp, 625))
		
		else:
			only_lit_imp.index = only_lit_imp.apply(lambda x: re.sub("^X", "", x.name), axis=1) # fix ORF IDs to match map_orfs
			only_lit_imp.index = only_lit_imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
			only_nlit_imp.index = only_nlit_imp.apply(lambda x: re.sub("^X", "", x.name), axis=1)
			only_nlit_imp.index = only_nlit_imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
			
			if data_type == "pav":
				x_bench = pd.Series(pav_train.loc[:,only_lit_imp.index].to_numpy().flatten()).astype(int).value_counts() # feature values
				x_nonbench = pd.Series(pav_train.loc[:,only_nlit_imp.index].to_numpy().flatten()).astype(int).value_counts()
				
				rb, pb = pearsonr(pav_train.loc[:,only_lit_imp.index].to_numpy().flatten(),
							  np.tile(only_lit_imp.mean_imp, 625))
				rnb, pnb = pearsonr(pav_train.loc[:,only_nlit_imp.index].to_numpy().flatten(),
								np.tile(only_nlit_imp.mean_imp, 625))
			if data_type == "cnv":
				x_bench = pd.Series(cnv_train.loc[:,only_lit_imp.index].to_numpy().flatten()).value_counts() # feature values
				x_nonbench = pd.Series(cnv_train.loc[:,only_nlit_imp.index].to_numpy().flatten()).value_counts()
				
				rb, pb = pearsonr(cnv_train.loc[:,only_lit_imp.index].to_numpy().flatten(),
							  np.tile(only_lit_imp.mean_imp, 625))
				rnb, pnb = pearsonr(cnv_train.loc[:,only_nlit_imp.index].to_numpy().flatten(),
								np.tile(only_nlit_imp.mean_imp, 625))
		
		results["Benchmark"][f"{data_type}_{env}"] = {"r": rb, "p-value":pb}
		results["Non-benchmark"][f"{data_type}_{env}"] = {"r": rnb, "p-value":pnb}
		
		# Are the feature values significantly different for benchmark genes vs non-benchmark genes?
		x_bench.name = "Benchmark" ; x_nonbench.name = "Non-Benchmark"
		contingency_feat_vals = pd.concat([x_bench, x_nonbench], axis=1).fillna(0).astype(int)
		Chix, px = chisquare(f_obs=contingency_feat_vals.Benchmark, f_exp=contingency_feat_vals["Non-Benchmark"]) # compare feature values
		
		KSy, py = ks_2samp(only_lit_imp.mean_imp, only_nlit_imp.mean_imp, alternative="two-sided", method="auto") # compare feature importances
		results["Comparison"][f"{data_type}_{env}"] = {"Chi2_feat_val": Chix,
			"p-value_feat_val": px, "KS_feat_imp": KSy, "p-value_feat_imp": py}
		
		# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4)) # this figure takes too long.
		# sns.kdeplot(x=x_bench, y=y_bench, ax=ax[0], cmap="Blues", fill=True,
		# 	thresh=False, cbar=True, cbar_kws={"label": "Benchmark gene density"})
		# ax[0].set_title(f"{data_type} {env}: Benchmark genes")
		# ax[0].set_xlabel("Feature value")
		# ax[0].set_ylabel("Feature importance")
		# sns.kdeplot(x=x_nonbench, y=y_nonbench, ax=ax[1], cmap="Reds", fill=True,
		# 	thresh=False, cbar=True, cbar_kws={"label": "Non-benchmark gene density"})
		# ax[1].set_title(f"{data_type} {env}: Important non-benchmark genes")
		# ax[1].set_xlabel("Feature value")
		# ax[1].set_ylabel("Feature importance")
		# plt.suptitle(f"r={rb:.2f}, p={pb:.2e} (Benchmark genes) | r={rnb:.2f}, p={pnb:.2e} (Non-benchmark genes)")
		# plt.tight_layout()
		# plt.savefig(
		# 	f"Scripts/Data_Vis/Section_6/Figure_{data_type}_{env}_imp_vs_feature_value.pdf")
		# plt.close()
		
		del only_lit_imp, only_nlit_imp, x_bench, x_nonbench

pd.DataFrame.from_dict({(i, j): results[i][j] for i in results.keys() for j in results[i].keys()},
	orient="index").to_csv(
	f"Scripts/Data_Vis/Section_6/Figure_imp_vs_feature_value_correlations.csv")
