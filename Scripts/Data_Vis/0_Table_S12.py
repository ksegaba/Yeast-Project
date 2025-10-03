#!/usr/bin/env python3

import glob, re, swifter
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

################################################################################
### TABLE S12
################################################################################
# For each data type:
		# total number of unique genes in gene1 and gene2 columns
			# overlap in gene1 and gene2
		# total number of interactions
		# how many literature genes were identified
			# how many interactions were literature genes found to have
			# how many of the known GIs of these literature genes were identified
		# are there "novel" interactions that make sense
	# Comparison of data types:
		# how many interactions overlap between data types

## Map summed SHAP interaction scores to genes

# Feature to gene maps (Mar. 10, 2025: I updated the paths, but I don't know what the whole code is for)
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t") # snp to gene map
map_snps.set_index("snp", inplace=True)
map_snps_dict = map_snps["gene"].to_dict()

map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_16_removed_expanded_benchmark.tsv", sep="\t") # orf to gene map
map_orfs.set_index("orf", inplace=True)
map_orfs_dict = map_orfs["gene"].to_dict()

################################# SHAP Values ##################################
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction"
top_features = {'snp':{}, 'pav':{}, 'cnv':{}} # for scatterplots in next section
kinship = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/kinship.csv", index_col=0)
test = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", sep="\t", header=None)
for env in target_envs:
	print(env)
	snp_shap = dt.fread(glob.glob(f"{d}/SNP/SHAP_values_sorted_average_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap = dt.fread(glob.glob(f"{d}/PAV/SHAP_values_sorted_average_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	cnv_shap = dt.fread(glob.glob(f"{d}/CNV/SHAP_values_sorted_average_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	# map to genes and lit genes
	snp_shap.rename(columns={"C0":"snp", "0":"mean_shap"}, inplace=True)
	pav_shap.rename(columns={"C0":"orf", "0":"mean_shap"}, inplace=True)
	cnv_shap.rename(columns={"C0":"orf", "0":"mean_shap"}, inplace=True)
	snp_shap = snp_shap.loc[snp_shap.mean_shap != 0.0,:] # drop unimportant features
	pav_shap = pav_shap.loc[pav_shap.mean_shap != 0.0,:]
	cnv_shap = cnv_shap.loc[cnv_shap.mean_shap != 0.0,:]
	pav_shap["orf"] = pav_shap.apply(lambda x: re.sub("^X", "", x.orf), axis=1) # fix orf IDs
	pav_shap["orf"] = pav_shap.apply(lambda x: re.sub("\.", "-", x.orf), axis=1)
	cnv_shap["orf"] = cnv_shap.apply(lambda x: re.sub("^X", "", x.orf), axis=1)
	cnv_shap["orf"] = cnv_shap.apply(lambda x: re.sub("\.", "-", x.orf), axis=1)
	snp_shap = snp_shap.set_index("snp").merge(map_snps.set_index("snp"), how="left", left_index=True, right_index=True) # mapping
	pav_shap = pav_shap.set_index("orf").merge(map_orfs.set_index("orf"), how="left", left_index=True, right_index=True)
	cnv_shap = cnv_shap.set_index("orf").merge(map_orfs.set_index("orf"), how="left", left_index=True, right_index=True)
	snp_shap["abs_mean_shap"] = snp_shap.mean_shap.abs() # absolute value of mean shap value
	pav_shap["abs_mean_shap"] = pav_shap.mean_shap.abs()
	cnv_shap["abs_mean_shap"] = cnv_shap.mean_shap.abs()
	snp_shap.sort_values("abs_mean_shap", ascending=False, inplace=True) # sort shap values
	pav_shap.sort_values("abs_mean_shap", ascending=False, inplace=True)
	cnv_shap.sort_values("abs_mean_shap", ascending=False, inplace=True)
	top_features['snp'][env] = snp_shap.index[0]
	top_features['pav'][env] = pav_shap.index[0]
	top_features['cnv'][env] = cnv_shap.index[0]
	# fill in NaNs in gene column with orf name
	pav_shap["gene"] = pav_shap.apply(lambda x: x.name if x.gene is np.nan else x.gene, axis=1)
	cnv_shap["gene"] = cnv_shap.apply(lambda x: x.name if x.gene is np.nan else x.gene, axis=1)
	pav_shap.set_index("gene", inplace=True)
	cnv_shap.set_index("gene", inplace=True)
	snp_shap.set_index("gene", inplace=True)
	# create a heatmap of shap values
	pav_shap.fillna(0, inplace=True)
	cnv_shap.fillna(0, inplace=True)
	# fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6,10)) ## SNP Figure
	# sns.heatmap(snp_shap.iloc[:20,:1], ax=ax[0], square=True, cbar_kws={"orientation":"horizontal"})
	# sns.heatmap(snp_shap.iloc[:20,1:], ax=ax[1], vmin=0, vmax=1, square=True, cbar_kws={"orientation":"horizontal"})
	# plt.tight_layout()
	# plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_snp_{env}_top_20_plus_comb_lit_genes_training.pdf")
	# plt.close()
	# fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6, 10)) ## PAV Figure
	# sns.heatmap(pav_shap.iloc[:20,:1], ax=ax[0], square=True, cbar_kws={"orientation":"horizontal"})
	# sns.heatmap(pav_shap.iloc[:20,1:], ax=ax[1], vmin=0, vmax=1, square=True, cbar_kws={"orientation":"horizontal"})
	# plt.tight_layout()
	# plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_pav_{env}_top_20_plus_comb_lit_genes_training.pdf")
	# plt.close()
	# fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(6, 10)) ## CNV Figure
	# sns.heatmap(cnv_shap.iloc[:20,:1], ax=ax[0], square=True, cbar_kws={"orientation":"horizontal"})
	# sns.heatmap(cnv_shap.iloc[:20,1:], ax=ax[1], vmin=0, vmax=1, square=True, cbar_kws={"orientation":"horizontal"})
	# plt.tight_layout()
	# plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_cnv_{env}_top_20_plus_comb_lit_genes_training.pdf")
	# plt.close()
	# Plot the SHAP value of the top genes vs kinship with W303
	snp_shap_all = pd.read_csv(glob.glob(f"{d}/SNP/SHAP_values_sorted_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
	pav_shap_all = pd.read_csv(glob.glob(f"{d}/PAV/SHAP_values_sorted_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
	cnv_shap_all = pd.read_csv(glob.glob(f"{d}/CNV/SHAP_values_sorted_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
	pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix orf IDs
	pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	cnv_shap_all.columns = cnv_shap_all.apply(lambda x: re.sub("^X", "", x.name), axis=0)
	cnv_shap_all.columns = cnv_shap_all.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	x = kinship.loc[~kinship.index.isin(test[0]), "SACE_GAV"]
	fig, ax = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(8.5,3.5))
	sns.regplot(x=x, y=snp_shap_all.loc[x.index, top_features['snp'][env]], ci=None, ax=ax[0])
	m, b, r, p, sd = stats.linregress(x, snp_shap_all.loc[x.index, top_features['snp'][env]])
	ax[0].annotate(f"Slope: {m:.2f}\nIntercept: {b:.2f}\nPCC: {r:.2f}\nP-value: {p:.2f}\nStd. Err.: {sd:.2f}",
		xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
	ax[0].set_xlabel("Kinship with W303")
	ax[0].set_ylabel(f"SHAP value of top gene {top_features['snp'][env]}")
	#
	sns.regplot(x=x, y=pav_shap_all.loc[x.index, top_features['pav'][env]], ci=None, ax=ax[1])
	m, b, r, p, sd = stats.linregress(x, pav_shap_all.loc[x.index, top_features['pav'][env]])
	ax[1].annotate(f"Slope: {m:.2f}\nIntercept: {b:.2f}\nPCC: {r:.2f}\nP-value: {p:.2f}\nStd. Err.: {sd:.2f}",
		xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
	ax[1].set_xlabel("Kinship with W303")
	ax[1].set_ylabel(f"SHAP value of top gene {top_features['pav'][env]}")
	#
	sns.regplot(x=x, y=cnv_shap_all.loc[x.index, top_features['cnv'][env]], ci=None, ax=ax[2])
	m, b, r, p, sd = stats.linregress(x, cnv_shap_all.loc[x.index, top_features['cnv'][env]])
	ax[2].annotate(f"Slope: {m:.2f}\nIntercept: {b:.2f}\nPCC: {r:.2f}\nP-value: {p:.2f}\nStd. Err.: {sd:.2f}",
		xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
	ax[2].set_xlabel("Kinship with W303")
	ax[2].set_ylabel(f"SHAP value of top gene {top_features['cnv'][env]}")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/{env}_top_gene_shap_vs_w303_kinship.pdf")
	plt.close()

# plot with the top gene for w303
	snp_shap = dt.fread(glob.glob(f"{d}/SNP/SHAP_values_sorted_average_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap = dt.fread(glob.glob(f"{d}/PAV/SHAP_values_sorted_average_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	cnv_shap = dt.fread(glob.glob(f"{d}/CNV/SHAP_values_sorted_average_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()

pd.DataFrame.from_dict(top_features).to_csv("Scripts/Data_Vis/Section_5/shap_value_top_features_per_data_type.csv")

###################### SHAP value trends of select genes #######################
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)
geno_snp = dt.fread("Data/Peter_2018/geno.csv").to_pandas().set_index("ID")
geno_pav = pd.read_csv("Data/Peter_2018/ORFs_pres_abs.csv", index_col=0)
geno_cnv = pd.read_csv("Data/Peter_2018/ORFs_no_NA.csv", index_col=0)
for env in target_envs:
	snp_shap = dt.fread(glob.glob(f"{d}/SNP/SHAP_values_sorted_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap = dt.fread(glob.glob(f"{d}/PAV/SHAP_values_sorted_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	cnv_shap = dt.fread(glob.glob(f"{d}/CNV/SHAP_values_sorted_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], header=True).to_pandas()
	pav_shap.columns = pav_shap.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix orf IDs
	pav_shap.columns = pav_shap.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	cnv_shap.columns = cnv_shap.apply(lambda x: re.sub("^X", "", x.name), axis=0)
	cnv_shap.columns = cnv_shap.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
	# calculate average and get the feature with the max average shap
	snp_top_most_feature = top_features['snp'][env]
	pav_top_most_feature = top_features['pav'][env]
	cnv_top_most_feature = top_features['cnv'][env]
	pav_top_most_feature2 = "X" + re.sub("-", ".", pav_top_most_feature)
	cnv_top_most_feature2 = "X" + re.sub("-", ".", cnv_top_most_feature)
	# plot data type values vs fitness scatterplot colored by shap value
	plt.figure(figsize=(3,3))
	sns.scatterplot(x=geno_snp.loc[snp_shap.ID, snp_top_most_feature].values,\
					y=pheno.loc[snp_shap.ID, env].values,\
					hue=snp_shap.loc[:,snp_top_most_feature].values,\
					palette="RdBu", alpha=0.5)
	plt.title(map_snps.loc[map_snps.snp==snp_top_most_feature,'gene'].values[0])
	plt.legend(frameon=False)
	plt.xlabel("Genotype")
	plt.ylabel("Fitness")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_snp_{env}_top_1_scatterplot.pdf") ## SNP Figure
	plt.close()
	plt.figure(figsize=(3,3))
	sns.scatterplot(x=geno_pav.loc[pav_shap.ID, pav_top_most_feature2].values,\
					y=pheno.loc[pav_shap.ID, env].values,\
					hue=pav_shap.loc[:,pav_top_most_feature].values,\
					palette="RdBu", alpha=0.5)
	try:
		plt.title(map_orfs.loc[map_orfs.orf==pav_top_most_feature,'gene'].values[0])
	except:
		plt.title(pav_top_most_feature)
	plt.legend(frameon=False)
	plt.xlabel("Presence (1) or Absence (0)")
	plt.ylabel("Fitness")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_pav_{env}_top_1_scatterplot.pdf") ## PAV Figure
	plt.close()
	plt.figure(figsize=(3,3))
	sns.scatterplot(x=geno_cnv.loc[cnv_shap.ID, cnv_top_most_feature2].values,\
					y=pheno.loc[cnv_shap.ID, env].values,\
					hue=cnv_shap.loc[:,cnv_top_most_feature].values,\
					palette="RdBu", alpha=0.5)
	try:
		plt.title(map_orfs.loc[map_orfs.orf==cnv_top_most_feature,'gene'].values[0])
	except:
		plt.title(cnv_top_most_feature)
	plt.legend(frameon=False)
	plt.xlabel("Copy number")
	plt.ylabel("Fitness")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_values_sorted_average_cnv_{env}_top_1_scatterplot.pdf") ## CNV Figure
	plt.close()

########################### SHAP Interaction Scores ############################
# Violin plots of interactions scores
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction"
for env in target_envs:	
	print(env)
	try:
		snp = dt.fread(glob.glob(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	except:
		pass # missing benomyl snp summed file
	pav = dt.fread(glob.glob(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	cnv = dt.fread(glob.glob(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	# violin plot of interaction scores
	try:
		sns.violinplot(pd.DataFrame({"snp":snp.Interaction, "pav":pav.Interaction, "cnv":cnv.Interaction}), fill=False)
	except:
		sns.violinplot(pd.DataFrame({"pav":pav.Interaction, "cnv":cnv.Interaction}), fill=False)
	plt.ylim(0,2) # for YPDSODIUMMETAARSENITE
	plt.ylabel("SHAP Interaction Scores")
	plt.tight_layout()
	plt.savefig(f"Scripts/Data_Vis/Section_5/shap_interaction_{env}_violin.pdf") #_ylim.pdf for YPDSODIUMMETAARSENITE
	plt.close()
	try: # for violin plot debugging
		del snp, pav, cnv
	except:
		del pav, cnv

# Count the number of experimentally verified lit genes were identified by SHAP
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt")
map_snps.columns = ["snp", "chr", "pos", "gene"]
map_snps.set_index("snp", inplace=True)
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs.set_index("orf", inplace=True)
for env in target_envs:
	try:
		snp = dt.fread(glob.glob(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	except:
		pass # missing benomyl snp summed file
	pav = dt.fread(glob.glob(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	cnv = dt.fread(glob.glob(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_*_plus_comb_lit_genes_summed.txt")[0]).to_pandas()
	## map feature names to genes
	try:
		snp["Gene1"] = snp.apply(lambda x: map_snps["gene"][x.Feature1], axis=1)
		snp["Gene2"] = snp.apply(lambda x: map_snps["gene"][x.Feature2], axis=1)
	except:
		pass
	# print(snp[["Gene1", "Gene2", "Interaction"]].groupby(["Gene1", "Gene2"]).max().shape)
	# print('snp genes', snp.Gene1.nunique(), snp.Gene2.nunique())
	#
	pav["Feature1"] = pav.apply(lambda x: re.sub("^X", "", x.Feature1), axis=1)
	pav["Feature1"] = pav.apply(lambda x: re.sub("\.", "-", x.Feature1), axis=1)
	pav["Feature2"] = pav.apply(lambda x: re.sub("^X", "", x.Feature2), axis=1)
	pav["Feature2"] = pav.apply(lambda x: re.sub("\.", "-", x.Feature2), axis=1)
	pav["Gene1"] = pav.apply(lambda x: map_orfs["gene"][x.Feature1] \
							if x.Feature1 in map_orfs["gene"].keys() else x.Feature1, axis=1)
	pav["Gene2"] = pav.apply(lambda x: map_orfs["gene"][x.Feature2] \
							if x.Feature2 in map_orfs["gene"].keys() else x.Feature2, axis=1)
	# print(pav[["Gene1", "Gene2", "Interaction"]].groupby(["Gene1", "Gene2"]).max().shape)
	# print('gene1 genes', pd.Series(pav.Gene1.unique()).isin(map_orfs["gene"].values()).sum())
	# print('gene2 genes', pd.Series(pav.Gene2.unique()).isin(map_orfs["gene"].values()).sum())
	# print('gene1 orfs', pav.Gene1.nunique() - pd.Series(pav.Gene1.unique()).isin(map_orfs["gene"].values()).sum())
	# print('gene2 orfs', pav.Gene2.nunique() - pd.Series(pav.Gene2.unique()).isin(map_orfs["gene"].values()).sum())
	cnv["Feature1"] = cnv.apply(lambda x: re.sub("^X", "", x.Feature1), axis=1)
	cnv["Feature1"] = cnv.apply(lambda x: re.sub("\.", "-", x.Feature1), axis=1)
	cnv["Feature2"] = cnv.apply(lambda x: re.sub("^X", "", x.Feature2), axis=1)
	cnv["Feature2"] = cnv.apply(lambda x: re.sub("\.", "-", x.Feature2), axis=1)
	cnv["Gene1"] = cnv.apply(lambda x: map_orfs["gene"][x.Feature1] \
							if x.Feature1 in map_orfs["gene"].keys() else x.Feature1, axis=1)
	cnv["Gene2"] = cnv.apply(lambda x: map_orfs["gene"][x.Feature2] \
							if x.Feature2 in map_orfs["gene"].keys() else x.Feature2, axis=1)
	#
	## number of literature genes identified
	try:
		snp["Gene1_benomyl"] = snp.apply(lambda x: 1 if x.Gene1 in beno["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_benomyl"] = snp.apply(lambda x: 1 if x.Gene2 in beno["Gene Systematic Name"] else 0, axis=1)
		snp["Gene1_caffeine"] = snp.apply(lambda x: 1 if x.Gene1 in caf["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_caffeine"] = snp.apply(lambda x: 1 if x.Gene2 in caf["Gene Systematic Name"] else 0, axis=1)
		snp["Gene1_cuso4"] = snp.apply(lambda x: 1 if x.Gene1 in cu["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_cuso4"] = snp.apply(lambda x: 1 if x.Gene2 in cu["Gene Systematic Name"] else 0, axis=1)
		snp["Gene1_sma"] = snp.apply(lambda x: 1 if x.Gene1 in sma["Gene Systematic Name"] else 0, axis=1) 
		snp["Gene2_sma"] = snp.apply(lambda x: 1 if x.Gene2 in sma["Gene Systematic Name"] else 0, axis=1)
		print('snp gene1 beno', pd.Series(snp.Gene1.unique()).isin(beno["Gene Systematic Name"]).sum(),\
			'snp gene2 beno', pd.Series(snp.Gene2.unique()).isin(beno["Gene Systematic Name"]).sum())
		print('snp gene1 caf', pd.Series(snp.Gene1.unique()).isin(caf["Gene Systematic Name"]).sum(),\
			'snp gene2 caf', pd.Series(snp.Gene2.unique()).isin(caf["Gene Systematic Name"]).sum())
		print('snp gene1 cu', pd.Series(snp.Gene1.unique()).isin(cu["Gene Systematic Name"]).sum(),\
			'snp gene2 cu', pd.Series(snp.Gene2.unique()).isin(cu["Gene Systematic Name"]).sum())
		print('snp gene1 sma', pd.Series(snp.Gene1.unique()).isin(sma["Gene Systematic Name"]).sum(),\
			'snp gene2 sma', pd.Series(snp.Gene2.unique()).isin(sma["Gene Systematic Name"]).sum())
	except:
		pass
	#
	pav["Gene1_benomyl"] = pav.apply(lambda x: 1 if x.Gene1 in beno["Gene Systematic Name"] else 0, axis=1) 
	pav["Gene2_benomyl"] = pav.apply(lambda x: 1 if x.Gene2 in beno["Gene Systematic Name"] else 0, axis=1)
	pav["Gene1_caffeine"] = pav.apply(lambda x: 1 if x.Gene1 in caf["Gene Systematic Name"] else 0, axis=1) 
	pav["Gene2_caffeine"] = pav.apply(lambda x: 1 if x.Gene2 in caf["Gene Systematic Name"] else 0, axis=1)
	pav["Gene1_cuso4"] = pav.apply(lambda x: 1 if x.Gene1 in cu["Gene Systematic Name"] else 0, axis=1) 
	pav["Gene2_cuso4"] = pav.apply(lambda x: 1 if x.Gene2 in cu["Gene Systematic Name"] else 0, axis=1)
	pav["Gene1_sma"] = pav.apply(lambda x: 1 if x.Gene1 in sma["Gene Systematic Name"] else 0, axis=1)
	pav["Gene2_sma"] = pav.apply(lambda x: 1 if x.Gene2 in sma["Gene Systematic Name"] else 0, axis=1)
	print('pav gene1 beno', pd.Series(pav.Gene1.unique()).isin(beno["Gene Systematic Name"]).sum(),\
		'pav gene2 beno', pd.Series(pav.Gene2.unique()).isin(beno["Gene Systematic Name"]).sum())
	print('pav gene1 caf', pd.Series(pav.Gene1.unique()).isin(caf["Gene Systematic Name"]).sum(),\
		'pav gene2 caf', pd.Series(pav.Gene2.unique()).isin(caf["Gene Systematic Name"]).sum())
	print('pav gene1 cu', pd.Series(pav.Gene1.unique()).isin(cu["Gene Systematic Name"]).sum(),\
		'pav gene2 cu', pd.Series(pav.Gene2.unique()).isin(cu["Gene Systematic Name"]).sum())
	print('pav gene1 sma', pd.Series(pav.Gene1.unique()).isin(sma["Gene Systematic Name"]).sum(),\
		'pav gene2 sma', pd.Series(pav.Gene2.unique()).isin(sma["Gene Systematic Name"]).sum())
	#
	cnv["Gene1_benomyl"] = cnv.apply(lambda x: 1 if x.Gene1 in beno["Gene Systematic Name"] else 0, axis=1) 
	cnv["Gene2_benomyl"] = cnv.apply(lambda x: 1 if x.Gene2 in beno["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene1_caffeine"] = cnv.apply(lambda x: 1 if x.Gene1 in caf["Gene Systematic Name"] else 0, axis=1) 
	cnv["Gene2_caffeine"] = cnv.apply(lambda x: 1 if x.Gene2 in caf["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene1_cuso4"] = cnv.apply(lambda x: 1 if x.Gene1 in cu["Gene Systematic Name"] else 0, axis=1) 
	cnv["Gene2_cuso4"] = cnv.apply(lambda x: 1 if x.Gene2 in cu["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene1_sma"] = cnv.apply(lambda x: 1 if x.Gene1 in sma["Gene Systematic Name"] else 0, axis=1)
	cnv["Gene2_sma"] = cnv.apply(lambda x: 1 if x.Gene2 in sma["Gene Systematic Name"] else 0, axis=1)
	#
	try:
		snp.to_csv(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt", index=False, sep="\t")
	except:
		pass
	pav.to_csv(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt", index=False, sep="\t")
	cnv.to_csv(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt", index=False, sep="\t")
	#
	## number of shap interactions for each literature gene
	try:
		snp_beno1 = snp.where(snp.Gene1.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_beno2 = snp.where(snp.Gene2.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_beno1, snp_beno2], ignore_index=False, axis=1)
		snp_caf1 = snp.where(snp.Gene1.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_caf2 = snp.where(snp.Gene2.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_int_count, snp_caf1, snp_caf2], ignore_index=False, axis=1)
		snp_cu1 = snp.where(snp.Gene1.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_cu2 = snp.where(snp.Gene2.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_int_count, snp_cu1, snp_cu2], ignore_index=False, axis=1)
		snp_sma1 = snp.where(snp.Gene1.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
		snp_sma2 = snp.where(snp.Gene2.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
		snp_int_count = pd.concat([snp_int_count, snp_sma1, snp_sma2], ignore_index=False, axis=1)
		snp_int_count.columns = ["NumInt_Gene1_beno", "NumInt_Gene2_beno", \
			"NumInt_Gene1_caffeine", "NumInt_Gene2_caffeine", "NumInt_Gene1_cuso4", \
			"NumInt_Gene2_cuso4", "NumInt_Gene1_sma", "NumInt_Gene2_sma"]
		snp_int_count.fillna(0, inplace=True)
		snp_int_count.index.name = "Gene"
		snp_int_count.to_csv(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes_int_counts.txt", index=False, sep="\t")
	except:
		pass
	pav_beno1 = pav.where(pav.Gene1.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_beno2 = pav.where(pav.Gene2.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_beno1 = cnv.where(cnv.Gene1.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_beno2 = cnv.where(cnv.Gene2.isin(beno["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_beno1, pav_beno2], ignore_index=False, axis=1)
	cnv_int_count = pd.concat([cnv_beno1, cnv_beno2], ignore_index=False, axis=1)
	#
	pav_caf1 = pav.where(pav.Gene1.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_caf2 = pav.where(pav.Gene2.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_caf1 = cnv.where(cnv.Gene1.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_caf2 = cnv.where(cnv.Gene2.isin(caf["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_int_count, pav_caf1, pav_caf2], ignore_index=False, axis=1)
	cnv_int_count = pd.concat([cnv_int_count, cnv_caf1, cnv_caf2], ignore_index=False, axis=1)
	#
	pav_cu1 = pav.where(pav.Gene1.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_cu2 = pav.where(pav.Gene2.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_cu1 = cnv.where(cnv.Gene1.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_cu2 = cnv.where(cnv.Gene2.isin(cu["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_int_count, pav_cu1, pav_cu2], ignore_index=False, axis=1)
	cnv_int_count = pd.concat([cnv_int_count, cnv_cu1, cnv_cu2], ignore_index=False, axis=1)
	#
	pav_sma1 = pav.where(pav.Gene1.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	pav_sma2 = pav.where(pav.Gene2.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	cnv_sma1 = cnv.where(cnv.Gene1.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene1").size()
	cnv_sma2 = cnv.where(cnv.Gene2.isin(sma["Gene Systematic Name"]), axis=0).dropna().groupby("Gene2").size()
	pav_int_count = pd.concat([pav_int_count, pav_sma1, pav_sma2], ignore_index=False, axis=1)
	pav_int_count.columns = snp_int_count.columns
	pav_int_count.fillna(0, inplace=True)
	pav_int_count.index.name = "Gene"
	cnv_int_count = pd.concat([cnv_int_count, cnv_sma1, cnv_sma2], ignore_index=False, axis=1)
	cnv_int_count.columns = snp_int_count.columns
	cnv_int_count.fillna(0, inplace=True)
	cnv_int_count.index.name = "Gene"
	#
	pav_int_count.to_csv(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes_int_counts.txt", index=False, sep="\t")
	cnv_int_count.to_csv(f"{d}/PAV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes_int_counts.txt", index=False, sep="\t")
	#
	try:
		del snp, pav, cnv
	except:
		del pav, cnv

## Compare SHAP-based interactions to known GIs from the literature
# BioGRID database genetic interactions
biogrid = dt.fread("Data/BioGRID/yeast_gi_biogrid.txt").to_pandas()
biogrid = biogrid.iloc[:,[5,6,7,8,11,13,14,17,36]]
biogrid.columns = ["Systematic Name Interactor A", "Systematic Name Interactor B",
				   "Standard Name Interactor A", "Standard name Interactor B",
				   "Evidence", "Author", "PMID", "Throughput", "Organism"]
biogrid = biogrid.loc[biogrid.Organism=="Saccharomyces cerevisiae (S288c)",:]
biogrid = biogrid.loc[biogrid.Evidence.isin(["Synthetic Growth Defect",\
			"Synthetic Lethality", "Synthetic Rescue", "Negative Genetic",\
			"Positive Genetic"]),:] # remove overexpression gene pairs
biogrid_gp = biogrid.apply(lambda x: set([x["Systematic Name Interactor A"], \
			x["Systematic Name Interactor B"]]), axis=1).values # gene pairs
# biogrid_gp = {tuple(value): index for index, value in enumerate(list(biogrid_gp))}
biogrid_gp = {frozenset(set) for set in biogrid_gp}
len(biogrid_gp) # 438546

# Costanzo 2021 benomyl genetic interactions
costanzo = pd.read_excel("Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx",\
	engine="openpyxl", sheet_name="Genome-scale_Benomyl")
# costanzo = dt.fread("Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx/Genome-scale_Benomyl") # xlrd no longer supports xlsx
# costanzo = costanzo.to_pandas()
costanzo = costanzo.iloc[:,0:10]
costanzo = costanzo.loc[(costanzo.mean_condition_epsilon.abs() > 0.12) & \
	(costanzo.condition_p_value < 0.05),:] # Stringent criteria filter
costanzo.insert(0, "array_gene", costanzo.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo.insert(0, "query_gene", costanzo.apply(lambda x: x.query_orf.split("_")[0], axis=1))
costanzo_gp = costanzo.apply(lambda x: set([x.query_gene, x.array_gene]), axis=1).values # gene pairs
costanzo_gp = {frozenset(set) for set in costanzo_gp}
len(costanzo_gp) # 3472
len(biogrid_gp.union(costanzo_gp)) # 440463

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
			   "YPDSODIUMMETAARSENITE"]
d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction"
num_int = {}
for env in target_envs:
	print(env)
	try:
		snp = dt.fread(glob.glob(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt")[0]).to_pandas()
	except:
		pass # missing benomyl snp summed file
	pav = dt.fread(glob.glob(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt")[0]).to_pandas()
	cnv = dt.fread(glob.glob(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes.txt")[0]).to_pandas()#
	#
	# Is there experimental evidence for these SHAP interactions?
	try:
		snp["Known_biogrid"] = snp.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
		snp["Known_costanzo_2021"] = snp.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
		print('snp biogrid', snp["Known_biogrid"].value_counts())
		print('snp costanzo', snp["Known_costanzo_2021"].value_counts())
		num_int[("biogrid", "snp", env)] = snp["Known_biogrid"].value_counts()
		num_int[("costanzo", "snp", env)] = snp["Known_costanzo_2021"].value_counts()
	except:
		pass
	pav["Known_biogrid"] = pav.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
	pav["Known_costanzo_2021"] = pav.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
	print('pav biogrid', pav["Known_biogrid"].value_counts())
	print('pav costanzo', pav["Known_costanzo_2021"].value_counts())
	num_int[("biogrid", "pav", env)] = pav["Known_biogrid"].value_counts()
	num_int[("costanzo", "pav", env)] = pav["Known_costanzo_2021"].value_counts()
	cnv["Known_biogrid"] = cnv.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in biogrid_gp else 0, axis=1)
	cnv["Known_costanzo_2021"] = cnv.swifter.apply(lambda x: 1 if set([x.Gene1, x.Gene2]) in costanzo_gp else 0, axis=1)
	print('cnv biogrid', cnv["Known_biogrid"].value_counts())
	print('cnv costanzo', cnv["Known_costanzo_2021"].value_counts())
	num_int[("biogrid", "cnv", env)] = cnv["Known_biogrid"].value_counts()
	num_int[("costanzo", "cnv", env)] = cnv["Known_costanzo_2021"].value_counts()
	#
	# make sure it's saved as 0s and 1s not a boolean
	try:
		snp[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
			"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
			snp[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
			"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
		snp.to_csv(f"{d}/SNP/shap_interaction_scores_snp_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	except:
		pass
	pav[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
		pav[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
	cnv[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]] = \
		cnv[["Gene1_benomyl", "Gene2_benomyl", "Gene1_caffeine", "Gene2_caffeine",\
		"Gene1_cuso4", "Gene2_cuso4", "Gene1_sma", "Gene2_sma"]].astype(int)
	pav.to_csv(f"{d}/PAV/shap_interaction_scores_pav_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	cnv.to_csv(f"{d}/CNV/shap_interaction_scores_cnv_{env}_top_plus_comb_lit_genes_summed_lit_genes_known.txt", index=False, sep="\t")
	# try:
	# 	del snp, pav, cnv
	# except:
	# 	del pav, cnv

num_int = pd.DataFrame(num_int).T
num_int.fillna(0, inplace=True)

biogrid_sub = num_int.loc[num_int.index.get_level_values(0)=="biogrid",:]
sns.barplot(y = biogrid_sub[1].values,\
	x = biogrid_sub.index.get_level_values(2).values,\
	hue = biogrid_sub.index.get_level_values(1).values,
	palette="viridis")
plt.xticks(rotation=45, size=7)
plt.ylabel("Number of known GIs identified")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/Num_known_GIs_biogrid_identified_barplot.pdf")
plt.close()

costanzo_sub = num_int.loc[num_int.index.get_level_values(0)=="costanzo",:]
sns.barplot(y = costanzo_sub[1].values,\
	x = costanzo_sub.index.get_level_values(2).values,\
	hue = costanzo_sub.index.get_level_values(1).values,
	palette="viridis")
plt.xticks(rotation=45, size=7)
plt.ylabel("Number of known GIs identified")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_5/Num_known_GIs_costanzo_identified_barplot.pdf")
plt.close()