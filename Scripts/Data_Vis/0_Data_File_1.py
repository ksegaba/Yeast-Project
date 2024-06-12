#!/usr/bin/env python3

import os
import re
import pandas as pd
import datatable as dt

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

################################################################################ #********* done for SNP, PAV, and CNV baseline and FS
### Data File 1
###############################################################################
## Combine RF FS & baseline model gini importance values for SNPs, PAVs, and CNVs individually
# read feature to gene map files
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)
map_snps.merge(map_orfs, how="inner", on="gene").gene.nunique() # 5356 shared genes
map_snps["gene_with_intergenic"] = map_snps.apply(lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1)

# paths to feature importance score files
dir = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs"
snp_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t")  # SNP FS results
snp_fs_files = [os.path.join(dir, f"{x}_imp") for x in snp_rf_res['ID']]
dir = "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/baseline"
snp_baseline_files = [os.path.join(dir, f"{x}_rf_baseline_imp") for x in mapping.keys()]
dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs"
pav_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t")  # ORF pres/abs FS results
pav_fs_files = [os.path.join(dir, f"{x}_imp") for x in pav_rf_res['ID']]
cnv_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t")  # CNV FS results
cnv_fs_files = [os.path.join(dir, f"{x}_imp") for x in cnv_rf_res['ID']]
dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline"
pav_baseline_files = [os.path.join(dir, f"{x}_pav_baseline_imp") for x in mapping.keys()]
cnv_baseline_files = [os.path.join(dir, f"{x}_cnv_baseline_imp") for x in mapping.keys()]

## combine gini importance for all envs per data type
def combine_imp_indiv(imp_files, map=map_snps, dtype="snp", save="", mapping=mapping):
	for i,env in enumerate(mapping.keys()):
		print(env)
		# Read gini importance file
		file = [f for f in imp_files if env in f]
		print(len(file)) # should be 1
		imp = dt.fread(file[0]).to_pandas()
		imp.set_index(imp.iloc[:,0], inplace=True) # feature names as index
		imp = imp.loc[:,"mean_imp"] # use mean gini importances
		imp.rename(env, inplace=True)
		imp = pd.DataFrame(imp)
		if dtype != "snp":
			imp.index = imp.apply(lambda x: re.sub("^X", "", x.name), axis=1) # rename index
			imp.index = imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
		if i == 0:
			merged = imp.copy(deep=True)
		else:
			merged = pd.concat([merged, imp], axis=1, ignore_index=False) # add to dictionary
		del imp
	print(merged.shape)
	# map to genes
	if dtype == "snp":
		merged = map_snps[["snp","gene","gene_with_intergenic"]].merge(merged, how="right", left_on="snp", right_index=True)
	else:
		merged = map_orfs[["orf","gene"]].merge(merged, how="right", left_on="orf", right_index=True)
	merged.to_csv(save, sep="\t", index=False)
	return merged

combine_imp_indiv(snp_fs_files, save="Scripts/Data_Vis/Section_4/RF_FS_imp_snp.tsv")
combine_imp_indiv(pav_fs_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_FS_imp_pav.tsv")
combine_imp_indiv(cnv_fs_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_FS_imp_cnv.tsv")
combine_imp_indiv(snp_baseline_files, save="Scripts/Data_Vis/Section_4/RF_baseline_imp_snp.tsv")
combine_imp_indiv(pav_baseline_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_baseline_imp_pav.tsv")
combine_imp_indiv(cnv_baseline_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_baseline_imp_cnv.tsv")

################################################################################ #********* done for SNP, PAV, and CNV baseline and FS
## Combine RF FS & baseline model SHAP values for SNPs, PAVs, and CNVs individually
# paths to feature SHAP value files
dir = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP"
snp_shap_files = [f"{dir}/SNP/fs/{file}" for file in os.listdir(dir + "/SNP/fs") if file.startswith("SHAP_values_sorted_average_Y")]
snp_shap_baseline_files = [f"{dir}/SNP/baseline/{file}" for file in os.listdir(dir + "/SNP/baseline") if file.startswith("SHAP_values_sorted_average_Y")]
pav_shap_files = [f"{dir}/PAV/fs/{file}" for file in os.listdir(dir + "/PAV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
pav_shap_baseline_files = [f"{dir}/PAV/baseline/{file}" for file in os.listdir(dir + "/PAV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_files = [f"{dir}/CNV/fs/{file}" for file in os.listdir(dir + "/CNV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_baseline_files = [f"{dir}/CNV/baseline/{file}" for file in os.listdir(dir + "/CNV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]

# read feature to gene map files
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)
map_snps["gene_with_intergenic"] = map_snps.apply(lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1) # label intergenic snps

# combine shap values for all envs per data type
def combine_shap_indiv(shap_files, map=map_snps, merged={}, dtype="snp", save="", mapping=mapping):
	for i,env in enumerate(mapping.keys()):
		print(env)
		# Read SHAP file
		file = [f for f in shap_files if env in f]
		print(len(file)) # should be 1
		shap = dt.fread(file[0]).to_pandas()
		shap.set_index(shap.iloc[:,0], inplace=True)
		shap = shap.iloc[:,1:]
		shap.rename(columns={"C1":env}, inplace=True)
		if dtype != "snp":
			shap.index = shap.apply(lambda x: re.sub("^X", "", x.name), axis=1) # rename index
			shap.index = shap.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
		if i == 0:
			merged = shap.copy(deep=True)
		else:
			merged = pd.concat([merged, shap], axis=1, ignore_index=False) # add to dictionary
		del shap
	print(merged.shape)
	# map to genes
	if dtype == "snp":
		merged = map_snps[["snp","gene","gene_with_intergenic"]].merge(merged, how="right", left_on="snp", right_index=True)
	else:
		merged = map_orfs[["orf","gene"]].merge(merged, how="right", left_on="orf", right_index=True)
	merged.to_csv(save, sep="\t", index=False)
	return merged

combine_shap_indiv(snp_shap_files, save="Scripts/Data_Vis/Section_4/RF_FS_shap_snp.tsv")
combine_shap_indiv(pav_shap_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_FS_shap_pav.tsv")
combine_shap_indiv(cnv_shap_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_FS_shap_cnv.tsv")
combine_shap_indiv(snp_shap_baseline_files, save="Scripts/Data_Vis/Section_4/RF_baseline_shap_snp.tsv")
combine_shap_indiv(pav_shap_baseline_files, map=map_orfs, dtype="pav", save="Scripts/Data_Vis/Section_4/RF_baseline_shap_pav.tsv")
combine_shap_indiv(cnv_shap_baseline_files, map=map_orfs, dtype="cnv", save="Scripts/Data_Vis/Section_4/RF_baseline_shap_cnv.tsv")




## Combine feature importance information from different environments into one table including all data types
# def combine_imp(imp_files, merged={}, merged_50={}, dtype="snp", mapping=mapping):
	# for env in mapping.keys():
	#     if env not in merged.keys():
	#         merged[env]={}
	#     if env not in merged_50.keys():
	#         merged_50[env]={}
	#     file = [f for f in imp_files if env in f] # Read feature importance score files
	#     imp = pd.read_csv(file[0], sep="\t", skiprows=1, names=["X", "mean_imp"])
	#     if dtype == "snp": # Map markers to genes and drop intergenic SNPs
	#         imp = map_snps[["gene", "snp"]].merge(imp, left_on="snp", right_on="X")
	#         imp = imp[imp["gene"] != "intergenic"]
	#     elif dtype != "snp": # Map ORFs to genes and drop those with no mapping
	#         imp.X = imp.X.apply(lambda x: re.sub("^X", "", x))
	#         imp.X = imp.X.apply(lambda x: re.sub("\.", "-", x))
	#         imp = map_orfs.merge(imp, left_on="orf", right_on="X")
	#     imp = imp.groupby("gene")["mean_imp"].max().reset_index()  # Get max importance score per gene
	#     imp.rename(columns={"mean_imp":"max_mean_imp"}, inplace=True)
	#     imp = imp.sort_values(by="max_mean_imp", ascending=False) # Rank genes
	#     imp["rank"] = imp.max_mean_imp.rank(method="average", ascending=False)
	#     if dtype == "snp": # Calculate rank percentiles
	#         imp["rank_per"] = imp["rank"] / num_genes_snps ### note 1/22/2024: I need to change to 1 - n/N because my bins below identify >=.95 as top genes. Right now, I was identifying low rank genes as top
	#     elif dtype != "snp":
	#         imp["rank_per"] = imp["rank"] / num_genes_orfs
	#     imp["rank_per_bin"] = pd.cut(imp["rank_per"],
	#         bins=[0, 0.01, 0.05, 0.1, 1],
	#         labels=["(0.99, 1]", "(0.95, 0.99]", "(0.9, 0.95]", "(0, 0.90]"]) # Bin percentiles into categories
	#     if dtype not in merged[env].keys():
	#         merged[env][dtype] = imp[["gene", "max_mean_imp", "rank_per", "rank_per_bin"]] # Merge feature importance measures
	#     if dtype not in merged_50[env].keys():
	#         imp_50 = imp.iloc[:50,:] # Top 50 genes
	#         new_idx = imp_50.gene.str.split(" // ", expand=True).stack().reset_index(level=1, drop=True).reset_index()
	#         new_idx.columns = ["index", "gene"]
	#         imp_50 = new_idx.merge(imp_50, left_on="gene", right_on="gene")
	#         imp_50.drop(columns=["index", "rank"], inplace=True)
	#         merged_50[env][dtype] = imp_50
#     return merged, merged_50


# def dict2df(nested_dict, idx_col="gene"):
#     """
#     Return a pandas dataframe with a column multi-index and and a row index
#     Arguments:
#             nested_dict: nested dictionary with pandas dataframes in the
#                             innermost values. The keys will become the column
#                             multi-index.
#             idx_col: column to set as row index.
#     """
#     new_dict = {}
#     for outer_key, inner_dict in nested_dict.items():# Iterate through the nested dictionary
#         for inner_key, df in inner_dict.items(): # Append the current dataframe and corresponding index information
#             new_dict[(outer_key, inner_key)] = df.set_index(idx_col)
#     result_df = pd.concat(new_dict, axis=1, ignore_index=False) # Concatenate the dataframes along the new multi-index
#     return result_df

# # combine feature importance data
# snp_baseline_merged, snp_base_merged_50 = combine_imp(snp_baseline_files, merged={}, merged_50={})
# orf_baseline_merged, orf_base_merged_50 = combine_imp(orf_baseline_files, merged=snp_baseline_merged, merged_50=snp_base_merged_50, dtype="orf")
# cnv_baseline_merged, cnv_base_merged_50 = combine_imp(cnv_baseline_files, merged=orf_baseline_merged, merged_50=orf_base_merged_50, dtype="cnv")
# snp_merged, snp_merged_50 = combine_imp(snp_files, merged={}, merged_50={})
# orf_merged, orf_merged_50 = combine_imp(orf_files, merged=snp_merged, merged_50=snp_merged_50, dtype="orf")
# cnv_merged, cnv_merged_50 = combine_imp(cnv_files, merged=orf_merged, merged_50=orf_merged_50, dtype="cnv")
# fs_combined = dict2df(cnv_merged) ; fs_combined_50 = dict2df(cnv_merged_50)
# baseline_combined = dict2df(cnv_baseline_merged) ; base_combined_50 = dict2df(cnv_base_merged_50)
# fs_combined.to_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types.tsv", sep="\t")
# fs_combined_50.to_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types_top50.tsv", sep="\t")
# baseline_combined.to_csv("Scripts/Data_Vis/RF_baseline_imp_all_data_types.tsv", sep="\t")
# base_combined_50.to_csv("Scripts/Data_Vis/RF_baseline_imp_all_data_types_top50.tsv", sep="\t")

# # correlations between baseline feature importances for different datatypes (how to do this for FS models?)
# imp = pd.read_csv("Scripts/Data_Vis/RF_baseline_imp_all_data_types.tsv", sep="\t", index_col=0, header=[0,1,2])
# corr_df = {"env":[], "comparison":[], "gene_num":[], "r":[], "p-val":[]}
# for env in imp.columns.get_level_values(0).unique():
#     snp_vs_pav = pd.concat([imp[env, "snp", "max_mean_imp"], imp[env, "orf", "max_mean_imp"]], axis=1).dropna()
#     snp_vs_cnv = pd.concat([imp[env, "snp", "max_mean_imp"], imp[env, "cnv", "max_mean_imp"]], axis=1).dropna()
#     pav_vs_cnv = pd.concat([imp[env, "orf", "max_mean_imp"], imp[env, "cnv", "max_mean_imp"]], axis=1).dropna()
#     r, p = pearsonr(snp_vs_pav.iloc[:,0], snp_vs_pav.iloc[:,1])
#     corr_df["env"].append(env)
#     corr_df["comparison"].append("snp_vs_pav")
#     corr_df["gene_num"].append(snp_vs_pav.shape[0])
#     corr_df["r"].append(r)
#     corr_df["p-val"].append(p)
#     r, p = pearsonr(snp_vs_cnv.iloc[:,0], snp_vs_cnv.iloc[:,1])
#     corr_df["env"].append(env)
#     corr_df["comparison"].append("snp_vs_cnv")
#     corr_df["gene_num"].append(snp_vs_cnv.shape[0])
#     corr_df["r"].append(r)
#     corr_df["p-val"].append(p)
#     r, p = pearsonr(pav_vs_cnv.iloc[:,0], pav_vs_cnv.iloc[:,1])
#     corr_df["env"].append(env)
#     corr_df["comparison"].append("pav_vs_cnv")
#     corr_df["gene_num"].append(pav_vs_cnv.shape[0])
#     corr_df["r"].append(r)
#     corr_df["p-val"].append(p)
#     del snp_vs_pav, snp_vs_cnv, pav_vs_cnv

# pd.DataFrame.from_dict(corr_df).sort_values('r', ascending=False).\
#     to_csv("Scripts/Data_Vis/RF_baseline_imp_data_type_correlations.tsv", sep="\t", index=False)






### Don't think I need below anymore 12/18/2023; at least I know it needs to be modified
## Combine RF FS and baseline model SHAP values for SNPs, PAVs, and CNVs by gene
# combine SHAP values into one table
# def combine_shap(shap_files, merged={}, merged_50={}, dtype="snp", drop_intergenic:bool=True, baseline:bool=False, mapping=mapping):
#     for env in mapping.keys():
#         print(env)
#         if env not in merged.keys():
#             merged[env]={}
#         if env not in merged_50.keys():
#             merged_50[env]={}
#         file = [f for f in shap_files if env in f] # Read SHAP files
#         if not file:
#             print(f"No file found for {env}")
#             return
#         else:
#             if baseline: # Read in data
#                 shap = dt.fread(file[0]).to_pandas()
#                 shap.set_index(shap.iloc[:,0], inplace=True)
#                 shap = shap.iloc[:,1:]
#                 shap = shap.T
#             else:
#                 shap = pd.read_csv(file[0], sep="\t", header=0, index_col=0).T
#             if dtype == "snp": # Map markers to genes and drop intergenic SNPs
#                 shap = map_snps[["gene", "snp"]].merge(shap, left_on="snp", right_index=True)
#                 if drop_intergenic:
#                     shap = shap[shap["gene"] != "intergenic"]
#                 else: # Rename intergenic snps to intergenic_<snp>
#                     shap["gene"] = shap.apply(lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1)
#                 shap.drop(columns="snp", inplace=True)
#             if dtype != "snp": # Map ORFs to genes and drop those with no mapping
#                 shap.index = shap.apply(lambda x: re.sub("^X", "", x.name), axis=1)
#                 shap.index = shap.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
#                 shap = map_orfs.merge(shap, left_on="orf", right_index=True)
#                 shap.drop(columns=["organism", "orf"], inplace=True)
#             shap.set_index("gene", inplace=True)
#             shap = shap.mean(axis=1)
#             shap = shap.groupby("gene").max().reset_index() # Get max importance score per gene
#             shap.rename(columns={0:"max_mean_shap"}, inplace=True)
#             shap = shap.sort_values(by="max_mean_shap", ascending=False)
#             shap["rank"] = shap.max_mean_shap.rank(method="average", ascending=False) # Rank genes
#             if dtype == "snp": # Calculate rank percentiles
#                 if drop_intergenic:
#                     shap["rank_per"] = shap["rank"] / num_genes_snps ### note 1/22/2024: I need to change to 1 - n/N because my bins below identify >=.95 as top genes. Right now, I was identifying low rank genes as top
#                 else:
#                     shap["rank_per"] = shap["rank"] / num_genes_snps_with_intergenic
#             elif dtype != "snp":
#                 shap["rank_per"] = shap["rank"] / num_genes_orfs
#             shap["rank_per_bin"] = pd.cut(shap["rank_per"],
#                 bins=[0, 0.01, 0.05, 0.1, 1],
#                 labels=["(0.99, 1]", "(0.95, 0.99]", "(0.9, 0.95]", "(0, 0.90]"]) # Bin percentiles into categories
#             if dtype not in merged[env].keys():
#                 merged[env][dtype] = shap.drop(columns="rank")
#             if dtype not in merged_50[env].keys():
#                 shap_50 = shap.iloc[:50,:] # Top 50 genes
#                 new_idx = shap_50.gene.str.split(" // ", expand=True).stack().reset_index(level=1, drop=True).reset_index()
#                 new_idx.columns = ["index", "gene"]
#                 shap_50 = new_idx.merge(shap_50, left_on="gene", right_on="gene")
#                 shap_50.drop(columns=["index", "rank"], inplace=True)
#                 merged_50[env][dtype] = shap_50
#     return merged, merged_50

# # Merge single-env RF FS shap results into a master table
# snp_shap_merged, snp_shap_merged_50 = combine_shap(snp_shap_files, merged={}, merged_50={}, dtype="snp")
# pav_shap_merged, pav_shap_merged_50 = combine_shap(pav_shap_files, merged=snp_shap_merged, merged_50=snp_shap_merged_50, dtype="pav")
# cnv_shap_merged, cnv_shap_merged_50 = combine_shap(cnv_shap_files, merged=pav_shap_merged, merged_50=pav_shap_merged_50, dtype="cnv")
# dict2df(cnv_shap_merged).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types.tsv", sep="\t") # function on line 252
# dict2df(cnv_shap_merged_50).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types_top50.tsv", sep="\t")
# snp_shap_merged_with_intergenic, snp_shap_merged_50_with_intergenic = combine_shap(snp_shap_files, merged={}, merged_50={}, dtype="snp", drop_intergenic=False)
# pav_shap_merged_with_intergenic, pav_shap_merged_50_with_intergenic = combine_shap(pav_shap_files, merged=snp_shap_merged_with_intergenic, merged_50=snp_shap_merged_50_with_intergenic, dtype="pav", drop_intergenic=False)
# cnv_shap_merged_with_intergenic, cnv_shap_merged_50_with_intergenic = combine_shap(cnv_shap_files, merged=pav_shap_merged_with_intergenic, merged_50=pav_shap_merged_50_with_intergenic, dtype="cnv", drop_intergenic=False)
# dict2df(cnv_shap_merged_with_intergenic).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types_with_intergenic.tsv", sep="\t")
# dict2df(cnv_shap_merged_50_with_intergenic).to_csv("Scripts/Data_Vis/RF_FS_shap_all_data_types_with_intergenic_top50.tsv", sep="\t")

# # Merge single-env baseline RF shap results into a master table
# snp_shap_merged_base, snp_shap_merged_base_50 = combine_shap(snp_shap_baseline_files, merged={}, merged_50={}, dtype="snp", baseline=True)
# pav_shap_merged_base, pav_shap_merged_base_50 = combine_shap(pav_shap_baseline_files, merged=snp_shap_merged_base, merged_50=snp_shap_merged_base_50, dtype="pav", baseline=True)
# cnv_shap_merged_base, cnv_shap_merged_base_50 = combine_shap(cnv_shap_baseline_files, merged=pav_shap_merged_base, merged_50=pav_shap_merged_base_50, dtype="cnv", baseline=True)
# dict2df(cnv_shap_merged_base).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types.tsv", sep="\t")
# dict2df(cnv_shap_merged_base_50).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types_top50.tsv", sep="\t")
# snp_shap_merged_base_with_intergenic, snp_shap_merged_base_50_with_intergenic = combine_shap(snp_shap_baseline_files, merged={}, merged_50={}, dtype="snp", drop_intergenic=False, baseline=True)
# pav_shap_merged_base_with_intergenic, pav_shap_merged_base_50_with_intergenic = combine_shap(pav_shap_baseline_files, merged=snp_shap_merged_base_with_intergenic, merged_50=snp_shap_merged_base_50_with_intergenic, dtype="pav", drop_intergenic=False, baseline=True)
# cnv_shap_merged_base_with_intergenic, cnv_shap_merged_base_50_with_intergenic = combine_shap(cnv_shap_baseline_files, merged=pav_shap_merged_base_with_intergenic, merged_50=pav_shap_merged_base_50_with_intergenic, dtype="cnv", drop_intergenic=False, baseline=True)
# dict2df(cnv_shap_merged_base_with_intergenic).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types_with_intergenic.tsv", sep="\t")
# dict2df(cnv_shap_merged_base_50_with_intergenic).to_csv("Scripts/Data_Vis/RF_baseline_shap_all_data_types_with_intergenic_top50.tsv", sep="\t")
