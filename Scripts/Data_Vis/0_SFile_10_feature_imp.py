#!/usr/bin/env python3

################################################################################
### SUPPLEMENTARY DATA FILE 10
###############################################################################

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

## Combine RF FS & baseline model gini importance values for SNPs, PAVs, and CNVs individually
# read feature to gene map files
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv", sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED.tsv", sep="\t")
map_orfs.drop_duplicates(subset="orf", keep=False, inplace=True) # drop orfs that mapped to multiple genes (16 orfs)
map_orfs.to_csv("~/Shiu_Lab/Project/Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t", index=False)
#old: map_orfs = map_orfs.loc[~map_orfs.gene.str.contains("//"),:] # drop orfs that mapped to multiple genes (16 orfs)
#old: map_orfs.to_csv("~/Shiu_Lab/Project/Data/Peter_2018/final_map_orf_to_gene_16_removed.txt", sep="\t", index=False)
map_snps.merge(map_orfs, how="inner", on="gene").gene.nunique() # 5370 shared genes
map_snps["gene_with_intergenic"] = map_snps.apply(lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1)

# paths to feature importance score files
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/fs"
snp_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t")  # SNP FS results
snp_fs_files = [os.path.join(dir, f"{x}_imp") for x in snp_rf_res['ID']]
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/baseline"
snp_baseline_files = [os.path.join(dir, f"{x}_rf_baseline_imp") for x in mapping.keys()]
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results/fs"
pav_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t")  # ORF pres/abs FS results
pav_fs_files = [os.path.join(dir, f"{x}_imp") for x in pav_rf_res['ID']]
cnv_rf_res = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t")  # CNV FS results
cnv_fs_files = [os.path.join(dir, f"{x}_imp") for x in cnv_rf_res['ID']]
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results/baseline"
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

## Combine RF FS & baseline model SHAP values for SNPs, PAVs, and CNVs individually
# paths to feature SHAP value files
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP"
snp_shap_files = [f"{dir}/SNP/fs/{file}" for file in os.listdir(dir + "/SNP/fs") if file.startswith("SHAP_values_sorted_average_Y")]
snp_shap_baseline_files = [f"{dir}/SNP/baseline/{file}" for file in os.listdir(dir + "/SNP/baseline") if file.startswith("SHAP_values_sorted_average_Y")]
pav_shap_files = [f"{dir}/PAV/fs/{file}" for file in os.listdir(dir + "/PAV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
pav_shap_baseline_files = [f"{dir}/PAV/baseline/{file}" for file in os.listdir(dir + "/PAV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_files = [f"{dir}/CNV/fs/{file}" for file in os.listdir(dir + "/CNV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_baseline_files = [f"{dir}/CNV/baseline/{file}" for file in os.listdir(dir + "/CNV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]

# read feature to gene map files
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv", sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv("~/Shiu_Lab/Project/Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")
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
