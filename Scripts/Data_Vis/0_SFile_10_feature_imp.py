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

''' What are the 16 ORFs that were removed?
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED.tsv", sep="\t")
map_orfs.loc[map_orfs.orf.duplicated(), "orf"]
                                   	Gene	Benomyl	Caffeine	CuSO4	Sodium_meta-arsenite
758       2391-YCR040W_NumOfGenes_2	YCR040W
                                   	YCL066W
1407      3079-YDR385W_NumOfGenes_2	YOR133W
                                   	YDR385W
1810                   3509-YER179W	YER179W
                                   	SPAR_E02530
2562      4306-YHR055C_NumOfGenes_4	YHR055C >> I'm surprised these aren't in the CuSO4 benchmark genes. They're CUP1-1 and CUP1-2
                                   	YHR053C
2564      4307-YHR056C_NumOfGenes_3	YHR056C
                                   	YHR054C
2745      4500-YIL018W_NumOfGenes_2	YIL018W	no	yes	no	no >> this orf is the only one that mapped to benchmark genes
                                   	YFR031C-A	no	yes	no	no
2888      4654-YIL172C_NumOfGenes_6	YJL221C
                                   	YIL172C
3526      5351-YKR059W_NumOfGenes_2	YKR059W
                                   	YJL138C
3777    5617-YLR154C-H_NumOfGenes_4	YLR159C-A
                                   	YLR157C-C
3779      5619-YLR157C_NumOfGenes_4	YLR160C
                                   	YLR158C
3781      5620-YLR161W_NumOfGenes_5	YLR161W
                                   	YLR159W
4059     5919-YLR467W_NumOfGenes_63	YOR396W
                                   	YLR467W
4868      6792-YNR073C_NumOfGenes_3	YNR073C
                                   	YEL070W
5390      7349-YOR394W_NumOfGenes_7	YPL282C
                                   	YOR394W
5660      7624-YPL281C_NumOfGenes_3	YPL281C
                                   	YOR393W
5741      7703-YPR080W_NumOfGenes_2	YPR080W
                                   	YBR118W
Name: orf, dtype: object
'''

###### What are the ORFs that got added to the corrected ORF to gene map? ######
og_map = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t")
new_map = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED.tsv", sep="\t")

# Orfs only in og_map
set(og_map.orf.tolist()) - set(new_map.orf.tolist())
# {'5859-YLR402W', '3317-YER007C-A', '1696-snap_masked-BIF_1-6636', '6462-YNL069C'} # 4 ORFs

# Orfs only in new_map
set(new_map.orf.tolist()) - set(og_map.orf.tolist())

len(set(new_map.orf.tolist()) - set(og_map.orf.tolist())) # 45 ORFs
# {'4883-YJL189W_NumOfGenes_2', '5349-YKR057W', '4311-YHR060W', '3732-YGL087C', '5537-YLR078C_NumOfGenes_2', '7580-YPL232W', '1897-YBL027W_NumOfGenes_2', '5228-YKL156W_NumOfGenes_2', '1922-YBL050W_NumOfGenes_2', '7086-YOR122C', '3055-YDR363W-A', '6976-YOR010C', '7143-YOR182C', '3166-YDR478W', '3122-YDR432W_NumOfGenes_2', '5746-YLR287C-A', '260-augustus_masked-5430-CDT_2', '2782-YDR092W', '5005-YJR099W', '5073-YKL006W', '4379-YHR123W', '6728-YNR010W', '2715-YDR025W', '7292-YOR332W_NumOfGenes_2', '2944-YDR252W', '6383-YMR311C', '6937-YOL139C_NumOfGenes_2', '5299-YKR006C', '7654-YPR028W', '1912-YBL040C_NumOfGenes_2', '2530-YDL083C', '5251-YKL180W', '2846-YDR156W', '6972-YOR008C', '5552-YLR093C', '3075-YDR381W', '2507-YDL059C', '6520-YNL130C', '6505-YNL113W_NumOfGenes_2', '5057-YJR150C', '2753-YDR064W', '4701-YJL001W', '7588-YPL241C', '2529-YDL082W', '7567-YPL218W'}

# 43 of the 45 ORFs in new_map represent 45 out of 75 genes unique to new_map
# this means, that the other 30 genes may have replaced previous gene annotations in og_map
new_map.loc[(new_map.orf.isin(set(new_map.orf.tolist()) - set(og_map.orf.tolist()))) &\
    (new_map.gene.isin(set(new_map.gene.tolist()) - set(og_map.explode("gene").gene.tolist()))),:].shape # (43, 3) 43 out of 45 ORFs

new_map.loc[(new_map.orf.isin(set(new_map.orf.tolist()) - set(og_map.orf.tolist()))) &\
    (new_map.gene.isin(set(new_map.gene.tolist()) - set(og_map.explode("gene").gene.tolist()))),:].gene.nunique() # 43 out of 75 genes

# Genes only in og_map
og_map.gene = og_map.gene.str.split("//") # separate geneA // geneB ...etc. values
og_map.shape # (5877, 3)
og_map.orf.nunique() # 5877
new_map.orf.nunique() # 5918
og_map.explode("gene").gene.nunique() # 5863
new_map.gene.nunique() # 5904
5918 - 5877 # 41 additional ORFs in new_map
5904 - 5863 # 41 additional genes in new_map

set(og_map.explode("gene").gene.tolist()) - set(new_map.gene.tolist())
# {' YHR054C', 'YJL221C ', 'YNL069C', 'YPL281C ', ' YOR393W', 'YLR159C-A ', 'YPR080W ', 'YER007C-A', 'YLR161W ', ' YOR394W', 'YNR073C ', 'YHR056C ', ' YDR545W', ' YHR053C', ' YDR385W', ' YCL066W', 'YOR133W ', 'YOR396W ', ' YEL070W', ' YIL172C', 'YLR160C ', 'YLR402W', 'YKR059W ', ' YLR467W ', ' YBR118W', ' YJL138C', 'YHR055C ', 'YIL018W ', ' YLR159W', ' YLR158C', ' YLR157C-C', 'YPL282C ', 'YCR040W ', ' YFR031C-A'}

len(set(og_map.explode("gene").gene.tolist()) - set(new_map.gene.tolist())) # 34 Genes

# Genes only in new_map
set(new_map.gene.tolist()) - set(og_map.explode("gene").gene.tolist())
# {'YBL027W', 'YIL172C', 'YDL059C', 'SPAR_G02370', 'YLR159C-A', 'YLR157C-C', 'YBR118W', 'YOR133W', 'YDL082W', 'YOR182C', 'YHR055C', 'YPR028W', 'YKL006W', 'YLR078C', 'YHR054C', 'YNR010W', 'YBL050W', 'YDR156W', 'YNL130C', 'YJR099W', 'YLR467W', 'YHR123W', 'YLR161W', 'YPL218W', 'YDR381W', 'YLR093C', 'YPL282C', 'YNL113W', 'YOR394W', 'YKL156W', 'YDR064W', 'YPL232W', 'YDR432W', 'YOR008C', 'YOR332W', 'YDL083C', 'YJL138C', 'YCR040W', 'YPL241C', 'YHR060W', 'SPAR_E02530', 'YDR363W-A', 'YBL040C', 'YCL066W', 'YIL018W', 'YOL139C', 'YDR252W', 'YEL070W', 'YOR010C', 'YGL087C', 'YFR031C-A', 'YLR159W', 'YOR122C', 'YNR073C', 'YLR287C-A', 'YDR025W', 'YJL189W', 'YJL001W', 'YLR158C', 'YPR080W', 'YMR311C', 'YKR059W', 'YJR150C', 'YKR057W', 'YDR385W', 'YOR393W', 'YHR056C', 'YKL180W', 'YKR006C', 'YLR160C', 'YJL221C', 'YPL281C', 'YDR478W', 'YOR396W', 'YHR053C'}

len(set(new_map.gene.tolist()) - set(og_map.explode("gene").gene.tolist())) # 75 Genes
################################################################################

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
