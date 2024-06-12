#!/usr/bin/env python3
"""
Map relevant genes for all environments (especifally Caffeine, Benomyl, and 
CuSO4) to RF importance scores (after FS), ranks, GO terms, pathways, and copy number 
and presence/absence across population.
"""

__author__ = "Kenia Segura Abá"


import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
from glob import glob
from tqdm import tqdm

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

if __name__ == "__main__":
    # Read in gene mapping files
    snp_gene_map = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", header=None) # SNP to gene map
    orf_gene_map = pd.read_csv("Data/Peter_2018/ORFs_gene_map.tsv", sep="\t") # ORF to gene map
    protein_map = pd.read_csv("Data/YeastMine/YeastMine_Protein_07252023.tsv", sep="\t")  # gene to protein map
    orf_presence = pd.read_csv("Data/Peter_2018/percent_iso_ORF_isin.csv") # ORF percent presence across population
    orf_copy = pd.read_csv("Data/Peter_2018/ORFs_no_NA.csv") # ORF copy number across population
    stats_copy_num = orf_copy.describe().T  # summary statistics for copy number
    pwy_gene_map = pd.read_csv("../Co-function/Data/MetaCyc/All-genes-pathways-S288c.txt", sep="\t") # pathway to gene map
    pwy_info = pd.read_csv("../Co-function/Data/MetaCyc/All_pathways_S288c_names.txt", sep="\t") # pathway names

    # Random Forest feature importance files (after feature selection)
    snp_files = [x for x in os.listdir("Scripts/Data_Vis/") if re.search(
        r"[A-Z0-9]+_(exp_rf|rf)_[0-9]+_imp.tsv", str(x))]  # SNP features
    orf_files = [x for x in os.listdir("Scripts/Data_Vis/") if re.search(
        r"[A-Z0-9]+_orf_[0-9]+_imp.tsv", str(x))]  # ORF presence/absence features
    cnv_files = [x for x in os.listdir("Scripts/Data_Vis/") if re.search(
        r"[A-Z0-9]+_cno_[0-9]+_imp.tsv", str(x))]  # ORF copy number features
    snp_files.sort()
    orf_files.sort()
    cnv_files.sort()
    
    # Create excel files
    snp_sheet = pd.ExcelWriter(Path("Scripts/Data_Vis") / "SNP_top_gene_ranks.xlsx")
    snp_sheet2 = pd.ExcelWriter(Path("Scripts/Data_Vis") / "SNP_top10_gene_ranks.xlsx")
    orf_sheet = pd.ExcelWriter(Path("Scripts/Data_Vis") / "ORF_top_gene_ranks.xlsx")
    orf_sheet2 = pd.ExcelWriter(Path("Scripts/Data_Vis") / "ORF_top10_gene_ranks.xlsx")
    cnv_sheet = pd.ExcelWriter(Path("Scripts/Data_Vis") / "CNV_top_gene_ranks.xlsx")
    cnv_sheet2 = pd.ExcelWriter(Path("Scripts/Data_Vis") / "CNV_top10_gene_ranks.xlsx")
    for i in tqdm(range(len(snp_files))):
        # sanity check
        trait = re.sub("_(exp_rf|rf)_[0-9]+_imp.tsv", "", snp_files[i])
        print(trait)
        if (trait==re.sub("_orf_[0-9]+_imp.tsv", "", orf_files[i])==re.sub("_cno_[0-9]+_imp.tsv", "", cnv_files[i])):
            # read feature importance files
            snp_imp = pd.read_csv(Path("Scripts/Data_Vis/") / snp_files[i], sep="\t")
            orf_imp = pd.read_csv(Path("Scripts/Data_Vis/") / orf_files[i], sep="\t")
            cnv_imp = pd.read_csv(Path("Scripts/Data_Vis/") / cnv_files[i], sep="\t")

            # read SHAP value files
            # NACL15M has no SNP SHAP file
            snp_shap = pd.read_csv(glob(
                f"Scripts/Data_Vis/SNP_Figures/SHAP/SHAP_values_sorted_{trait}*")[0], sep="\t")
            orf_shap = pd.read_csv(glob(
                f"Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_presence_absence/SHAP_values_sorted_{trait}*")[0], sep="\t")
            cnv_shap = pd.read_csv(glob(
                f"Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_copy_number/SHAP_values_sorted_{trait}*")[0], sep="\t")
            
            # get SHAP value stats
            snp_shap = snp_shap.describe().T
            orf_shap = orf_shap.describe().T
            cnv_shap = cnv_shap.describe().T
            
            # map SNPs to genes
            snp_shap = snp_gene_map[[0, 3]].merge(snp_shap, left_on=0, right_index=True)
            snp_shap.rename(columns={3: "gene", "mean":"mean_shap", "std":"std_shap"}, inplace=True)
            orf_shap.rename(columns={"mean":"mean_shap", "std":"std_shap"}, inplace=True)
            cnv_shap.rename(columns={"mean":"mean_shap", "std":"std_shap"}, inplace=True)

            # add SHAP values
            snp_imp = snp_imp.merge(snp_shap[["gene", "mean_shap", "std_shap"]], left_on="gene", right_on="gene")
            orf_imp = orf_imp.merge(orf_shap[["mean_shap", "std_shap"]], left_on="X", right_index=True)
            cnv_imp = cnv_imp.merge(cnv_shap[["mean_shap", "std_shap"]], left_on="X", right_index=True)

            # drop SNPs and ORFs with no gene mapping
            snp_imp = snp_imp.loc[snp_imp.gene!="intergenic", :]
            orf_imp = orf_imp.loc[~orf_imp.gene.isnull(), :]
            cnv_imp = cnv_imp.loc[~cnv_imp.gene.isnull(), :]

            # add ranks
            snp_imp["rank_percent_imp"] = snp_imp["max"].rank(method="first", ascending=False, pct=True)
            orf_imp["rank_percent_imp"] = orf_imp["mean_imp"].rank(method="first", ascending=False, pct=True)
            cnv_imp["rank_percent_imp"] = cnv_imp["mean_imp"].rank(method="first", ascending=False, pct=True)

            # rank percent by SHAP values
            snp_imp["rank_percent_shap"] = snp_imp["mean_shap"].rank(method="first", ascending=False, pct=True)
            orf_imp["rank_percent_shap"] = orf_imp["mean_shap"].rank(method="first", ascending=False, pct=True)
            cnv_imp["rank_percent_shap"] = cnv_imp["mean_shap"].rank(method="first", ascending=False, pct=True)
            
            # sort feature importance rank percent
            snp_imp.sort_values("rank_percent_imp", inplace=True)
            orf_imp.sort_values("rank_percent_imp", inplace=True)
            cnv_imp.sort_values("rank_percent_imp", inplace=True)

            # combine
            # combined_outer = snp_imp.merge(orf_imp, how="outer", on="gene")
            # combined_outer = combined_outer = combined_outer.merge(cnv_imp, how="outer", on="gene")
            # combined_outer.columns = ["count_snp", "max_imp_snp", "rank_percent_imp_snp", "X_orf",
            #                           "mean_imp_orf", "rank_percent_imp_orf", "X_cnv", "mean_imp_cnv", "rank_percent_imp_cnv"]
            # combined_outer.reset_index(inplace=True)
            # combined_outer.shape

            # map to protein
            # combined_outer = combined_outer.merge(
            #     protein_map[["Genes Systematic Name", "Protein Standard Name"]],
            #     how="left", left_on="gene", right_on="Genes Systematic Name")
            snp_imp = snp_imp.merge(
                protein_map[["Genes Systematic Name", "Protein Standard Name"]],
                how="left", left_on="gene", right_on="Genes Systematic Name")
            orf_imp = orf_imp.merge(
                protein_map[["Genes Systematic Name", "Protein Standard Name"]],
                how="left", left_on="gene", right_on="Genes Systematic Name")
            cnv_imp = cnv_imp.merge(
                protein_map[["Genes Systematic Name", "Protein Standard Name"]],
                how="left", left_on="gene", right_on="Genes Systematic Name")
            
            # map to pathways
            # combined_outer = combined_outer.merge(
            #     pwy_gene_map[["Accession-1", "Pathways of gene"]], how="left",
            #     left_on="Genes Systematic Name", right_on="Accession-1")
            snp_imp = snp_imp.merge(
                pwy_gene_map[["Accession-1", "Pathways of gene"]], how="left",
                left_on="gene", right_on="Accession-1")
            orf_imp = orf_imp.merge(
                pwy_gene_map[["Accession-1", "Pathways of gene"]], how="left",
                left_on="gene", right_on="Accession-1")
            cnv_imp = cnv_imp.merge(
                pwy_gene_map[["Accession-1", "Pathways of gene"]], how="left",
                left_on="gene", right_on="Accession-1")

            # map genes to ORFs
            # combined_outer = combined_outer.merge(
            #     orf_gene_map[["qacc", "gene"]], how="left", on="gene")
            snp_imp = snp_imp.merge(
                orf_gene_map[["gene", "qacc"]], how="left", on="gene")

            # add percent presence across population for each orf
            # combined_outer = combined_outer.merge(
            #     orf_presence, how="left", left_on="qacc", right_on="orf")
            snp_imp = snp_imp.merge(
                orf_presence, how="left", left_on="qacc", right_on="orf")
            orf_imp = orf_imp.merge(
                orf_presence, how="left", left_on="X", right_on="orf")
            cnv_imp = cnv_imp.merge(
                orf_presence, how="left", left_on="X", right_on="orf")
            
            # add copy number across population for each orf
            stats_copy_num.rename(columns={"mean": "mean_copy_num",
                                    "std": "std_copy_num", "min": "min_copy_num",
                                    "25%": "25_copy_num", "50%": "50_copy_num",
                                    "75%": "75_copy_num", "max": "max_copy_num"},
                                    inplace=True)
            # combined_outer = combined_outer.merge(
            #     stats_copy_num.iloc[:, 1:], how="left", left_on="qacc", right_index=True)
            snp_imp = snp_imp.merge(
                stats_copy_num.iloc[:, 1:], how="left", left_on="qacc", right_index=True)
            orf_imp = orf_imp.merge(
                stats_copy_num.iloc[:, 1:], how="left", left_on="X", right_index=True)
            cnv_imp = cnv_imp.merge(
                stats_copy_num.iloc[:, 1:], how="left", left_on="X", right_index=True)
            
            # drop unecessary columns, reorder, & rename columns
            snp_imp.drop(["count", "Genes Systematic Name", "Accession-1", "qacc"],
                        axis=1, inplace=True)
            snp_imp.rename(columns={"max": "max_imp",
                                    "Protein Standard Name":"protein",
                                    "Pathways of gene":"pathways"}, inplace=True)
            orf_imp.drop(["Genes Systematic Name", "Accession-1", "orf"],
                        axis=1, inplace=True)
            orf_imp.rename(columns={"X": "orf", "Protein Standard Name":"protein",
                                    "Pathways of gene":"pathways"}, inplace=True)
            orf_imp = orf_imp[["orf", "mean_imp", "rank_percent_imp", "mean_shap",
                            "std_shap", "rank_percent_shap", "gene", "protein",
                            "pathways", "percent_presence", "mean_copy_num",
                            "std_copy_num", "min_copy_num", "25_copy_num",
                            "50_copy_num", "75_copy_num", "max_copy_num"]]
            cnv_imp.drop(["Genes Systematic Name", "Accession-1", "orf"],
                        axis=1, inplace=True)
            cnv_imp.rename(columns={"X": "orf", "Protein Standard Name":"protein",
                                    "Pathways of gene":"pathways"}, inplace=True)
            cnv_imp = cnv_imp[["orf", "mean_imp", "rank_percent_imp", "mean_shap",
                            "std_shap", "rank_percent_shap", "gene", "protein",
                                "pathways", "percent_presence", "mean_copy_num",
                                "std_copy_num", "min_copy_num", "25_copy_num",
                                "50_copy_num", "75_copy_num", "max_copy_num"]]

            # round numbers
            # tmp = snp_imp.select_dtypes(include=[np.number])
            # snp_imp.loc[:, tmp.columns] = np.round(tmp, decimals=3)
            # tmp = orf_imp.select_dtypes(include=[np.number])
            # orf_imp.loc[:, tmp.columns] = np.round(tmp, decimals=3)
            # tmp = cnv_imp.select_dtypes(include=[np.number])
            # cnv_imp.loc[:, tmp.columns] = np.round(tmp, decimals=3)

            # combine into one excel file and save
            # combined_outer.drop("qacc", axis=1, inplace=True)
            # combined_outer.to_csv(
            #     Path("Scripts/Data_Vis/") /
            #     re.sub("(exp_rf|rf)_[0-9]+_imp.tsv", "snp_orf_combined_metadata.tsv",
            #     snp_files[i]),sep="\t", index=False)
            
            # save sheet to excel file
            snp_imp.to_excel(snp_sheet, sheet_name=trait, float_format="%.3f", index=False)
            snp_imp.nsmallest(10, "rank_percent_imp").to_excel(snp_sheet2, sheet_name=trait, float_format="%.3f", index=False)
            orf_imp.to_excel(orf_sheet, sheet_name=trait, float_format="%.3f", index=False)
            orf_imp.nsmallest(10, "rank_percent_imp").to_excel(orf_sheet2, sheet_name=trait,float_format="%.3f", index=False)
            cnv_imp.to_excel(cnv_sheet, sheet_name=trait, float_format="%.3f", index=False)
            cnv_imp.nsmallest(10, "rank_percent_imp").to_excel(cnv_sheet2, sheet_name=trait, float_format="%.3f", index=False)
        else:
            print("Error: trait names do not match in all files.")

    snp_sheet.save()
    snp_sheet2.save()
    orf_sheet.save()
    orf_sheet2.save()
    cnv_sheet.save()
    cnv_sheet2.save()