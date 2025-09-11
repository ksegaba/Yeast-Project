#!/usr/bin/env python3

################################################################################
# TABLE S5: SPEARMAN's RHO FEATURE RANK PERCENTILE CORRELATIONS BETWEEN VARIANTS
################################################################################

import os
import pandas as pd
import datatable as dt
from scipy.stats import spearmanr

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Isolate growth condition labels; will be used throughout the script
mapping = {"YPACETATE": "YP Acetate 2%", "YPD14": "YPD 14ºC", "YPD40": "YPD 40ºC",
           "YPD42": "YPD 42ºC", "YPD6AU": "YPD 6-Azauracile 600 µg/ml",
           "YPDANISO10": "YPD Anisomycin 10 µg/ml", "YPDANISO20": "YPD Anisomycin 20 µg/ml",
           "YPDANISO50": "YPD Anisomycin 50 µg/ml", "YPDBENOMYL200": "YPD Benomyl 200 µg/ml",
           "YPDBENOMYL500": "YPD Benomyl 500 µg/ml", "YPDCAFEIN40": "YPD Caffeine 40 mM",
           "YPDCAFEIN50": "YPD Caffeine 50 mM", "YPDCHX05": "YPD Cycloheximide 0.5 µg/ml",
           "YPDCHX1": "YPD Cycloheximide 1 µg/ml", "YPDCUSO410MM": "YPD CuSO4 10 mM",
           "YPDDMSO": "YPD DMSO 6%", "YPDETOH": "YPD Ethanol 15%",
           "YPDFLUCONAZOLE": "YPD Fluconazole 20 µg/ml", "YPDFORMAMIDE4": "YPD Formamide 4%",
           "YPDFORMAMIDE5": "YPD Formamide 5%", "YPDHU": "YPD Hydroxyurea 30 mg/ml",
           "YPDKCL2M": "YPD KCL 2 M", "YPDLICL250MM": "YPD LiCl 250 mM",
           "YPDMV": "YPD Methylviologen 20 mM", "YPDNACL15M": "YPD NaCl 1.5 M",
           "YPDNACL1M": "YPD NaCl 1 M", "YPDNYSTATIN": "YPD Nystatin 10 µg/ml",
           "YPDSDS": "YPD SDS 0.2%", "YPDSODIUMMETAARSENITE": "YPD Sodium metaarsenite 2.5 mM",
           "YPETHANOL": "YP Ethanol 2%", "YPGALACTOSE": "YP Galactose 2%",
           "YPRIBOSE": "YP Ribose 2%", "YPGLYCEROL": "YP Glycerol 2%",
           "YPXYLOSE": "YP Xylose 2%", "YPSORBITOL": "YP Sorbitol 2%"}

# SPEARMAN'S RHO CORRELATIONS BEWTEEN DATA TYPES
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                       sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")

res = [["Model Type", "Importance Type", "Env",
        "Comparison", "rho", "pval", "Num Shared Genes"]]

for mod_type in ["baseline", "FS"]:
    for imp_type in ["imp", "shap"]:
        # _rank_per.tsv").to_pandas()
        df_snp = dt.fread(
            f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_snp.tsv").to_pandas()
        # _rank_per.tsv").to_pandas()
        df_pav = dt.fread(
            f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_pav.tsv").to_pandas()
        # _rank_per.tsv").to_pandas()
        df_cnv = dt.fread(
            f"Scripts/Data_Vis/Section_4/RF_{mod_type}_{imp_type}_cnv.tsv").to_pandas()
        #
        # Calculate spearman's rho
        if imp_type == "imp":
            imp_score_type = "mean gini importance"
        else:
            imp_score_type = "mean absolute SHAP value"
        #
        # Prep data for SNP to PAV/CNV comparisons
        tmp_snp = df_snp.loc[df_snp.gene != "intergenic", :]
        # drop snps that mapped to multiple genes
        tmp_snp = df_snp.loc[~df_snp.gene.str.contains(","), :]
        # drop ORFs with no gene mapping
        tmp_pav = df_pav.loc[df_pav.gene != "", :]
        tmp_cnv = df_cnv.loc[df_cnv.gene != "", :]
        #
        for env in mapping.keys():
            # For sanity check, the SHAP values are average absolute values.
            if tmp_snp[env].any() < 0:
                print("Negative values found in SNP data!")
            if tmp_pav[env].any() < 0:
                print("Negative values found in PAV data!")
            if tmp_cnv[env].any() < 0:
                print("Negative values found in CNV data!")
            #
            # SNP vs PAV
            try:
                df = pd.concat([tmp_snp.loc[:, ["gene", env]].groupby("gene").max(),
                                tmp_pav.loc[:, ["gene", env]].groupby("gene").max()],
                               ignore_index=False, axis=1).dropna()
                df = df.loc[~df.eq(0).any(axis=1)]
                df = df.rank(pct=True)
                if len(df) > 2:
                    rho = df.corr(method=lambda x,
                                  y: spearmanr(x, y).statistic)
                    pval = df.corr(method=lambda x, y: spearmanr(x, y).pvalue)
                    res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs",
                                rho.iloc[0, 1], pval.iloc[0, 1], len(df)])
                else:
                    res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs",
                                None, None, len(df)])
                del df
            except:
                res.append([mod_type, "max " + imp_score_type, env, "SNPs vs PAVs",
                            None, None, None])
            # SNP vs CNV
            try:
                df = pd.concat([tmp_snp.loc[:, ["gene", env]].groupby("gene").max(),
                                tmp_cnv.loc[:, ["gene", env]].groupby("gene").max()],
                               ignore_index=False, axis=1).dropna()
                df = df.loc[~df.eq(0).any(axis=1)]
                df = df.rank(pct=True)
                if len(df) > 2:
                    rho = df.corr(method=lambda x,
                                  y: spearmanr(x, y).statistic)
                    pval = df.corr(method=lambda x, y: spearmanr(x, y).pvalue)
                    res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs",
                                rho.iloc[0, 1], pval.iloc[0, 1], len(df)])
                else:
                    res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs",
                                None, None, len(df)])
                del df
            except:
                res.append([mod_type, "max " + imp_score_type, env, "SNPs vs CNVs",
                            None, None, None])
            # PAV vs CNV (keep all ORFs regardless of gene mapping)
            try:
                pav = df_pav.loc[df_pav.orf != "", ["orf", env]
                                 ].set_index("orf")
                cnv = df_cnv.loc[df_cnv.orf != "", ["orf", env]
                                 ].set_index("orf")
                df = pd.concat([pav, cnv], ignore_index=False, axis=1).dropna()
                df = df.loc[~df.eq(0).any(axis=1)]
                df = df.rank(pct=True)
                if len(df) > 2:
                    rho = df.corr(method=lambda x,
                                  y: spearmanr(x, y).statistic)
                    pval = df.corr(method=lambda x, y: spearmanr(x, y).pvalue)
                    res.append([mod_type, imp_score_type, env, "PAVs vs CNVs",
                                rho.iloc[0, 1], pval.iloc[0, 1], len(df)])
                else:
                    res.append([mod_type, imp_score_type, env, "PAVs vs CNVs",
                                None, None, len(df)])
                del df
            except:
                res.append([mod_type, imp_score_type, env, "PAVs vs CNVs",
                            None, None, None])
        del df_snp, df_pav, df_cnv, tmp_snp, tmp_pav, tmp_cnv

res = pd.DataFrame(res)
res.columns = res.iloc[0, :]
res = res.iloc[1:, :]
res.sort_values(by='rho', ascending=False, inplace=True)
# res.to_csv("Scripts/Data_Vis/Section_4/Table_S5_rank_per_corr_btwn_data_types.tsv", sep="\t", index=False)
res.to_csv("Scripts/Data_Vis/Section_4/Table_S5_rank_per_corr_btwn_data_types_with_n_shared_genes_CORRECTED.tsv",
           sep="\t", index=False)


# Get the average overlap of shared genes between variant types across all 35 environments
res = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Table_S5_rank_per_corr_btwn_data_types_with_n_shared_genes_CORRECTED.tsv", sep="\t")
res.insert(1, "imp_type", res.apply(
    lambda x: "gini" if "gini" in x["Importance Type"] else "SHAP", axis=1))
res.groupby(['Model Type', 'imp_type', 'Comparison'])[
    'Num Shared Genes'].describe()  # all envs
#                                   count         mean          std     min     25%     50%     75%     max
# Model Type imp_type Comparison
# FS         SHAP     PAVs vs CNVs   35.0    33.485714    31.307180     0.0    12.0    28.0    40.5   121.0
#                     SNPs vs CNVs   35.0     2.971429     4.859713     0.0     0.0     1.0     3.0    21.0
#                     SNPs vs PAVs   35.0     2.057143     3.271471     0.0     0.0     1.0     2.0    13.0
#            gini     PAVs vs CNVs   35.0    33.485714    31.307180     0.0    12.0    28.0    40.5   121.0
#                     SNPs vs CNVs   35.0     3.000000     4.970501     0.0     0.0     1.0     3.0    22.0
#                     SNPs vs PAVs   35.0     2.057143     3.271471     0.0     0.0     1.0     2.0    13.0
# baseline   SHAP     PAVs vs CNVs   35.0   906.400000   394.296186    85.0   699.5   922.0  1097.0  1755.0
#                     SNPs vs CNVs   35.0  1515.885714   939.633463    28.0   932.5  1715.0  2086.0  3734.0
#                     SNPs vs PAVs   35.0   296.342857   165.704014     8.0   188.0   313.0   396.0   625.0
#            gini     PAVs vs CNVs   35.0  1990.571429   470.569573   575.0  1889.0  2145.0  2277.0  2473.0
#                     SNPs vs CNVs   35.0  4550.142857  1210.454731  1009.0  4789.0  5091.0  5242.5  5359.0
#                     SNPs vs PAVs   35.0   655.857143   157.088050   158.0   682.5   728.0   738.5   747.0

# Get the average overlap for target environments only
res_sub = res[res.Env.isin(['YPDCAFEIN40', 'YPDCAFEIN50', 'YPDBENOMYL500',
                            'YPDCUSO410MM', 'YPDSODIUMMETAARSENITE'])]  # target envs only
res_sub.groupby(['Model Type', 'imp_type', 'Comparison'])[
    'Num Shared Genes'].describe()
#                                   count    mean          std     min     25%     50%     75%     max
# Model Type imp_type Comparison
# FS         SHAP     PAVs vs CNVs    5.0    38.4    22.941229     8.0    28.0    39.0    47.0    70.0
#                     SNPs vs CNVs    5.0     5.4     4.037326     1.0     3.0     4.0     8.0    11.0
#                     SNPs vs PAVs    5.0     7.2     5.761944     1.0     2.0     7.0    13.0    13.0
#            gini     PAVs vs CNVs    5.0    38.4    22.941229     8.0    28.0    39.0    47.0    70.0
#                     SNPs vs CNVs    5.0     5.4     4.037326     1.0     3.0     4.0     8.0    11.0
#                     SNPs vs PAVs    5.0     7.2     5.761944     1.0     2.0     7.0    13.0    13.0
# baseline   SHAP     PAVs vs CNVs    5.0  1012.0   315.167416   514.0   959.0  1090.0  1129.0  1368.0
#                     SNPs vs CNVs    5.0  1565.8   841.145469   544.0  1299.0  1362.0  1783.0  2841.0
#                     SNPs vs PAVs    5.0   384.4    79.178911   265.0   376.0   391.0   404.0   486.0
#            gini     PAVs vs CNVs    5.0  2097.6   321.458862  1575.0  2123.0  2145.0  2188.0  2457.0
#                     SNPs vs CNVs    5.0  4654.0  1025.237777  2832.0  4976.0  5055.0  5113.0  5294.0
#                     SNPs vs PAVs    5.0   738.6     3.577709   735.0   736.0   738.0   740.0   744.0
