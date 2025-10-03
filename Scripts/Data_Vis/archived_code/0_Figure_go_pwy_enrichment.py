#!/usr/bin/env python3
############################################################################
# Figure 5: GO and pathway Enrichment
############################################################################

import os
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

################################  GO Enrichment ################################
# Read excel spreadsheets
snp_go = pd.read_excel("Scripts/Data_Vis/Section_4/SNP_GO_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
snp_go.columns = ["Environment", "GO", "Gene_top_has_GO", "Gene_not_top_has_GO",
                    "Gene_top_no_GO", "Gene_not_top_no_GO", "direction",
                    "p.val", "odds ratio", "qvalues"]
snp_go["Data"] = "SNP"
pav_go = pd.read_excel("Scripts/Data_Vis/Section_4/PAV_GO_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
pav_go.columns = ["Environment", "GO", "Gene_top_has_GO", "Gene_not_top_has_GO",
                    "Gene_top_no_GO", "Gene_not_top_no_GO", "direction",
                    "p.val", "odds ratio", "qvalues"]
pav_go["Data"] = "ORF"
cnv_go = pd.read_excel("Scripts/Data_Vis/Section_4/CNV_GO_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
cnv_go.columns = ["Environment", "GO", "Gene_top_has_GO", "Gene_not_top_has_GO",
                    "Gene_top_no_GO", "Gene_not_top_no_GO", "direction",
                    "p.val", "odds ratio", "qvalues"]
cnv_go["Data"] = "CNV"

# Combine
go_combined = pd.concat([snp_go, pav_go, cnv_go], axis=0)
go_combined = go_combined[go_combined["qvalues"]<=0.05] # filter for significant GO terms

# Add BP, MF, and CC information
sgd_go = pd.read_csv("Data/yeast_GO/sgd_GO_BP.tsv", sep="\t")

go_combined = go_combined.merge(sgd_go[["GO.ID", "BP", "CC", "MF"]], left_on="GO", right_on="GO.ID", how="left")
go_combined = go_combined[["Data", "Environment", "BP", "CC", "MF", "qvalues", "odds ratio", "direction"]]
go_combined["log10(odds ratio)"] = np.log10(go_combined["odds ratio"]) # log10 transform odds ratio
# note: all odds ratios are 0, so log10(0) = -Inf
# thus, will not do the below:
go_combined.loc[np.isinf(go_combined["odds ratio"]),"log10(odds ratio)"] = 1 # set infinite odds ratios to 1
go_combined = go_combined.melt(id_vars=["Data", "Environment", "log10(odds ratio)"],
                                value_vars=["BP", "CC", "MF"], var_name="GO Type",
                                value_name="GO Description") # melt GO related columns
go_combined.dropna(inplace=True) # drop rows with missing values (no GO terms)
go_combined.set_index(keys=["GO Type", "GO Description"], inplace=True) # create multi-index
go_combined = go_combined.pivot(columns=["Data", "Environment"]) # unstack multi-index
go_combined = go_combined[~(go_combined.isin([np.nan, np.NINF]).all(axis=1))] # remove rows with a mixture of NaNs and -Inf values in all columns
go_combined.to_csv("Scripts/Data_Vis/Section_4/Figure_5_data_GO_enrichment_all.txt", sep="\t")

# Create heatmaps
fig, axes = plt.subplots(nrows=go_combined.index.levels[0].nunique(),
                            ncols=go_combined.columns.levels[1].nunique(),
                            figsize=(25, 30), sharex='col', sharey='row')
cbar_ax = fig.add_axes([.91,.3,.03,.4])
for row, group in enumerate(go_combined.index.levels[0]):
    for col, category in enumerate(go_combined.columns.levels[1]):
        ax = axes[row, col]
        subset = go_combined.loc[go_combined.index.get_level_values(0)==group,
                                    go_combined.columns.get_level_values(1)==category]
        sns.heatmap(subset, ax=ax, fmt=".2f", cmap='RdBu_r', center=0,
                    vmin=-1, vmax=np.ceil(go_combined.max().max()),
                    cbar_ax=cbar_ax, xticklabels=True, yticklabels=True)
        ax.set_title(f'Group {group} - Category {category}')

plt.savefig("Scripts/Data_Vis/Section_4/Figure_5_GO_enrichment_all_heatmap_v3.pdf")
plt.close()

############################## Pathway Enrichment ##############################
# Read excel spreadsheets
snp_pwy = pd.read_excel("Scripts/Data_Vis/Section_4/SNP_pathway_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
snp_pwy.columns = ["Environment", "PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY",
                    "Gene_top_no_PWY", "Gene_not_top_no_PWY", "direction",
                    "p.val", "odds ratio", "qvalues",  "PWY name", "PWY ID"]
snp_pwy["Data"] = "SNP"
pav_pwy = pd.read_excel("Scripts/Data_Vis/Section_4/PAV_pathway_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
pav_pwy.columns = ["Environment", "PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY",
                    "Gene_top_no_PWY", "Gene_not_top_no_PWY", "direction",
                    "p.val", "odds ratio", "qvalues",  "PWY name", "PWY ID"]
pav_pwy["Data"] = "ORF"
cnv_pwy = pd.read_excel("Scripts/Data_Vis/Section_4/CNV_pathway_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
cnv_pwy.columns = ["Environment", "PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY",
                    "Gene_top_no_PWY", "Gene_not_top_no_PWY", "direction",
                    "p.val", "odds ratio", "qvalues",  "PWY name", "PWY ID"]
cnv_pwy["Data"] = "CNV"

# Combine
pwy_combined = pd.concat([snp_pwy, pav_pwy, cnv_pwy], axis=0)
pwy_combined = pwy_combined[pwy_combined["qvalues"]<=0.05] # filter for significant pathways, None were found
pwy_combined = pwy_combined[["Data", "Environment", "PWY name", "qvalues", "odds ratio", "direction"]]
pwy_combined["log10(odds ratio)"] = np.log10(pwy_combined["odds ratio"]) # log10 transform odds ratio
pwy_combined.loc[np.isinf(pwy_combined["odds ratio"]), "log10(odds ratio)"] = 1 # set infinite odds ratios to 1
pwy_combined.set_index(keys="PWY name", inplace=True) # create index
pwy_combined = pwy_combined.pivot(columns=["Data", "Environment"]) # unstack Data and Environment variables
pwy_combined = pwy_combined[~(pwy_combined.isin([np.nan, np.NINF]).all(axis=1))] # remove rows with a mixture of NaNs and -Inf values in all columns
pwy_combined = pwy_combined.loc[:, pwy_combined.columns.get_level_values(0)=="log10(odds ratio)"]
pwy_combined.to_csv("Scripts/Data_Vis/Section_4/Figure_5_data_pathway_enrichment_all.txt", sep="\t")

# Create heatmaps
# fig, axes = plt.subplots(nrows=1, ncols=pwy_combined.columns.levels[1].nunique(),
#                          figsize=(25, 30), sharex='col', sharey='row')
# cbar_ax = fig.add_axes([.91,.3,.03,.4])
# for row, group in enumerate(pwy_combined.index.unique()):
# for col, category in enumerate(pwy_combined.columns.levels[1]):
#     ax = axes[col]
#     subset = pwy_combined.loc[:, pwy_combined.columns.get_level_values(1)==category]
#     sns.heatmap(subset, fmt=".2f", cmap='RdBu_r', center=0,
#                 vmin=-1, vmax=np.ceil(pwy_combined.max().max()),
#                 cbar_ax=cbar_ax, xticklabels=True, yticklabels=True)
#     ax.set_title(f'Category {category}')

fig = plt.figure()
sns.heatmap(pwy_combined, fmt=".2f", cmap='RdBu_r', center=0,
            vmin=-1, vmax=np.ceil(pwy_combined.max().max()),
            xticklabels=True, yticklabels=True)

plt.savefig("Scripts/Data_Vis/Section_4/Figure_5_pathway_enrichment_all_heatmap_v1.pdf")
plt.close()

