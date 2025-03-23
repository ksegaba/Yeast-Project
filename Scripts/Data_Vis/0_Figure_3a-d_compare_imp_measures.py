#!/usr/bin/env python3
############################################################################
# Figure 3
############################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]

## 3a. Gini vs SHAP kernel density plots for baseline models
# Based on average values
res = {"baseline": {"snp": {}, "pav": {}, "cnv": {}},
       "optimized": {"snp": {}, "pav": {}, "cnv": {}}}
for data_type in ["snp", "pav", "cnv"]:
    gini_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_imp_{data_type}.tsv", sep="\t", index_col=0)
    shap_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_shap_{data_type}.tsv", sep="\t", index_col=0)
    
    for env in target_envs:
        df = pd.concat([gini_base_rank.loc[:, env], shap_base_rank.loc[:,env]],
            ignore_index=False, axis=1).dropna()
        df = df.abs().rank(pct=True)
        
        # Calculate spearman correlation
        r, p = spearmanr(df.iloc[:,0], df.iloc[:,1])
        res["baseline"][data_type][env] = {"r": r, "p": p}
        
        # Plotting
        sns.kdeplot(x=df.iloc[:,0], y=df.iloc[:,1], cmap="viridis", fill=True,
            cbar=True, bw_adjust=.5)
        plt.xlabel("Gini rank percentile")
        plt.ylabel("SHAP rank percentile")
        plt.title(f"{data_type}: {env}")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_3_baseline_gini_vs_shap_rank_per_{data_type}_{env}.pdf")
        plt.close()
    
    del gini_base_rank, shap_base_rank

## 3b. Gini vs SHAP kernel density plots for optimized models
for data_type in ["snp", "pav", "cnv"]:
    gini_opt_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_imp_{data_type}.tsv", sep="\t", index_col=0)
    shap_opt_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_shap_{data_type}.tsv", sep="\t", index_col=0)
    
    for env in target_envs:
        df = pd.concat([gini_opt_rank.loc[:, env], shap_opt_rank.loc[:,env]],
            ignore_index=False, axis=1).dropna()
        df = df.abs().rank(pct=True)
        
        # Calculate spearman correlation
        r, p = spearmanr(df.iloc[:,0], df.iloc[:,1])
        res["optimized"][data_type][env] = {"r": r, "p": p}
        
        # Plotting
        sns.kdeplot(x=df.iloc[:,0], y=df.iloc[:,1], cmap="viridis", fill=True,
            cbar=True, bw_adjust=.5)
        plt.xlabel("Gini rank percentile")
        plt.ylabel("SHAP rank percentile")
        plt.title(f"{data_type}: {env}")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_3_optimized_gini_vs_shap_rank_per_{data_type}_{env}.pdf")
        plt.close()
    
    del gini_opt_rank, shap_opt_rank

# Save the pandas rank percentile results
out = pd.DataFrame.from_dict({(i, j, k, h): res[i][j][k][h]
                              for i in res.keys()
                              for j in res[i].keys()
                              for k in res[i][j].keys()
                              for h in res[i][j][k].keys()},
                       orient='index')
out.index = pd.MultiIndex.from_tuples(out.index, names=["model", "variant", "env", "stat"])
out = out.pivot_table(index=["model", "variant", "env"], columns="stat", values=0)
out.to_csv("Scripts/Data_Vis/Section_4/Figure_3_gini_vs_shap_rank_per_pd_spearman.tsv",
           sep="\t")

## 3c. Genetic variant comparisons: kernel density plots for baseline models
# Based on rank percentiles
res = {"imp": {}, "shap": {}}
for imp_type in ["imp", "shap"]:
    snp_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_snp.tsv", sep="\t")
    pav_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_pav.tsv", sep="\t")
    cnv_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_cnv.tsv", sep="\t")
    
    for env in target_envs:
        ## SNP vs PAV
        snp_pav = pd.concat([snp_base_rank.loc[:,["gene", env]].groupby("gene").max(),
                            pav_base_rank.loc[:,["gene", env]].groupby("gene").max()],
                            ignore_index=False, axis=1).dropna()
        snp_pav = snp_pav.abs().rank(pct=True)
        # snp_pav = snp_pav / len(snp_pav) # Calculate the rank percentiles
        # snp_pav.columns = ["SNP", "PAV"]
        
        ## SNP vs CNV
        snp_cnv = pd.concat([snp_base_rank.loc[:,["gene", env]].groupby("gene").max(),
                            cnv_base_rank.loc[:,["gene", env]].groupby("gene").max()],
                            ignore_index=False, axis=1).dropna()
        snp_cnv = snp_cnv.abs().rank(pct=True)
        # snp_cnv = snp_cnv / len(snp_cnv)
        # snp_cnv.columns = ["SNP", "CNV"]
        
        ## PAV vs CNV
        pav_cnv = pd.concat([pav_base_rank.loc[:, ["orf", env]].set_index("orf"),
                            cnv_base_rank.loc[:, ["orf", env]].set_index("orf")],
                            ignore_index=False, axis=1).dropna()
        pav_cnv = pav_cnv.abs().rank(pct=True)
        # pav_cnv = pav_cnv / len(pav_cnv)
        # pav_cnv.columns = ["PAV", "CNV"]
        
        # print(snp_pav.sort_values(by="SNP"),
        #       snp_cnv.sort_values(by="SNP"),
        #       pav_cnv.sort_values(by="PAV"))
        
        # Calculate spearman correlation
        snp_pav_r, snp_pav_p = spearmanr(snp_pav.iloc[:,0], snp_pav.iloc[:,1])
        snp_cnv_r, snp_cnv_p = spearmanr(snp_cnv.iloc[:,0], snp_cnv.iloc[:,1])
        pav_cnv_r, pav_cnv_p = spearmanr(pav_cnv.iloc[:,0], pav_cnv.iloc[:,1])
        res[imp_type][env] = {"SNP vs PAV": {"r": snp_pav_r, "p": snp_pav_p},
                        "SNP vs CNV": {"r": snp_cnv_r, "p": snp_cnv_p},
                        "PAV vs CNV": {"r": pav_cnv_r, "p": pav_cnv_p}}
        
        # Plotting
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        sns.kdeplot(x=snp_pav.iloc[:,0], y=snp_pav.iloc[:,1], cmap="viridis",
                    fill=True, cbar=True, bw_adjust=.5, ax=ax[0])
        ax[0].set_xlabel(f"SNP {imp_type} rank percentile")
        ax[0].set_ylabel(f"PAV {imp_type} rank percentile")
        ax[0].set_title("SNP vs PAV")
        sns.kdeplot(x=snp_cnv.iloc[:,0], y=snp_cnv.iloc[:,1], cmap="viridis",
                    fill=True, cbar=True, bw_adjust=.5, ax=ax[1])
        ax[1].set_xlabel(f"SNP {imp_type} rank percentile")
        ax[1].set_ylabel(f"CNV {imp_type} rank percentile")
        ax[1].set_title("SNP vs CNV")
        sns.kdeplot(x=pav_cnv.iloc[:,0], y=pav_cnv.iloc[:,1], cmap="viridis",
                    fill=True, cbar=True, bw_adjust=.5, ax=ax[2])
        ax[2].set_xlabel(f"PAV {imp_type} rank percentile")
        ax[2].set_ylabel(f"CNV {imp_type} rank percentile")
        ax[2].set_title("PAV vs CNV")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_3_{imp_type}_rank_per_{env}_variant_comparison_rank_per_pd.pdf") # rank per pandas
        # plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_3_{imp_type}_rank_per_{env}_variant_comparison.pdf") # rank
        # plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_3_{imp_type}_rank_per_{env}_variant_comparison_rank_per.pdf") # rank per manual
        plt.close()
    
    del snp_base_rank, pav_base_rank, cnv_base_rank

# Save the pandas rank percentile results
out = pd.DataFrame.from_dict({(i, j, k, h): res[i][j][k][h]
                              for i in res.keys()
                              for j in res[i].keys()
                              for k in res[i][j].keys()
                              for h in res[i][j][k].keys()},
                       orient='index')
out.index = pd.MultiIndex.from_tuples(out.index, names=["imp_type", "env", "comparison", "stat"])
out = out.pivot_table(index=["imp_type", "env", "comparison"], columns="stat", values=0)
out.rename(index={"imp": "gini"}, inplace=True)
out.to_csv("Scripts/Data_Vis/Section_4/Table_S5_figure_3_baseline_rf_rank_per_pd_variant_comparison_spearman.tsv",
           sep="\t")

## 3c. Genetic variant comparisons: kernel density plots for optimized models
# Based on rank percentiles
res = {"imp": {}, "shap": {}}
for imp_type in ["imp", "shap"]:
    snp_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_snp.tsv", sep="\t")
    pav_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_pav.tsv", sep="\t")
    cnv_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_cnv.tsv", sep="\t")
    
    for env in target_envs:
        ## SNP vs PAV
        snp_pav = pd.concat([snp_base_rank.loc[:,["gene", env]].groupby("gene").max(),
                            pav_base_rank.loc[:,["gene", env]].groupby("gene").max()],
                            ignore_index=False, axis=1).dropna()
        snp_pav = snp_pav.abs().rank(pct=True)
        
        ## SNP vs CNV
        snp_cnv = pd.concat([snp_base_rank.loc[:,["gene", env]].groupby("gene").max(),
                            cnv_base_rank.loc[:,["gene", env]].groupby("gene").max()],
                            ignore_index=False, axis=1).dropna()
        snp_cnv = snp_cnv.abs().rank(pct=True)
        
        ## PAV vs CNV
        pav_cnv = pd.concat([pav_base_rank.loc[:, ["orf", env]].set_index("orf"),
                            cnv_base_rank.loc[:, ["orf", env]].set_index("orf")],
                            ignore_index=False, axis=1).dropna()
        pav_cnv = pav_cnv.abs().rank(pct=True)
        
        # Calculate spearman correlation
        snp_pav_r, snp_pav_p = spearmanr(snp_pav.iloc[:,0], snp_pav.iloc[:,1])
        snp_cnv_r, snp_cnv_p = spearmanr(snp_cnv.iloc[:,0], snp_cnv.iloc[:,1])
        pav_cnv_r, pav_cnv_p = spearmanr(pav_cnv.iloc[:,0], pav_cnv.iloc[:,1])
        res[imp_type][env] = {"SNP vs PAV": {"r": snp_pav_r, "p": snp_pav_p, "n_genes": len(snp_pav)},
                        "SNP vs CNV": {"r": snp_cnv_r, "p": snp_cnv_p, "n_genes": len(snp_cnv)},
                        "PAV vs CNV": {"r": pav_cnv_r, "p": pav_cnv_p, "n_genes": len(pav_cnv)}}
        
        # Plotting
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        sns.kdeplot(x=snp_pav.iloc[:,0], y=snp_pav.iloc[:,1], cmap="viridis",
                    fill=True, cbar=True, bw_adjust=.5, ax=ax[0])
        ax[0].set_xlabel(f"SNP {imp_type} rank percentile")
        ax[0].set_ylabel(f"PAV {imp_type} rank percentile")
        ax[0].set_title("SNP vs PAV")
        sns.kdeplot(x=snp_cnv.iloc[:,0], y=snp_cnv.iloc[:,1], cmap="viridis",
                    fill=True, cbar=True, bw_adjust=.5, ax=ax[1])
        ax[1].set_xlabel(f"SNP {imp_type} rank percentile")
        ax[1].set_ylabel(f"CNV {imp_type} rank percentile")
        ax[1].set_title("SNP vs CNV")
        sns.kdeplot(x=pav_cnv.iloc[:,0], y=pav_cnv.iloc[:,1], cmap="viridis",
                    fill=True, cbar=True, bw_adjust=.5, ax=ax[2])
        ax[2].set_xlabel(f"PAV {imp_type} rank percentile")
        ax[2].set_ylabel(f"CNV {imp_type} rank percentile")
        ax[2].set_title("PAV vs CNV")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_3_{imp_type}_optimized_rank_per_{env}_variant_comparison_rank_per_pd.pdf") # rank per pandas
        plt.close()
    
    del snp_base_rank, pav_base_rank, cnv_base_rank

# Save the pandas rank percentile results
out = pd.DataFrame.from_dict({(i, j, k, h): res[i][j][k][h]
                              for i in res.keys()
                              for j in res[i].keys()
                              for k in res[i][j].keys()
                              for h in res[i][j][k].keys()},
                       orient='index')
out.index = pd.MultiIndex.from_tuples(out.index, names=["imp_type", "env", "comparison", "stat"])
out = out.pivot_table(index=["imp_type", "env", "comparison"], columns="stat", values=0)
out.rename(index={"imp": "gini"}, inplace=True)
out.to_csv("Scripts/Data_Vis/Section_4/Table_S5_figure_3_optimized_rf_rank_per_pd_variant_comparison_spearman.tsv",
           sep="\t")


### 3A-C Heatmaps
ab_data = pd.read_csv("Scripts/Data_Vis/Section_4/Figure_3_gini_vs_shap_rank_per_pd_spearman.tsv", sep="\t")

# 3A
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
sns.heatmap(ab_data.loc[ab_data.model=="baseline", :].pivot_table(index="variant", columns="env", values="r"),
            cmap="RdBu_r", annot=True, annot_kws={"size": 7}, fmt=".2f", ax=ax[0], cbar=True, square=True)
sns.heatmap(ab_data.loc[ab_data.model=="baseline", :].pivot_table(index="variant", columns="env", values="p"),
            cmap="Reds", annot=True, annot_kws={"size": 7}, fmt=".2e", ax=ax[1], cbar=True, square=True)
ax[0].set_title("Gini vs SHAP Rank Percentiles Baseline Models")
ax[1].set_title("Gini vs SHAP Rank Percentiles Baseline Models")
cbar = ax[0].collections[0].colorbar
cbar.set_label("rho")
cbar = ax[1].collections[0].colorbar
cbar.set_label("p-value")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Figure_3a_baseline_gini_vs_shap_rank_per_pd_spearman.pdf")
plt.close()

# 3B
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
sns.heatmap(ab_data.loc[ab_data.model=="optimized", :].pivot_table(index="variant", columns="env", values="r"),
            cmap="RdBu_r", annot=True, annot_kws={"size": 7}, fmt=".2f", ax=ax[0], cbar=True, square=True)
sns.heatmap(ab_data.loc[ab_data.model=="optimized", :].pivot_table(index="variant", columns="env", values="p"),
            cmap="Reds", annot=True, annot_kws={"size": 7}, fmt=".2e", ax=ax[1], cbar=True, square=True)
ax[0].set_title("Gini vs SHAP Rank Percentiles Optimized Models")
ax[1].set_title("Gini vs SHAP Rank Percentiles Optimized Models")
cbar = ax[0].collections[0].colorbar
cbar.set_label("rho")
cbar = ax[1].collections[0].colorbar
cbar.set_label("p-value")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Figure_3b_optimized_gini_vs_shap_rank_per_pd_spearman.pdf")
plt.close()

# 3C
c_data = pd.read_csv("Scripts/Data_Vis/Section_4/Table_S5_figure_3_baseline_rf_rank_per_pd_variant_comparison_spearman.tsv", sep="\t")
fig, ax = plt.subplots(2, 2, figsize=(10, 10))
sns.heatmap(c_data.loc[c_data.imp_type=="gini", :].pivot_table(index="comparison", columns="env", values="r"),
            cmap="RdBu_r", annot=True, annot_kws={"size": 7}, fmt=".2f", ax=ax[0][0], cbar=True, square=True)
sns.heatmap(c_data.loc[c_data.imp_type=="gini", :].pivot_table(index="comparison", columns="env", values="p"),
            cmap="Reds", annot=True, annot_kws={"size": 7}, fmt=".2e", ax=ax[0][1], cbar=True, square=True)
sns.heatmap(c_data.loc[c_data.imp_type=="shap", :].pivot_table(index="comparison", columns="env", values="r"),
            cmap="RdBu_r", annot=True, annot_kws={"size": 7}, fmt=".2f", ax=ax[1][0], cbar=True, square=True)
sns.heatmap(c_data.loc[c_data.imp_type=="shap", :].pivot_table(index="comparison", columns="env", values="p"),
            cmap="Reds", annot=True, annot_kws={"size": 7}, fmt=".2e", ax=ax[1][1], cbar=True, square=True)
ax[0][0].set_title("Gini, Variant comparison, Baseline Models")
ax[1][0].set_title("SHAP, Variant comparison, Baseline Models")
cbar = ax[0][0].collections[0].colorbar
cbar.set_label("rho")
cbar = ax[0][1].collections[0].colorbar
cbar.set_label("p-value")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Figure_3c_baseline_rf_rank_per_pd_variant_comparison_spearman.pdf")
plt.close()

c_data = pd.read_csv("Scripts/Data_Vis/Section_4/Table_S5_figure_3_optimized_rf_rank_per_pd_variant_comparison_spearman.tsv", sep="\t")
fig, ax = plt.subplots(2, 2, figsize=(10, 10))
sns.heatmap(c_data.loc[c_data.imp_type=="gini", :].pivot_table(index="comparison", columns="env", values="r"),
            cmap="RdBu_r", annot=True, annot_kws={"size": 7}, fmt=".2f", ax=ax[0][0], cbar=True, square=True)
sns.heatmap(c_data.loc[c_data.imp_type=="gini", :].pivot_table(index="comparison", columns="env", values="p"),
            cmap="Reds", annot=True, annot_kws={"size": 7}, fmt=".2e", ax=ax[0][1], cbar=True, square=True)
sns.heatmap(c_data.loc[c_data.imp_type=="shap", :].pivot_table(index="comparison", columns="env", values="r"),
            cmap="RdBu_r", annot=True, annot_kws={"size": 7}, fmt=".2f", ax=ax[1][0], cbar=True, square=True)
sns.heatmap(c_data.loc[c_data.imp_type=="shap", :].pivot_table(index="comparison", columns="env", values="p"),
            cmap="Reds", annot=True, annot_kws={"size": 7}, fmt=".2e", ax=ax[1][1], cbar=True, square=True)
ax[0][0].set_title("Gini, Variant comparison, Optimized Models")
ax[1][0].set_title("SHAP, Variant comparison, Optimized Models")
cbar = ax[0][0].collections[0].colorbar
cbar.set_label("rho")
cbar = ax[0][1].collections[0].colorbar
cbar.set_label("p-value")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Figure_3c_optimized_rf_rank_per_pd_variant_comparison_spearman.pdf")
plt.close()