#!/usr/bin/env python3
############################################################################
# Figure 2
############################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]

## Gini vs SHAP kernel density plots for baseline models
for data_type in ["snp", "pav", "cnv"]:
    gini_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_imp_{data_type}_rank_per.tsv")
    shap_base_rank = pd.read_csv(f"Scripts/Data_Vis/Section_4/RF_baseline_shap_{data_type}_rank_per.tsv")
    for env in target_envs:
        df = pd.concat([gini_base_rank.loc[:,env], shap_base_rank.loc[:,env]],
            ignore_index=False, axis=1).dropna()
        sns.kdeplot(x=df.iloc[:,0], y=df.iloc[:,1], cmap="viridis", fill=True,
            cbar=True, bw_adjust=.5)
        plt.xlabel("Gini rank percentile")
        plt.ylabel("SHAP rank percentile")
        plt.title(f"{data_type}: {env}")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_4/Figure_4_gini_vs_shap_rank_per_{data_type}_{env}.pdf")
        plt.close()
    del gini_base_rank, shap_base_rank

## Genetic variant comparisons: kernel density plots for baseline models
