#!/usr/bin/env python3
################################################################################
# Table S6 (needed for Figure 3)
################################################################################

from scipy.stats import ks_2samp
from tqdm import tqdm
import os
import re
import pandas as pd
import numpy as np
import datatable as dt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
from scipy.stats import spearmanr, linregress, pearsonr
from matplotlib.colors import ListedColormap, BoundaryNorm

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

map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                       sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")

################################################################################
# 3e. ORF (gene) importance vs percent presence density plots
# Read presence/absence variation data and calculate percent orf presence in the population
pav = pd.read_csv("Data/Peter_2018/ORFs_pres_abs.csv", index_col=0)
percent_presence = pd.DataFrame(
    (pav.sum()/pav.shape[0])*100, index=pav.columns).reset_index()
percent_presence.columns = ["orf", "percent_presence"]
percent_presence.set_index("orf", inplace=True)
percent_presence.percent_presence = percent_presence.percent_presence.astype(
    "float64")
percent_presence.index = percent_presence.apply(
    lambda x: re.sub("^X", "", x.name), axis=1)
percent_presence.index = percent_presence.apply(
    lambda x: re.sub("\.", "-", x.name), axis=1)

# Read in orf (gene) copy number data
cnv = pd.read_csv(
    "~/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv", index_col=0)
cnv_stats = cnv.describe().T
cnv_stats.index = cnv_stats.apply(lambda x: re.sub("^X", "", x.name), axis=1)
cnv_stats.index = cnv_stats.apply(lambda x: re.sub("\.", "-", x.name), axis=1)

# Collect results for ORF gene importance vs percent presence
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]
res = {"baseline": {"imp": {"snp": {}, "pav": {}, "cnv": {}},
                    "shap": {"snp": {}, "pav": {}, "cnv": {}}},
       "optimized": {"imp": {"snp": {}, "pav": {}, "cnv": {}},
                     "shap": {"snp": {}, "pav": {}, "cnv": {}}}}

presence_vs_cnv_res = {"baseline": {"imp": {"snp": {}, "pav": {}, "cnv": {}},
                                    "shap": {"snp": {}, "pav": {}, "cnv": {}}},
                       "optimized": {"imp": {"snp": {}, "pav": {}, "cnv": {}},
                                     "shap": {"snp": {}, "pav": {}, "cnv": {}}}}

for imp_type in ["imp", "shap"]:
    for data_type in ["pav", "cnv"]:
        if data_type == "snp":  # baseline model feature importances
            baseline = pd.read_csv(
                f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_snp.tsv",
                sep="\t", index_col=0)
        else:
            baseline = pd.read_csv(
                f"Scripts/Data_Vis/Section_4/RF_baseline_{imp_type}_{data_type}.tsv",
                sep="\t", index_col=0)
        opt = pd.read_csv(
            f"Scripts/Data_Vis/Section_4/RF_FS_{imp_type}_{data_type}.tsv",
            sep="\t", index_col=0)  # optimized model feature importances
        #
        # Plot gene importance vs percent presence scatter plots
        # colored by average copy number
        fig, axes = plt.subplots(5, 2, figsize=(5.6, 11))
        # fig2, axes2 = plt.subplots(5, 2, figsize=(5.6, 11)) # density plots
        for i, env in enumerate(target_envs):
            # Plot baseline model performances
            # colored by average copy number
            toplot = percent_presence.merge(baseline[env], how="inner",
                                            left_index=True, right_index=True)
            toplot = toplot.merge(cnv_stats, how="left",
                                  left_index=True, right_index=True)
            #
            # Calculate spearman's rho and p-value between % presence and importance
            r, p = spearmanr(toplot["percent_presence"], toplot[env])
            #
            # Calculate pearson's r between % presence and copy number
            pcc, pcc_p = pearsonr(toplot["percent_presence"], toplot["mean"])
            presence_vs_cnv_res["baseline"][imp_type][data_type][env] = {
                "pcc": pcc, "pcc p-val": pcc_p}
            #
            sns.scatterplot(data=toplot, x="percent_presence", y=env,
                            c=cm.plasma(toplot["mean"]), alpha=0.3, ax=axes[i, 0])

            m, b, r2, p_val, std_err = linregress(
                toplot["percent_presence"], toplot[env])  # linear fit
            axes[i, 0].plot(toplot["percent_presence"], m *
                            toplot["percent_presence"] + b)  # plot line
            res["baseline"][imp_type][data_type][env] = {"rho": r, "rho p-val": p, "m": m,
                                                         "b": b, "r2": r, "p-val": p_val, "std_err": std_err}  # save results

            # normalize mean gene copy number values
            norm = plt.Normalize(toplot["mean"].min(), toplot["mean"].max())
            cbar = plt.cm.ScalarMappable(norm=norm, cmap="plasma")  # colorbar
            axes[i, 0].figure.colorbar(cbar, ax=axes[i, 0])  # set colorbar
            axes[i, 0].set_title(f"Baseline: {env}")
            #
            # density plot
            # sns.kdeplot(x=toplot["percent_presence"], y=toplot[env],
            # 	cmap="viridis", fill=True, cbar=True, bw_adjust=.5, ax=axes2[i,0])
            # axes2[i,0].set_title(f"Baseline: {env}")
            del toplot
            #
            # Plot optimized model performances
            # colored by average copy number
            toplot = percent_presence.merge(opt[env].dropna(), how="inner",
                                            left_index=True, right_index=True)
            toplot = toplot.merge(cnv_stats, how="left",
                                  left_index=True, right_index=True)
            #
            # Calculate spearman's rho and p-value
            r, p = spearmanr(toplot["percent_presence"], toplot[env])
            #
            # Calculate pearson's r between % presence and copy number
            pcc, pcc_p = pearsonr(toplot["percent_presence"], toplot["mean"])
            presence_vs_cnv_res["optimized"][imp_type][data_type][env] = {
                "pcc": pcc, "pcc p-val": pcc_p}
            #
            sns.scatterplot(data=toplot, x="percent_presence", y=env,
                            c=cm.plasma(toplot["mean"]), alpha=0.3, ax=axes[i, 1])

            m, b, r2, p_val, std_err = linregress(
                toplot["percent_presence"], toplot[env])  # linear fit
            axes[i, 1].plot(toplot["percent_presence"], m *
                            toplot["percent_presence"] + b)  # plot line
            res["optimized"][imp_type][data_type][env] = {"rho": r, "rho p-val": p, "m": m,
                                                          "b": b, "r2": r, "p-val": p_val, "std_err": std_err}  # save results

            # normalize mean gene copy number values
            norm = plt.Normalize(toplot["mean"].min(), toplot["mean"].max())
            cbar = plt.cm.ScalarMappable(norm=norm, cmap="plasma")  # colorbar
            axes[i, 1].figure.colorbar(cbar, ax=axes[i, 1])  # set colorbar
            axes[i, 1].set_title(f"Optimized: {env}")
            #
            # density plot
            # sns.kdeplot(x=toplot["percent_presence"], y=toplot[env],
            # 	cmap="viridis", fill=True, cbar=True, bw_adjust=.5, ax=axes2[i,1])
            # axes2[i,1].set_title(f"Optimized: {env}")
            del toplot

        plt.tight_layout()  # fig
        plt.savefig(
            f'Scripts/Data_Vis/Section_4/Figure_3_%_presence_vs_{imp_type}_{data_type}.pdf')  # fig
        # plt.savefig(f'Scripts/Data_Vis/Section_4/Figure_3_%_presence_vs_{imp_type}_{data_type}_density.pdf') # fig2
        plt.close()

out = pd.DataFrame.from_dict({(i, j, k, h): res[i][j][k][h]
                              for i in res.keys()
                              for j in res[i].keys()
                              for k in res[i][j].keys()
                              for h in res[i][j][k].keys()},
                             orient="index")
out.index.names = ["model", "imp_type", "data_type", "env"]
out.to_csv("Scripts/Data_Vis/Section_4/Table_S6_figure_3_%_presence_vs_gene_importance_spearman.tsv", sep="\t")

out2 = pd.DataFrame.from_dict({(i, j, k, h): presence_vs_cnv_res[i][j][k][h]
                               for i in presence_vs_cnv_res.keys()
                               for j in presence_vs_cnv_res[i].keys()
                               for k in presence_vs_cnv_res[i][j].keys()
                               for h in presence_vs_cnv_res[i][j][k].keys()},
                              orient="index")
out2.index.names = ["model", "imp_type", "data_type", "env"]
out2.to_csv(
    "Scripts/Data_Vis/Section_4/Table_S6_figure_3_%_presence_vs_average_gene_cnv_pearson.tsv", sep="\t")

# old code (makes baseline and optimized model density SCATTER plots, not kde plots):
# Plot density scatter plots of percent presence vs gene importance
# res = []
# fig, axes = plt.subplots(5, 2, figsize=(5.6, 11))
# for i, col in enumerate(imp.columns):
# 	# Density scatter plots
# 	spearman_corr, p = spearmanr(percent_presence["percent_presence"], imp[col]) # spearman's rho and p-value
# 	values = np.vstack([percent_presence["percent_presence"], imp[col]])
# 	kernel = stats.gaussian_kde(values)(values) # point kernel density
# 	row_idx, col_idx = divmod(i, 5)  # row and column index for subplots
# 	ax = axes[row_idx, col_idx] # set subplot
# 	sns.scatterplot(x=percent_presence["percent_presence"], y=imp[col],
# 					c=kernel, edgecolor = 'none', size=-5, alpha=.3,
# 					cmap="viridis", legend=False, ax=ax) # density scatter plot
# 	ax.set_title(f'{mapping[col]}')
# 	norm = plt.Normalize(kernel.min(), kernel.max()) # normalize kde values
# 	cbar = plt.cm.ScalarMappable(norm=norm, cmap="viridis") # colorbar
# 	ax.figure.colorbar(cbar, ax=ax) # set colorbar
# 	m, b, r, p_val, std_err = linregress(percent_presence["percent_presence"], imp[col]) # linear fit
# 	ax.plot(percent_presence["percent_presence"], m*percent_presence["percent_presence"] + b) # plot line
# 	res.append([col, spearman_corr, p, m, b, r, p_val, std_err]) # save stats results
# 	del spearman_corr, p, values, kernel, row_idx, col_idx, ax, norm, cbar, m, b, r, p_val, std_err
# plt.tight_layout()
# plt.savefig("Scripts/Data_Vis/Figure2c_baseline.pdf")
# plt.close()
# res = pd.DataFrame(res, columns=["Environment", "Spearman's rho", "p-value",
# 									"lm slope", "lm intercept", "lm pearson r",
# 									"lm p-value", "lm std_err"])
# res.to_csv("Scripts/Data_Vis/Table_S_baseline_importance_vs_percent_presence.csv", index=False)
#
# # Read in orf importance scores data for RF FS models
# dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs"
# pav_rf_res = pd.read_csv("Results/RESULTS_RF_ORFs_FS.txt", sep="\t")  # ORF pres/abs FS results
# pav_files = [os.path.join(dir, f"{x}_imp") for x in pav_rf_res['ID']]
#
# # Plot density scatter plots
# res = []
# fig, axes = plt.subplots(7, 5, figsize=(15, 15))
# fig2, axes2 = plt.subplots(7, 5, figsize=(15, 15))
# for i, col in enumerate(mapping.keys()):
# 	file = [f for f in pav_files if col in f] # Read feature importance score files
# 	imp = pd.read_csv(file[0], sep="\t", skiprows=1, names=["orf", "mean_imp"], index_col=0)
# 	imp = imp.merge(percent_presence, how="left", on="orf") # merge with percent_presence data
# 	imp = imp.merge(cnv_stats, how="left", left_on="orf", right_index=True) # merge with cnv_stats data
# 	spearman_corr, p = spearmanr(imp["percent_presence"], imp["mean_imp"]) # spearman's rho and p-value
# 	values = np.vstack([imp["percent_presence"], imp["mean_imp"]])
# 	kernel = stats.gaussian_kde(values)(values) # point kernel density
# 	row_idx, col_idx = divmod(i, 5)  # row and column index for subplots
# 	ax = axes[row_idx, col_idx] # set subplot
# 	ax2 = axes2[row_idx, col_idx] # set subplot
# 	# color by point density
# 	sns.scatterplot(x=imp["percent_presence"], y=imp["mean_imp"],
# 					c=kernel, edgecolor = 'none', size=-5, alpha=.3,
# 					cmap="viridis", legend=False, ax=ax) # density scatter plot
# 	ax.set_title(f"{mapping[col]}")
# 	norm = plt.Normalize(kernel.min(), kernel.max()) # normalize kde values
# 	cbar = plt.cm.ScalarMappable(norm=norm, cmap="viridis") # colorbar
# 	ax.figure.colorbar(cbar, ax=ax) # set colorbar
# 	m, b, r, p_val, std_err = linregress(imp["percent_presence"], imp["mean_imp"]) # linear fit
# 	ax.plot(imp["percent_presence"], m*imp["percent_presence"] + b) # plot line
# 	res.append([col, spearman_corr, p, m, b, r, p_val, std_err]) # save stats results
# 	# sizes by gene copy number
# 	sns.scatterplot(x=imp["percent_presence"], y=imp["mean_imp"], c=np.log10(imp["mean"]),
# 					size=np.log10(imp["mean"]), edgecolor = 'gray', alpha=.15,
# 					cmap="plasma", legend=True, ax=ax2) # density scatter plot
# 	ax2.set_title(f"{mapping[col]}")
# 	norm = plt.Normalize(np.log10(imp["mean"]).min(), np.log10(imp["mean"]).max()) # normalize gene copy number values
# 	cbar = plt.cm.ScalarMappable(norm=norm, cmap="plasma") # colorbar
# 	ax2.figure.colorbar(cbar, ax=ax2) # set colorbar
# 	ax2.plot(imp["percent_presence"], m*imp["percent_presence"] + b) # plot line
# 	del spearman_corr, p, values, kernel, row_idx, col_idx, ax, norm, cbar, m, b, r, p_val, std_err
# fig.tight_layout()
# fig.savefig("Scripts/Data_Vis/Figure2c_fs.pdf")
# fig2.tight_layout()
# fig2.savefig(f"Scripts/Data_Vis/Figure2d_fs.pdf")
# plt.close(); plt.close()
# res = pd.DataFrame(res, columns=["Environment", "Spearman's rho", "p-value",
# 									"lm slope", "lm intercept", "lm pearson r",
# 									"lm p-value", "lm std_err"])
# res.to_csv("Scripts/Data_Vis/Table_S_fs_importance_vs_percent_presence.csv", index=False)

################################################################################
# 3# . SNP frequency vs importance
snp = dt.fread('Data/Peter_2018/geno.csv').to_pandas()

# ... need to finish


################################################################################
# 3f. Heatmap of spearman correlation between envs
f_data = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Table_S7_rank_per_corr_btwn_envs.tsv", sep="\t")
f_data.set_index(["Env1", "Env2"], inplace=True)
count_data = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Table_S7_rank_per_corr_btwn_envs_shared_feat_counts.tsv", sep="\t")
count_data.set_index(["Unnamed: 0", "Unnamed: 1"], inplace=True)
count_data.index.names = ["Env1", "Env2"]

env_order = ['YPDCHX05', 'YPDCHX1', 'YPDANISO50', 'YPDANISO10', 'YPDANISO20',
                         'YPDDMSO', 'YPDMV', 'YPDSDS', 'YPD40', 'YPD42', 'YPDKCL2M',
                         'YPDCAFEIN40', 'YPDCAFEIN50', 'YPDBENOMYL200', 'YPDBENOMYL500',
                         'YPDETOH', 'YPDNYSTATIN', 'YPACETATE', 'YPXYLOSE', 'YPRIBOSE',
                         'YPSORBITOL', 'YPGLYCEROL', 'YPETHANOL', 'YPGALACTOSE',
                         'YPDLICL250MM',  'YPDNACL15M', 'YPDNACL1M', 'YPDFORMAMIDE4',
                         'YPDFORMAMIDE5', 'YPDHU', 'YPD14', 'YPDFLUCONAZOLE', 'YPDSODIUMMETAARSENITE',
                         'YPD6AU', 'YPDCUSO410MM']  # According to Fig. 1A


def df2triangle(df, upper=True):
    """Convert a dataframe to an upper triangle matrix or a lower triangle matrix"""
    if df.shape[0] != df.shape[1]:
        raise ValueError("Dataframe must be square")

    if upper:
        df_upper = df.copy(deep=True)
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                if i == j:
                    df_upper.iloc[i, j] = np.nan
                elif i < j:
                    if pd.isna(df.iloc[i, j]):
                        df_upper.iloc[i, j] = df.iloc[j, i]
                    else:
                        df_upper.iloc[i, j] = df.iloc[i, j]
                elif i > j:
                    df_upper.iloc[i, j] = np.nan
        return df_upper

    else:
        df_lower = df.copy(deep=True)
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                if i == j:
                    df_lower.iloc[i, j] = np.nan
                elif i > j:
                    if pd.isna(df.iloc[i, j]):
                        df_lower.iloc[i, j] = df.iloc[j, i]
                    else:
                        df_lower.iloc[i, j] = df.iloc[i, j]
                elif i < j:
                    df_lower.iloc[i, j] = np.nan

        return df_lower


def add_missing_row_col(rho_dat):
    """Add missing row and column to a 34x34 square dataframe formed from unstack()"""
    df = rho_dat.copy(deep=True)
    missing_in_rows = [e for e in env_order if e not in df.index.values][0]
    missing_in_cols = [e for e in env_order if e not in df.columns.values][0]
    df.insert(0, missing_in_cols, np.nan)
    df = pd.concat([df, pd.Series(
        [np.nan]*35, name=missing_in_rows).to_frame().T], axis=0).iloc[:, :35]
    return df


for data_type in ["snp", "pav", "cnv"]:
    toplot = f_data.loc[:, f_data.columns.str.contains(data_type)]

    for i in range(0, 8, 2):
        # 34 by 34, missing one env
        rho_dat = add_missing_row_col(toplot.iloc[:, i].unstack())
        # 34 by 34, missing one env
        p_dat = add_missing_row_col(toplot.iloc[:, i+1].unstack())
        c_dat = add_missing_row_col(
            count_data[toplot.columns[i].replace('rho', 'shared_feats')].unstack())

        # sort according to Fig 1A
        rho_dat = df2triangle(rho_dat.loc[env_order, env_order])
        p_dat = df2triangle(p_dat.loc[env_order, env_order])
        c_dat = df2triangle(c_dat.loc[env_order, env_order])
        c_dat = c_dat.fillna(0).astype(int)

        # plot heatmaps (rho and number of shared features)
        fig, ax = plt.subplots(1, 2, figsize=(30, 15))
        f = sns.heatmap(rho_dat, ax=ax[0],
                        cmap="RdBu_r", annot=True, annot_kws={"size": 6}, center=0,
                        fmt=".2f", cbar=True, square=True, vmin=-1, vmax=1)  # ordered according to Fig 1A
        ax[0].set_title(toplot.columns[i])

        # sns.heatmap(df2triangle(p_dat.loc[env_order, env_order], upper=False),
        # 			ax=ax[1], cmap="Reds", annot=False, cbar=True, square=True) # ordered according to Fig 1A
        sns.heatmap(c_dat, ax=ax[1], cmap="YlGnBu", annot=True,
                    annot_kws={"size": 6}, fmt="d", cbar=True, square=True)
        ax[1].set_title(toplot.columns[i].replace('rho', 'shared_feats'))
        plt.tight_layout()
        plt.savefig(
            f"Scripts/Data_Vis/Section_4/Figure_S5_{toplot.columns[i]}_rank_per_corr_btwn_envs.pdf")
        plt.close()

    del toplot


################################################################################
# 3f. distribution of how many genes are shared in n number of environments


def make_fig3f(df, data_type, imp_type, ax, ax_ridx, ax_cidx):
    """Make Figure 3f: distribution of how many genes are shared in n number of
    environments"""
    #
    # Randomize environment columns individually
    n = 10000  # number of iterations to randomize counts
    # empty dictionary to store new env counts for each gene after randomization
    counts_rand = {str(i): None for i in range(n)}
    ks_res = {str(i): None for i in range(n)}  # KS test results
    #
    # original number of environments per gene
    orig_env_counts = df.sum(axis=1)
    #
    for i in tqdm(range(n)):
        df_rand = {}
        for j, env in enumerate(df.columns):
            # randomized sample of if gene is a top gene or not (1/0)
            df_rand[env] = np.random.default_rng(
                seed=i * 100 + j).permutation(df[env])
        #
        df_rand = pd.DataFrame(df_rand)
        # Recalculate number of environments a gene is a top feature for
        env_counts = df_rand.sum(axis=1)
        counts_rand[str(i)] = env_counts  # Save results
        # 2-sample KS test (orig_env_counts is greater than random chance?)
        ks_res[str(i)] = ks_2samp(orig_env_counts,
                                  env_counts, alternative="greater")
        #
    # Plot KS test results
    ks_df = pd.DataFrame.from_dict(ks_res, orient="index")
    f = ks_df.hist(figsize=(8.5, 4), bins=50)
    f[0][1].set_xlim(0, 0.05)
    plt.savefig(
        f"Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_{n}permutations_{data_type}_{imp_type}_CORRECTED.pdf")
    plt.close()
    ks_df.to_csv(
        f"Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_{n}permutations_{data_type}_{imp_type}_CORRECTED.csv")
    del f
    #
    # Plot histograms of counts
    # Calculate the median number of environments across all randomizations and the original counts
    # get frequency of env counts at each randomization iteration
    counts_rand_freq = pd.DataFrame(counts_rand).apply(pd.value_counts)
    counts_rand_freq_median = counts_rand_freq.median(
        axis=1).astype("int")  # calculate median counts for each bin
    # KS test between original and median randomized counts
    ks_test = ks_2samp(
        orig_env_counts, counts_rand_freq_median, alternative="greater")
    print(f"KS test between original and median randomized counts: {ks_test}")
    # counted environments per gene
    counted_envs = df.apply(lambda x: " // ".join(x.index), axis=1)
    #
    # Plot the distribution of median counts for each random iteration
    orig_env_counts.plot.hist(bins=35, range=[1, 35], alpha=0.5, color="blue",
                              label=imp_type, ax=ax[ax_ridx][ax_cidx])
    pd.concat([orig_env_counts, counted_envs], ignore_index=False, axis=1).to_csv(
        f"Scripts/Data_Vis/Figure_3f_{data_type}_{imp_type}_data_CORRECTED.csv")
    ax[ax_ridx][ax_cidx].axvline(orig_env_counts.median(
    ), ls="-", color="blue", linewidth=1)  # median env count
    counts_rand_freq_median.plot.line(alpha=0.5, color="black",
                                      label="median randomized counts",
                                      ax=ax[ax_ridx][ax_cidx])  # randomized counts
    ax[ax_ridx][ax_cidx].set_title(f"{data_type}: {imp_type}")
    ax[ax_ridx][ax_cidx].legend()
    #
    return ks_df, ax


# Calculate the Gini feature importance rank percentiles
# snp_gini_rank = pd.read_csv("Scripts/Data_Vis/Section_4/RF_baseline_imp_snp_rank_per.tsv", sep="\t", index_col=0)
# pav_gini_rank = pd.read_csv("Scripts/Data_Vis/Section_4/RF_baseline_imp_pav_rank_per.tsv", sep="\t", index_col=0)
# cnv_gini_rank = pd.read_csv("Scripts/Data_Vis/Section_4/RF_baseline_imp_cnv_rank_per.tsv", sep="\t", index_col=0)
snp_gini = pd.read_csv(
    "Scripts/Data_Vis/Section_4/RF_FS_imp_snp.tsv", sep="\t", index_col=0)
pav_gini = pd.read_csv(
    "Scripts/Data_Vis/Section_4/RF_FS_imp_pav.tsv", sep="\t", index_col=0)
cnv_gini = pd.read_csv(
    "Scripts/Data_Vis/Section_4/RF_FS_imp_cnv.tsv", sep="\t", index_col=0)

snp_gini_env_rank_out = pd.DataFrame()
pav_gini_env_rank_out = pd.DataFrame()
cnv_gini_env_rank_out = pd.DataFrame()
for env in mapping.keys():
    snp_gini_env = snp_gini[[env]].dropna()
    pav_gini_env = pav_gini[[env]].dropna()
    cnv_gini_env = cnv_gini[[env]].dropna()
    #
    # drop features with 0 importance
    snp_gini_env = snp_gini_env[snp_gini_env[env] > 0]
    pav_gini_env = pav_gini_env[pav_gini_env[env] > 0]
    cnv_gini_env = cnv_gini_env[cnv_gini_env[env] > 0]
    #
    # calculate the ranks
    snp_gini_env_rank = snp_gini_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    pav_gini_env_rank = pav_gini_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    cnv_gini_env_rank = cnv_gini_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    #
    # concatenate all env ranks
    snp_gini_env_rank_out = pd.concat(
        [snp_gini_env_rank_out, snp_gini_env_rank], axis=1)
    pav_gini_env_rank_out = pd.concat(
        [pav_gini_env_rank_out, pav_gini_env_rank], axis=1)
    cnv_gini_env_rank_out = pd.concat(
        [cnv_gini_env_rank_out, cnv_gini_env_rank], axis=1)

# Calculate the average absolute SHAP feature importance rank percentiles
# snp_shap_rank = pd.read_csv("Scripts/Data_Vis/Section_4/RF_baseline_shap_snp_rank_per.tsv", sep="\t", index_col=0)
# pav_shap_rank = pd.read_csv("Scripts/Data_Vis/Section_4/RF_baseline_shap_pav_rank_per.tsv", sep="\t", index_col=0)
# cnv_shap_rank = pd.read_csv("Scripts/Data_Vis/Section_4/RF_baseline_shap_cnv_rank_per.tsv", sep="\t", index_col=0)
snp_shap = pd.read_csv(
    "Scripts/Data_Vis/Section_4/RF_FS_shap_snp.tsv", sep="\t", index_col=0)  # these are average absolute shap values
pav_shap = pd.read_csv(
    "Scripts/Data_Vis/Section_4/RF_FS_shap_pav.tsv", sep="\t", index_col=0)
cnv_shap = pd.read_csv(
    "Scripts/Data_Vis/Section_4/RF_FS_shap_cnv.tsv", sep="\t", index_col=0)

# Calculate the ranks
snp_shap_env_rank_out = pd.DataFrame()
pav_shap_env_rank_out = pd.DataFrame()
cnv_shap_env_rank_out = pd.DataFrame()
for env in mapping.keys():
    snp_shap_env = snp_shap[[env]].dropna()
    pav_shap_env = pav_shap[[env]].dropna()
    cnv_shap_env = cnv_shap[[env]].dropna()
    #
    # drop features with 0 importance
    snp_shap_env = snp_shap_env[snp_shap_env[env] != 0]
    pav_shap_env = pav_shap_env[pav_shap_env[env] != 0]
    cnv_shap_env = cnv_shap_env[cnv_shap_env[env] != 0]
    #
    # calculate the ranks based on the absolute average SHAP values
    snp_shap_env_rank = snp_shap_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    pav_shap_env_rank = pav_shap_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    cnv_shap_env_rank = cnv_shap_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    #
    # concatenate all env ranks
    snp_shap_env_rank_out = pd.concat(
        [snp_shap_env_rank_out, snp_shap_env_rank], axis=1)
    pav_shap_env_rank_out = pd.concat(
        [pav_shap_env_rank_out, pav_shap_env_rank], axis=1)
    cnv_shap_env_rank_out = pd.concat(
        [cnv_shap_env_rank_out, cnv_shap_env_rank], axis=1)


# Map features to genes
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                       sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")

snp_gini_rank = map_snps[["snp", "gene"]].merge(
    snp_gini_env_rank_out, left_on="snp", right_index=True, how="right")
pav_gini_rank = map_orfs[["orf", "gene"]].merge(
    pav_gini_env_rank_out, left_on="orf", right_index=True, how="right")
cnv_gini_rank = map_orfs[["orf", "gene"]].merge(
    cnv_gini_env_rank_out, left_on="orf", right_index=True, how="right")

snp_shap_rank = map_snps[["snp", "gene"]].merge(
    snp_shap_env_rank_out, left_on="snp", right_index=True, how="right")
pav_shap_rank = map_orfs[["orf", "gene"]].merge(
    pav_shap_env_rank_out, left_on="orf", right_index=True, how="right")
cnv_shap_rank = map_orfs[["orf", "gene"]].merge(
    cnv_shap_env_rank_out, left_on="orf", right_index=True, how="right")

# Replace gene names with ORFs that aren't in S288C
pav_gini_rank["gene"] = pav_gini_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
cnv_gini_rank["gene"] = cnv_gini_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
pav_shap_rank["gene"] = pav_shap_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
cnv_shap_rank["gene"] = cnv_shap_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
pav_gini_rank = pav_gini_rank.drop(columns="orf").set_index("gene")
cnv_gini_rank = cnv_gini_rank.drop(columns="orf").set_index("gene")
pav_shap_rank = pav_shap_rank.drop(columns="orf").set_index("gene")
cnv_shap_rank = cnv_shap_rank.drop(columns="orf").set_index("gene")

# Get the SNP feature with the best rank for each gene
snp_gini_rank_max = snp_gini_rank.groupby("gene").min()
snp_shap_rank_max = snp_shap_rank.groupby("gene").min()
snp_gini_rank_max.drop(index="intergenic", inplace=True)
snp_shap_rank_max.drop(index="intergenic", inplace=True)
snp_gini_rank_max = snp_gini_rank_max.loc[~snp_gini_rank_max.index.str.contains(
    ","), :]  # drop snps that mapped to multiple genes
snp_shap_rank_max = snp_shap_rank_max.loc[~snp_shap_rank_max.index.str.contains(
    ","), :]
snp_gini_rank_max.drop(columns="snp", inplace=True)
snp_shap_rank_max.drop(columns="snp", inplace=True)

# Get the counts of how many environments each gene is in the optimized feature sets for
snp_gini_counts = snp_gini_rank_max.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 2
pav_gini_counts = pav_gini_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 3
cnv_gini_counts = cnv_gini_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 1
snp_shap_counts = snp_shap_rank_max.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 2
pav_shap_counts = pav_shap_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 3
cnv_shap_counts = cnv_shap_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 1

# Make Figure 3f
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(
    12, 7.5), sharex=True)  # env 1-35 x-axis
ks_df, ax = make_fig3f(
    snp_gini_rank_max.applymap(lambda x: 1 if x > 0 else 0), "snp", "gini", ax, 0, 0)
ks_df1, ax = make_fig3f(
    pav_gini_rank.applymap(lambda x: 1 if x > 0 else 0), "pav", "gini", ax, 0, 1)
ks_df2, ax = make_fig3f(
    cnv_gini_rank.applymap(lambda x: 1 if x > 0 else 0), "cnv", "gini", ax, 0, 2)

ks_df3, ax = make_fig3f(
    snp_shap_rank_max.applymap(lambda x: 1 if x > 0 else 0), "snp", "shap", ax, 1, 0)
ks_df4, ax = make_fig3f(
    pav_shap_rank.applymap(lambda x: 1 if x > 0 else 0), "pav", "shap", ax, 1, 1)
ks_df5, ax = make_fig3f(
    cnv_shap_rank.applymap(lambda x: 1 if x > 0 else 0), "cnv", "shap", ax, 1, 2)

plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_4/Figure_3f_gini_shap_FS_genes_shared_with_rand_CORRECTED.pdf")
plt.close('all')

# snp, gini: KS test between original and median randomized counts: KstestResult(statistic=0.5714285714285714, pvalue=3.868330132700767e-05, statistic_location=18, statistic_sign=1)
# pav, gini: KS test between original and median randomized counts: KstestResult(statistic=0.4444444444444444, pvalue=0.0005590074707047819, statistic_location=18, statistic_sign=1)
# cnv, gini: KS test between original and median randomized counts: KstestResult(statistic=0.5593355481727574, pvalue=6.704828589356284e-05, statistic_location=19, statistic_sign=1)
# snp, shap: KS test between original and median randomized counts: KstestResult(statistic=0.5714285714285714, pvalue=3.863356359698062e-05, statistic_location=18, statistic_sign=1)
# pav, shap: KS test between original and median randomized counts: KstestResult(statistic=0.4444444444444444, pvalue=0.0005590074707047819, statistic_location=18, statistic_sign=1)
# cnv, shap: KS test between original and median randomized counts: KstestResult(statistic=0.5212403100775194, pvalue=0.00014579537229090375, statistic_location=19, statistic_sign=1)

ks_df = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_10000permutations_snp_gini_CORRECTED.csv", index_col=0)
ks_df1 = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_10000permutations_pav_gini_CORRECTED.csv", index_col=0)
ks_df2 = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_10000permutations_cnv_gini_CORRECTED.csv", index_col=0)
ks_df3 = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_10000permutations_snp_shap_CORRECTED.csv", index_col=0)
ks_df4 = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_10000permutations_pav_shap_CORRECTED.csv", index_col=0)
ks_df5 = pd.read_csv(
    "Scripts/Data_Vis/Section_4/Figure_3f_ks2samp_rand_10000permutations_cnv_shap_CORRECTED.csv", index_col=0)

ks_df.pvalue.median()  # 1.8580137674662276e-35
ks_df1.pvalue.median()  # 1.524324440450238e-31
ks_df2.pvalue.median()  # 4.707639945685864e-71
ks_df3.pvalue.median()  # 1.1891384826327052e-35
ks_df4.pvalue.median()  # 1.524324440450238e-31
ks_df5.pvalue.median()  # 2.0975573145092256e-71
