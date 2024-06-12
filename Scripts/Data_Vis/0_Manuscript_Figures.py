#!/usr/bin/env python3
""" Yeast Fitness GxE Manuscript Figures """
import os
import re
import math
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy
from matplotlib import cm
from matplotlib.colors import Normalize, ListedColormap, rgb2hex
from matplotlib.patches import Rectangle
from scipy.interpolate import interpn
from scipy.stats import binned_statistic, spearmanr, linregress, ks_2samp
from tqdm import tqdm

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

### Functions for figure S1
def normalized_cbar(z):
    """ Normalize colorbar to density values """
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    return norm


def density_scatter(x , y, ax=None, fig=None, sort=True, bins=20):
    """
    Calculate density between x and y vectors and plot as a scatter plot
    source: https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density/53865762#53865762
    """
    if ax is None :
        fig , ax = plt.subplots()
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn((0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1])), data, np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)
    z[np.where(np.isnan(z))] = 0.0 # to be sure to plot all data
    if sort: # sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    # return z
    if ax is not None:
        ax.scatter(x, y, c=z, s=0.3)
        fig.colorbar(cm.ScalarMappable(norm=normalized_cbar(z), cmap="viridis"), ax=ax)
        return ax
    else:
        plt.scatter(x, y, c=z, s=0.3)
        fig.colorbar(cm.ScalarMappable(norm=normalized_cbar(z), cmap="viridis"))
        return fig


def map_clades(df_melt, clades):
    """ Add clade information to the isolates in the correlation matrices """
    try:
        df_melt = df_melt.merge(clades[["Standardized name", "Clades"]], left_on="level_0", right_on="Standardized name", how="left", copy=False)
        df_melt = df_melt.merge(clades[["Standardized name", "Clades"]], left_on="level_1", right_on="Standardized name", how="left")
        df_melt = df_melt[["level_0", "level_1", 0, "Clades_x", "Clades_y"]]
    except KeyError:
        df_melt = df_melt.merge(clades[["Standardized name", "Clades"]], left_on="ID", right_on="Standardized name", how="left", copy=False)
        df_melt = df_melt.merge(clades[["Standardized name", "Clades"]], left_on="level_1", right_on="Standardized name", how="left")
        df_melt = df_melt[["ID", "level_1", 0, "Clades_x", "Clades_y"]]
    df_melt.columns = ["level_0", "level_1", "PCC", "Clades_lvl0", "Clades_lvl1"]
    return df_melt


### Functions for figure 1
def plot_colortable(colors, *, ncols=4, sort_colors=True):
    """
    Plot a table of colors with their names
    source: https://matplotlib.org/stable/gallery/color/named_colors.html 
    """
    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    n = len(colors)
    nrows = math.ceil(n / ncols)
    width = cell_width * ncols + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72
    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-margin)/height)
    ax.set_xlim(0, cell_width * ncols)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    for i, name in enumerate(colors):
        row = i % nrows
        col = i // nrows
        y = row * cell_height
        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7
        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')
        ax.add_patch(
            Rectangle(xy=(swatch_start_x, y-9), width=swatch_width,
                      height=18, facecolor=colors[i], edgecolor='0.7')
        )
    return fig


### Functions for figure 2a
def make_fig2a(df, lvl_vals, fig1_ax, fig2_ax, fig3_ax, fig4_ax, ax_row_idx, legend_label):
    """ Make Figure 2a """
    ###### Randomize environment columns individually
    n = 1000 # number of iterations to randomize counts
    counts_rand = {lvl_val:{str(i):None for i in range(n)} for lvl_val in lvl_vals} # empty dictionary to store new env counts for each gene after randomization
    ks_res = {lvl_val:{str(i):None for i in range(n)} for lvl_val in lvl_vals} # KS test results
    for lvl_val in lvl_vals:
        lvl_df = df.loc[:, df.columns.get_level_values("Data")==lvl_val] # subset df for each data type
        lvl_df = lvl_df.drop_duplicates(keep=False) # drop genes with 0 in all columns (these genes are top genes for other data types)
        orig_env_counts = lvl_df.sum(axis=1) # original number of environments per gene
        for i in tqdm(range(n)):
            df_rand = {}
            for j,col in enumerate(lvl_df.columns):
                df_rand[col] = np.random.default_rng(seed=j+i).permutation(lvl_df[col]) # randomized sample of if gene is a top gene or not (1/0)
            df_rand = pd.DataFrame(df_rand)
            df_rand.columns.names = ["Env", "Data", "Info"]
            # Realculate number of environments a gene is a top feature for
            env_counts = df_rand.sum(axis=1) # randomized number of environments per gene
            counts_rand[lvl_val][str(i)] = env_counts # Save results
            # Calculate KS-statistic
            ks_res[lvl_val][str(i)] = ks_2samp(env_counts, orig_env_counts) # 2-sample KS test
    # Plot KS test results
    ks_res_df = pd.DataFrame() # results of KS test for all data type env count frequency distributions
    for lvl_val in ks_res.keys():
        ks_df = pd.DataFrame(ks_res[lvl_val], index=["D", "p-val"]).T
        ks_df["Data"] = lvl_val
        ks_res_df = pd.concat([ks_res_df, ks_df], axis=0)
        ks_df.hist(figsize=(8.5,4))
        plt.savefig(f"Scripts/Data_Vis/Figure_2a_ks2samp_{n}permutations_{lvl_val}_{legend_label}.pdf")
        plt.close()
    ##### Randomize groups of the same environments (but different concentrations) together and all others individually
    env_groups = pd.DataFrame({"YPD14":1, "YPD40":1, "YPD42":1, "YPDANISO10":2,
                               "YPDANISO20":2, "YPDANISO50":2, "YPDBENOMYL200":3,
                               "YPDBENOMYL500":3, "YPDCAFEIN40":4, "YPDCAFEIN50":4,
                               "YPDCHX05":5, "YPDCHX1":5, "YPDFORMAMIDE4":6,
                               "YPDFORMAMIDE5":6, "YPDNACL15M":7, "YPDNACL1M":7,
                               "YPDETOH":8, "YPETHANOL":8}, index=["Group"]).T
    df_gps = df.copy(deep=True)
    counts_rand_gps = {lvl_val:{str(i):None for i in range(n)} for lvl_val in lvl_vals}
    ks_res_gps = {lvl_val:{str(i):None for i in range(n)} for lvl_val in lvl_vals}
    ks_rand_vs_gps = {lvl_val:{str(i):None for i in range(n)} for lvl_val in lvl_vals}
    for lvl_val in lvl_vals:
        lvl_df = df_gps.loc[:, df_gps.columns.get_level_values("Data")==lvl_val] # subset df for each data type
        lvl_df = lvl_df.drop_duplicates(keep=False) # drop genes with 0 in all columns (these genes are top genes for other data types)
        orig_env_counts = lvl_df.sum(axis=1) # original number of environments per gene
        for i in tqdm(range(n)):
            df_rand_gps = pd.DataFrame()
            for group in env_groups.Group.unique(): # randomize rows in groups of columns together
                group_df = lvl_df.loc[:,lvl_df.columns.get_level_values("Env").isin(env_groups.loc[env_groups.Group==group].index)]
                df_rand_gps = pd.concat([df_rand_gps, pd.DataFrame(np.random.default_rng(seed=group+i).permutation(group_df))], axis=1, ignore_index=True)
            for j,col in enumerate(lvl_df.columns): # randomize the remaining columns individually
                if col[0] not in env_groups.index:
                    df_rand_gps = pd.concat([df_rand_gps, pd.DataFrame(np.random.default_rng(seed=j+i+9).permutation(lvl_df[col]))], axis=1, ignore_index=True)
            # Recalculate number of environments a gene is a top feature for & the KS statistic
            env_counts = df_rand_gps.sum(axis=1)
            counts_rand_gps[lvl_val][str(i)] = env_counts
            ks_res_gps[lvl_val][str(i)] = ks_2samp(env_counts, orig_env_counts)
            ks_rand_vs_gps[lvl_val][str(i)] = ks_2samp(counts_rand[lvl_val][str(i)], env_counts)
    ks_res_gps_df = pd.DataFrame() # results of KS test for all data type env count frequency distributions
    ks_rand_vs_gps_df = pd.DataFrame() # results of KS test between non-grouped permutations and grouped permutations
    for lvl_val in ks_res_gps.keys():
        ks_df = pd.DataFrame(ks_res_gps[lvl_val], index=["D", "p-val"]).T
        ks_df["Data"] = lvl_val
        ks_res_gps_df = pd.concat([ks_res_gps_df, ks_df], axis=0)
        ks_df.hist(figsize=(8.5,4))
        plt.savefig(f"Scripts/Data_Vis/Figure_2a_ks2samp_{n}permutations_{lvl_val}_{legend_label}_grouped.pdf")
        plt.close()
        ks_df2 = pd.DataFrame(ks_rand_vs_gps[lvl_val], index=["D", "p-val"]).T
        ks_df2["Data"] = lvl_val
        ks_rand_vs_gps_df = pd.concat([ks_rand_vs_gps_df, ks_df2], axis=0)
        ks_df2.hist(figsize=(8.5,4))
        plt.savefig(f"Scripts/Data_Vis/Figure_2a_ks2samp_{n}permutations_{lvl_val}_{legend_label}_nongrouped_vs_grouped.pdf")
        plt.close()
    ##### Plot histograms of counts
    for i,lvl_val in enumerate(lvl_vals):
        # Calculate the median number of environments across all randomizations and the original counts
        counts_rand_freq = pd.DataFrame(counts_rand[lvl_val]).apply(pd.value_counts) # get frequency of env counts at each randomization iteration
        counts_rand_freq_median = counts_rand_freq.median(axis=1).astype("int") # calculate median counts for each bin
        counts_rand_freq_median.drop(index=0) # drop the 0 env count
        log_idx_counts_rand_median = np.log10(counts_rand_freq_median.index.values) # logged x-axis
        counts_rand_gps_freq = pd.DataFrame(counts_rand_gps[lvl_val]).apply(pd.value_counts)
        counts_rand_gps_freq_median = counts_rand_gps_freq.median(axis=1).astype("int").drop(index=0)
        log_idx_counts_rand_gps = np.log10(counts_rand_gps_freq_median.index.values)
        lvl_df = df.loc[:, df.columns.get_level_values("Data")==lvl_val] # subset df for each data type
        lvl_df = lvl_df.drop_duplicates(keep=False) # drop genes with 0 in all columns (these genes are top genes for other data types)
        counts = lvl_df.sum(axis=1) # original number of environments per gene
        counted_envs = lvl_df.apply(lambda x: " // ".join(x.loc[:,x != 0].index.get_level_values(0)), axis=1) # counted environments per gene
        ## fig1 logged x-axis
        np.log10(counts).plot.hist(bins=35, alpha=0.5, color="blue", label=legend_label, ax=fig1_ax[ax_row_idx][i]) # original counts
        pd.concat([np.log10(counts), counted_envs], ignore_index=False, axis=1).to_csv(f"Scripts/Data_Vis/Figure2b_logged_x-axis_{lvl_val}_{ax_row_idx}_{i}_data.csv")
        fig1_ax[ax_row_idx][i].axvline(np.log10(counts.median()), ls="-", color="blue", linewidth=1)
        fig1_ax[ax_row_idx][i].plot(log_idx_counts_rand_median, counts_rand_freq_median, ls="-",
                                   alpha=0.5, color="black", label="median randomized counts") # randomized counts
        fig1_ax[ax_row_idx][i].plot(log_idx_counts_rand_gps, counts_rand_gps_freq_median, ls=":",
                                   alpha=0.5, color="black", label="median grouped randomized counts") # grouped randomized counts
        fig1_ax[ax_row_idx][i].set_title(f"{lvl_val}: logged x-axis")
        fig1_ax[ax_row_idx][i].legend()
        ## fig2 regular x-axis (environment counts 1 to 35)
        counts.astype("int").plot.hist(bins=35, range=[1,35], alpha=0.5, color="blue",
                                       label=legend_label, ax=fig2_ax[ax_row_idx][i])
        pd.concat([counts, counted_envs], ignore_index=False, axis=1).to_csv(f"Scripts/Data_Vis/Figure2b_{lvl_val}_{ax_row_idx}_{i}_data.csv")
        fig2_ax[ax_row_idx][i].axvline(counts.median(), ls="-", color="blue", linewidth=1)
        counts_rand_freq_median.plot.line(alpha=0.5, color="black",
                                          label="median randomized counts",
                                          ax=fig2_ax[ax_row_idx][i]) # randomized counts
        counts_rand_gps_freq_median.plot.line(ls=":", alpha=0.5, color="black",
                                              label="median grouped randomized counts",
                                              ax=fig2_ax[ax_row_idx][i]) # grouped randomized counts
        fig2_ax[ax_row_idx][i].set_title(f"{lvl_val}: regular x-axis")
        fig2_ax[ax_row_idx][i].legend()
        ## fig3 Restrict range to range=[1,20] for logged x-axis
        np.log10(counts).plot.hist(bins=35, alpha=0.5, range=[np.log10(1),np.log10(20)],
                                      color="blue", label=legend_label, ax=fig3_ax[ax_row_idx][i])
        fig3_ax[ax_row_idx][i].axvline(np.log10(counts.median()), ls="-", color="blue", linewidth=1)
        fig3_ax[ax_row_idx][i].plot(log_idx_counts_rand_median, counts_rand_freq_median, ls="-",
                                   alpha=0.5, color="black", label="median randomized counts") # randomized counts
        fig3_ax[ax_row_idx][i].plot(log_idx_counts_rand_gps, counts_rand_gps_freq_median, ls=":",
                                   alpha=0.5, color="black", label="median grouped randomized counts") # grouped randomized counts
        fig3_ax[ax_row_idx][i].set_title(f"{lvl_val}: logged x-axis; environment counts (1-20)")
        fig3_ax[ax_row_idx][i].legend()
        ## fig4 Restrict range to range=[1,20], regular x-axis
        counts.astype("int").plot.hist(bins=35, range=[1,20], alpha=0.5, color="blue",
                                          label=legend_label, ax=fig4_ax[ax_row_idx][i])
        fig4_ax[ax_row_idx][i].axvline(counts.median(), ls="-", color="blue", linewidth=1)
        counts_rand_freq_median.plot.line(alpha=0.5, color="black",
                                          label="median randomized counts", ax=fig4_ax[ax_row_idx][i]) # randomized counts
        counts_rand_freq_median.plot.line(ls=":", alpha=0.5, color="black",
                                          label="median grouped randomized counts",ax=fig4_ax[ax_row_idx][i]) # grouped randomized counts
        fig4_ax[ax_row_idx][i].set_title(f"{lvl_val}: regular x-axis; environment counts (1-20)")
        fig4_ax[ax_row_idx][i].legend()
    return ks_res_df, ks_res_gps_df, ks_rand_vs_gps_df, fig1_ax, fig2_ax, fig3_ax, fig4_ax  


if __name__ == "__main__":
    # Environment truncated and full names mapping
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

    ################################################################################
    # Miscellaneous figures
    ################################################################################
    #%%
    import pandas as pd
    import matplotlib.pyplot as plt
    # heterozygosity (F-statistic) for each isolate
    het = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/haplotypes/geno_11292023.het", sep="\t")
    het["F"].hist(bins=20)
    plt.savefig("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/heterozygosity_hist.pdf")
    plt.close()
    #%%
    # Hardy-Weinberg Equilibrium measures for each marker
    hwe = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/haplotypes/geno_11292023.hwe", sep="\t")
    het_counts = hwe["OBS(HOM1/HET/HOM2)"].apply(lambda x: int(x.split("/")[1])) # heterozygous isolates count
    print(sum(het_counts==0))
    fig = plt.figure(figsize=(8, 6))
    plt.hist(het_counts, bins=20)
    plt.xticks(rotation=90)
    plt.grid(visible=True, which="both")
    plt.savefig("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/heterozygosity_counts_hist.pdf")
    #%%
    ################################################################################
    # Figure 1
    ################################################################################
    pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # fitness data
    pheno.rename(columns=mapping, inplace=True) # full environment names
    pCorEnvs = pheno.corr(method="pearson") # fitness correlations between environments
    pCorEnvs.to_csv("Data/Peter_2018/fitness_correlations_env_pair.csv")

    ## S?. Distributions of fitness values for training and testing sets
    test_ids = pd.read_csv("Data/Peter_2018/Test.txt", header=None) # test set ids
    train_y = pheno.loc[~pheno.index.isin(test_ids[0]),:] # training set
    test_y = pheno.loc[pheno.index.isin(test_ids[0]),:] # test set
    
    fig, axis = plt.subplots(7, 5, figsize=(15,21))
    train_y.hist(bins=25, ax=axis)
    plt.savefig("Scripts/Data_Vis/Section_3/Fitness_training_set_distributions.pdf")
    plt.close()
    
    fig, axis = plt.subplots(7, 5, figsize=(15,21))
    test_y.hist(bins=25, ax=axis)
    plt.savefig("Scripts/Data_Vis/Section_3/Fitness_test_set_distributions.pdf")
    plt.close()
    
    ## S1a. Violin plot of fitness distributions in each environment ### will move this to somewhere else
    sns.violinplot(data=pheno)
    violin_ax.set_xticks([]) # hide x-axis ticks
    violin_ax.set_ylabel("Fitness")
    plt.savefig("Scripts/Data_Vis/Figure_S1ab.pdf", width=7, height=8, bbox_inches="tight")
    plt.close()

    ## S1b. Heatmap of fitness correlations between environments
    # Even though the method/metric is the same, the figure I fixed in Adobe was made in R and it's different from the Python version.
    sns.clustermap(pCorEnvs, cmap="RdBu_r", vmin=-1, vmax=1, cbar_kws={"shrink":0.5}, xticklabels=True, method="complete", metric="euclidean")
    plt.savefig("Scripts/Data_Vis/Figure_S1b_v2.pdf") # keep the R code figure instead. See Manuscript_Figures.R
    plt.close()

    ## S1c. Heatmaps of isolate pair data correlations
    kinship = pd.read_csv("Data/Peter_2018/kinship.csv", index_col=0)
    snps = dd.read_csv("Data/Peter_2018/geno.csv", sample=10000000, blocksize=None)
    snps = snps.compute()
    snps = pd.read_csv("Data/Peter_2018/geno.csv", chunksize=100, iterator=True)
    snps = pd.concat(snps, ignore_index=True)
    pavs = pd.read_csv("Data/Peter_2018/ORFs_pres_abs.csv", index_col=0)
    pavs = pavs.loc[kinship.index,:] # sort ORFs by kinship isolate order
    cnvs = pd.read_csv("Data/Peter_2018/ORFs_no_NA.csv", index_col=0)
    cnvs = cnvs.loc[kinship.index,:] # sort CNVs by kinship isolate order
    clades = pd.read_excel("Data/Peter_2018/Peter_2018_Supplementary_Tables.xls", sheet_name="Table S1", skiprows=3, nrows=1011) # isolate clades

    # summary statistics for each isolate and each clade
    # fitness, snp heterozygosity, pav stats, cnv stats, number of isolates
    clades = clades[["Standardized name", "Clades"]] # subset relevant columns
    pheno.merge(clades, left_index=True, right_on="Standardized name").groupby("Clades").describe().to_csv("Scripts/Data_Vis/pheno_clades_stats.csv")
    snps.merge(clades, left_index=True, right_on="Standardized name").groupby("Clades").value_counts()
    pavs.merge(clades, left_index=True, right_on="Standardized name").groupby("Clades").describe()
    cnvs.merge(clades, left_index=True, right_on="Standardized name").groupby("Clades").value_counts()

    # cluster data
    kin_dendrogram = hierarchy.dendrogram(hierarchy.linkage(kinship, method="average"))
    kin_row_order = kin_dendrogram["leaves"]
    pheno_corr = pheno.T.corr(method="pearson")
    pavs_corr = pavs.T.corr(method="pearson")
    cnvs_corr = cnvs.T.corr(method="pearson")
    pheno_corr.to_csv("Data/Peter_2018/fitness_correlations_isolate_pair.csv")
    pavs_corr.to_csv("Data/Peter_2018/ORFs_pres_abs_correlations_isolate_pair.csv")
    cnvs_corr.to_csv("Data/Peter_2018/ORFs_copy_num_correlations_isolate_pair.csv")
    sum(kinship.index==pheno_corr.index)
    sum(kinship.index==pavs_corr.index)
    sum(kinship.index==cnvs_corr.index)
    sum(kinship.columns==pheno_corr.columns)
    sum(kinship.columns==pavs_corr.columns)
    sum(kinship.columns==cnvs_corr.columns)

    # Map clades to isolates and create a colormap
    clades = clades[["Standardized name", "Clades"]] # subset relevant columns
    clades = clades.loc[clades["Standardized name"].isin(kinship.index)] # diploid isolates
    clades.set_index("Standardized name", inplace=True)
    clades = clades.loc[kinship.index,:] # reorder by original order
    clades = clades.iloc[kin_row_order,:] # reorder by dendrogram order
    clades.loc[clades.Clades.isna(),"Clades"] = "Unknown" # replace NaN with Unknown
    
    cmap = cm.get_cmap('hsv', clades.Clades.nunique()) # sample hsv color palette
    color_list = {clades.Clades.unique()[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
    row_colors = clades['Clades'].map(color_list) # assign colors to isolates
    
    plot_colortable(row_colors.unique()) # plot color table
    plt.savefig("Scripts/Data_Vis/clade_colors.png")
    plt.close()
    with open("Scripts/Data_Vis/clade_colors.json", "w") as f:
        json.dump(color_list, f, indent=4) # save color dictionary
    
    # plot
    fig, ax = plt.subplots(2, 2)
    a=sns.heatmap(kinship.iloc[kin_row_order, kin_row_order], cmap="RdBu_r", center=0, square=True, ax=ax[0][0], xticklabels=False, yticklabels=False)
    ax[0][0].tick_params(axis="y", which="major", pad=20, length=0) # extra padding for row colors
    for i, color in enumerate(row_colors):
        ax[0][0].add_patch(plt.Rectangle(xy=(-0.05, i), width=0.05, height=1, color=color, lw=0,
                                transform=ax[0][0].get_yaxis_transform(), clip_on=False))
    b=sns.heatmap(pheno_corr.iloc[kin_row_order, kin_row_order], cmap="RdBu_r", center=pheno_corr.mean().mean(), square=True, ax=ax[0][1], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=pheno_melt_freqs.index.min().left, vmax=pheno_melt_freqs.index.max().right)
    c=sns.heatmap(pavs_corr.iloc[kin_row_order, kin_row_order], cmap="RdBu_r", center=0.8, square=True, ax=ax[1][0], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=pavs_melt_freqs.index.min().left, vmax=pavs_melt_freqs.index.max().right),
    d=sns.heatmap(cnvs_corr.iloc[kin_row_order, kin_row_order], cmap="RdBu_r", center=0.8, square=True, ax=ax[1][1], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=cnvs_melt_freqs.index.min().left, vmax=cnvs_melt_freqs.index.max().right),
    fig.tight_layout()
    # fig.savefig("Scripts/Data_Vis/Figure_S1c_data_correlations2.png", dpi=300)
    fig.savefig("Scripts/Data_Vis/Figure_S1c_data_correlations3.png", dpi=300) # with row clade colors
    plt.close()
    
    ################################################################################
    # Figure S1
    ################################################################################
    ## S1d. data pair correlations
    clades = pd.read_excel("Data/Peter_2018/Peter_2018_Supplementary_Tables.xls", sheet_name="Table S1", skiprows=3, nrows=1011)
    kinship = pd.read_csv("Data/Peter_2018/kinship.csv", index_col=0)
    pheno_corr = pd.read_csv("Data/Peter_2018/fitness_correlations_isolate_pair.csv", index_col=0)
    pavs_corr = pd.read_csv("Data/Peter_2018/ORFs_pres_abs_correlations_isolate_pair.csv", index_col=0)
    cnvs_corr = pd.read_csv("Data/Peter_2018/ORFs_copy_num_correlations_isolate_pair.csv", index_col=0)
    
    # Remove duplicate pairs ex. AAC AAC; AVD APE and APE AVD
    kin_melt = kinship.where(np.triu(kinship, k=1).astype(bool)).stack()
    pheno_melt = pheno_corr.where(np.triu(pheno_corr, k=1).astype(bool)).stack()
    pavs_melt = pavs_corr.where(np.triu(pavs_corr, k=1).astype(bool)).stack()
    cnvs_melt = cnvs_corr.where(np.triu(cnvs_corr, k=1).astype(bool)).stack()

    # Sanity check (order matters for spearman's rank correlation)
    sum(kin_melt.index==pheno_melt.index)
    sum(kin_melt.index==pavs_melt.index)
    sum(kin_melt.index==cnvs_melt.index)
    kin_melt = kin_melt.reset_index()
    pheno_melt = pheno_melt.reset_index()
    pavs_melt = pavs_melt.reset_index()
    cnvs_melt = cnvs_melt.reset_index()
    kin_melt["rank"] = kin_melt.rank(numeric_only=True, ascending=False)
    pheno_melt["rank"] = pheno_melt.rank(numeric_only=True, ascending=False)
    pavs_melt["rank"] = pavs_melt.rank(numeric_only=True, ascending=False)
    cnvs_melt["rank"] = cnvs_melt.rank(numeric_only=True, ascending=False)
    kin_melt.to_csv("Scripts/Data_Vis/kin_melt.csv")
    pheno_melt.to_csv("Scripts/Data_Vis/pheno_melt.csv")
    pavs_melt.to_csv("Scripts/Data_Vis/pavs_melt.csv")
    cnvs_melt.to_csv("Scripts/Data_Vis/cnvs_melt.csv")

    data_pair_corr = pd.DataFrame(columns=["comparison", "rho", "p-value", "note"])
    # Plot
    fig, ax = plt.subplots(2, 3, figsize=(12, 6))
    ax[0][0] = density_scatter(kin_melt[0], pheno_melt[0], ax[0][0], fig, bins=50)
    rho, p = spearmanr(kin_melt['rank'], pheno_melt['rank']) # calculate spearman's rho
    stats.ttest_ind(kin_melt[0], pheno_melt[0], equal_var=False, alternative='two-sided') # p-values are all 0
    print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p)) # rho = .269, p-val = 0.000E+00 ; in R, cor.test p-values are all < 2.2e-16
    
    ax[0][1] = density_scatter(pavs_melt[0], pheno_melt[0], ax[0][1], fig, bins=50)
    rho, p = spearmanr(pavs_melt['rank'], pheno_melt['rank'])
    stats.ttest_ind(pavs_melt[0], pheno_melt[0], equal_var=False, alternative='two-sided')
    print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p)) # rho = .299, p-val = 0.000E+00
    
    ax[0][2] = density_scatter(cnvs_melt[0], pheno_melt[0], ax[0][2], fig, bins=50)
    rho, p = spearmanr(cnvs_melt['rank'], pheno_melt['rank'])
    stats.ttest_ind(cnvs_melt[0], pheno_melt[0], equal_var=False, alternative='two-sided')
    print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p)) # rho = 0.138, p-val = 0.000E+00

    ax[1][0] = density_scatter(pavs_melt[0], kin_melt[0], ax[1][0], fig, bins=50)
    rho, p = spearmanr(pavs_melt['rank'], kin_melt['rank'])
    stats.ttest_ind(pavs_melt[0], kin_melt[0], equal_var=False, alternative='two-sided')
    print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p)) # rho = 0.475, p-val = 0.000E+00

    ax[1][1] = density_scatter(cnvs_melt[0], kin_melt[0], ax[1][1], fig, bins=50)
    rho, p = spearmanr(cnvs_melt['rank'], kin_melt['rank'])
    stats.ttest_ind(cnvs_melt[0], kin_melt[0], equal_var=False, alternative='two-sided')
    print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p)) # rho = 0.096, p-val = 0.000E+00

    ax[1][2] = density_scatter(cnvs_melt[0], pavs_melt[0], ax[1][2], fig, bins=50)
    rho, p = spearmanr(cnvs_melt['rank'], pavs_melt['rank'])
    stats.ttest_ind(cnvs_melt[0], pavs_melt[0], equal_var=False, alternative='two-sided')
    print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p)) # rho = 0.172, p-val = 0.000E+00
    fig.tight_layout()
    #fig.savefig("Scripts/Data_Vis/Figure_S1d_data_pair_correlations.png", dpi=300) # has issue with isolate pairs not matching between data types
    #fig.savefig("Scripts/Data_Vis/Figure_S1d_data_pair_correlations3.png", dpi=300) # after fixing issue with rows not matching for the rho calculation
    fig.savefig("Scripts/Data_Vis/Figure_S1d_data_pair_correlations3.png", dpi=300) # removed duplicate pairs, final version
    plt.close()

    # Map isolates to clades
    kin_melt = map_clades(kin_melt, clades)
    pheno_melt = map_clades(pheno_melt, clades)
    pavs_melt = map_clades(pavs_melt, clades)
    cnvs_melt = map_clades(cnvs_melt, clades)
    
    # Create individual plots for each clade (pairs in different clades will be excluded for now)
    for clade in kin_melt.Clades_lvl0.unique():
        kin_sub = kin_melt.loc[(kin_melt.Clades_lvl0==clade) & (kin_melt.Clades_lvl1==clade),]
        pheno_sub = pheno_melt.loc[(pheno_melt.Clades_lvl0==clade) & (pheno_melt.Clades_lvl1==clade),]
        pavs_sub = pavs_melt.loc[(pavs_melt.Clades_lvl0==clade) & (pavs_melt.Clades_lvl1==clade),]
        cnvs_sub = cnvs_melt.loc[(cnvs_melt.Clades_lvl0==clade) & (cnvs_melt.Clades_lvl1==clade),]
        sum(kin_sub.level_0==pheno_sub.level_0) ; sum(kin_sub.level_1==pheno_sub.level_1)
        sum(kin_sub.level_0==pavs_sub.level_0) ; sum(kin_sub.level_1==pavs_sub.level_1)
        sum(kin_sub.level_0==cnvs_sub.level_0) ; sum(kin_sub.level_1==cnvs_sub.level_1)
        print('shape', kin_sub.shape, pheno_sub.shape, pavs_sub.shape, cnvs_sub.shape)
        if (kin_sub.shape[0] > 0) and (clade is not None):
            clade_lbl = clade.replace("/", "_")
            fig, ax = plt.subplots(2, 3, figsize=(12, 6))
            ax[0][0].scatter(x=kin_sub["PCC"], y=pheno_sub["PCC"], alpha=0.3) ; ax[0][0].set_xlim(-1.5,4.2); ax[0][0].set_ylim(-.5,1)
            ax[0][1].scatter(x=pavs_sub["PCC"], y=pheno_sub["PCC"], alpha=0.3) ; ax[0][1].set_xlim(.6,1) ; ax[0][1].set_ylim(-.5,1)
            ax[0][2].scatter(x=cnvs_sub["PCC"], y=pheno_sub["PCC"], alpha=0.3) ; ax[0][2].set_xlim(.2,1) ; ax[0][2].set_ylim(-.5,1)
            ax[1][0].scatter(x=pavs_sub["PCC"], y=kin_sub["PCC"], alpha=0.3) ; ax[1][0].set_xlim(.6,1) ; ax[1][0].set_ylim(-1.5,4.2)
            ax[1][1].scatter(x=cnvs_sub["PCC"], y=kin_sub["PCC"], alpha=0.3) ; ax[1][1].set_xlim(.2,1) ; ax[1][1].set_ylim(-1.5,4.2)
            ax[1][2].scatter(x=cnvs_sub["PCC"], y=pavs_sub["PCC"], alpha=0.3) ; ax[1][2].set_xlim(.2,1) ; ax[1][2].set_ylim(.6,1)
            fig.tight_layout()
            fig.savefig(f"Scripts/Data_Vis/Figure_S1d_data_pair_correlations_{clade_lbl}.png", dpi=300) # removed duplicate pairs, final version
            plt.close()

    ## How do spearman's rho values change when kinship row order is not considered? They don't change
    clades = pd.read_excel("Data/Peter_2018/Peter_2018_Supplementary_Tables.xls", sheet_name="Table S1", skiprows=3, nrows=1011)
    kinship = pd.read_csv("Data/Peter_2018/kinship.csv", index_col=0)
    pheno_corr = pd.read_csv("Data/Peter_2018/fitness_correlations_isolate_pair.csv", index_col=0)
    pavs_corr = pd.read_csv("Data/Peter_2018/ORFs_pres_abs_correlations_isolate_pair.csv", index_col=0)
    cnvs_corr = pd.read_csv("Data/Peter_2018/ORFs_copy_num_correlations_isolate_pair.csv", index_col=0)
    # consider pav clustered order
    pav_dendrogram = hierarchy.dendrogram(hierarchy.linkage(pavs_corr, method="average")) # cluster
    pav_row_order = pav_dendrogram["leaves"]
    sum(kinship.index==pavs_corr.index) # sanity check
    sum(kinship.columns==pavs_corr.columns)
    kinship = kinship.iloc[pav_row_order, pav_row_order] # reorder by pav cluster order
    kin_melt = kinship.where(np.triu(kinship, k=1).astype(bool)).stack()
    sum(pheno_corr.index==pavs_corr.index)
    sum(pheno_corr.columns==pavs_corr.columns)
    pheno_corr = pheno_corr.iloc[pav_row_order, pav_row_order]
    pheno_melt = pheno_corr.where(np.triu(pheno_corr, k=1).astype(bool)).stack()
    sum(cnvs_corr.index==pavs_corr.index)
    sum(cnvs_corr.columns==pavs_corr.columns)
    cnvs_corr = cnvs_corr.iloc[pav_row_order, pav_row_order]
    cnvs_melt = cnvs_corr.where(np.triu(cnvs_corr, k=1).astype(bool)).stack()
    pavs_corr = pavs_corr.iloc[pav_row_order, pav_row_order]
    pavs_melt = pavs_corr.where(np.triu(pavs_corr, k=1).astype(bool)).stack()
    kin_melt["rank"] = kin_melt.rank(numeric_only=True, ascending=False) # rank values
    pheno_melt["rank"] = pheno_melt.rank(numeric_only=True, ascending=False)
    pavs_melt["rank"] = pavs_melt.rank(numeric_only=True, ascending=False)
    cnvs_melt["rank"] = cnvs_melt.rank(numeric_only=True, ascending=False)
    spearmanr(kin_melt, pheno_melt) # rho = 0.269, p-val = 0.000E+00
    spearmanr(kin_melt, pavs_melt) # rho = 0.475, p-val = 0.000E+00
    spearmanr(kin_melt, cnvs_melt) # rho = 0.096, p-val = 0.000E+00
    spearmanr(pheno_melt, pavs_melt) # rho = 0.299, p-val = 0.000E+00
    spearmanr(pheno_melt, cnvs_melt) # rho = 0.138, p-val = 0.000E+00
    spearmanr(pavs_melt, cnvs_melt) # rho = 0.172, p-val = 0.000E+00
    
    # Map clades to isolates and create a colormap
    clades = clades[["Standardized name", "Clades"]] # subset relevant columns
    clades = clades.loc[clades["Standardized name"].isin(kinship.index)] # diploid isolates
    clades.set_index("Standardized name", inplace=True)
    kin_dendrogram = hierarchy.dendrogram(hierarchy.linkage(kinship, method="average")) # cluster
    kin_row_order = kin_dendrogram["leaves"]
    clades = clades.loc[kinship.index,:] # reorder by original order
    clades = clades.loc[kinship.iloc[kin_row_order, kin_row_order].index,:] # reorder by dendrogram order
    clades.loc[clades.Clades.isna(),"Clades"] = "Unknown" # replace NaN with Unknown
    cmap = cm.get_cmap('hsv', clades.Clades.nunique()) # sample hsv color palette
    color_list = {clades.Clades.unique()[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
    row_colors = clades['Clades'].map(color_list) # assign colors to isolates
    row_colors = row_colors[pavs_corr.index]
    
    # instead, redraw correlation heatmaps (Figure 1)
    fig, ax = plt.subplots(2, 2)
    a=sns.heatmap(kinship, cmap="RdBu_r", center=0, square=True, ax=ax[0][0], xticklabels=False, yticklabels=False)
    b=sns.heatmap(pheno_corr, cmap="RdBu_r", center=pheno_corr.mean().mean(), square=True, ax=ax[0][1], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=pheno_melt_freqs.index.min().left, vmax=pheno_melt_freqs.index.max().right)
    c=sns.heatmap(pavs_corr, cmap="RdBu_r", center=0.8, square=True, ax=ax[1][0], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=pavs_melt_freqs.index.min().left, vmax=pavs_melt_freqs.index.max().right),
    ax[1][0].tick_params(axis="y", which="major", pad=20, length=0) # extra padding for row colors
    for i, color in enumerate(row_colors):
        ax[1][0].add_patch(plt.Rectangle(xy=(-0.05, i), width=0.05, height=1, color=color, lw=0,
                                transform=ax[1][0].get_yaxis_transform(), clip_on=False))
    d=sns.heatmap(cnvs_corr, cmap="RdBu_r", center=0.8, square=True, ax=ax[1][1], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=cnvs_melt_freqs.index.min().left, vmax=cnvs_melt_freqs.index.max().right),
    fig.tight_layout()
    fig.savefig("Scripts/Data_Vis/Figure_S1c_data_correlations_pav_order.png", dpi=300) # with row clade colors
    plt.close()

    # consider cnv cluster order
    cnv_dendrogram = hierarchy.dendrogram(hierarchy.linkage(cnvs_corr, method="average")) # cluster
    cnv_row_order = cnv_dendrogram["leaves"]
    kinship = kinship.iloc[cnv_row_order, cnv_row_order] # reorder by pav cluster order
    pheno_corr = pheno_corr.iloc[cnv_row_order, cnv_row_order]
    cnvs_corr = cnvs_corr.iloc[cnv_row_order, cnv_row_order]
    pavs_corr = pavs_corr.iloc[cnv_row_order, cnv_row_order]

    # redraw correlation heatmaps (Figure 1)
    row_colors = row_colors[cnvs_corr.index]
    fig, ax = plt.subplots(2, 2)
    a=sns.heatmap(kinship, cmap="RdBu_r", center=0, square=True, ax=ax[0][0], xticklabels=False, yticklabels=False)
    b=sns.heatmap(pheno_corr, cmap="RdBu_r", center=pheno_corr.mean().mean(), square=True, ax=ax[0][1], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=pheno_melt_freqs.index.min().left, vmax=pheno_melt_freqs.index.max().right)
    c=sns.heatmap(pavs_corr, cmap="RdBu_r", center=0.8, square=True, ax=ax[1][0], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=pavs_melt_freqs.index.min().left, vmax=pavs_melt_freqs.index.max().right),
    d=sns.heatmap(cnvs_corr, cmap="RdBu_r", center=0.8, square=True, ax=ax[1][1], xticklabels=False, yticklabels=False) # norm=plt.Normalize(vmin=cnvs_melt_freqs.index.min().left, vmax=cnvs_melt_freqs.index.max().right),
    ax[1][1].tick_params(axis="y", which="major", pad=20, length=0) # extra padding for row colors
    for i, color in enumerate(row_colors): # see lines 311-320
        ax[1][1].add_patch(plt.Rectangle(xy=(-0.05, i), width=0.05, height=1, color=color, lw=0,
                                transform=ax[1][1].get_yaxis_transform(), clip_on=False))
    fig.tight_layout()
    fig.savefig("Scripts/Data_Vis/Figure_S1c_data_correlations_cnv_order.png", dpi=300) # with row clade colors
    plt.close()

    ############################################################################
    # Figure S4 
    ############################################################################
    # Read feature importance rankings
    fs_combined = pd.read_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types.tsv", sep="\t", index_col=0, header=[0,1,2])
    fs_combined.columns.names = ["Env", "Data", "Info"]
    fs_combined_50 = pd.read_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types_top50.tsv", sep="\t", index_col=0, header=[0,1,2])
    fs_combined_50.columns.names = ["Env", "Data", "Info"]

    # Heatmaps of all top genes and top 50 genes in each environment
    sub = fs_combined.loc[:,fs_combined.columns.get_level_values("Info")=="rank_per_bin"]
    sub.replace({"(0.99, 1]":4, "(0.95, 0.99]":3, "(0.9, 0.95]":2, "(0, 0.90]":1}, inplace=True)
    sub.fillna(0, inplace=True)
    cmap = ListedColormap(["whitesmoke", "gold", "orange", "red", "darkred"])
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 8.5), sharex=True)
    sns.heatmap(sub.loc[:,sub.columns.get_level_values("Data")=="snp"], cmap=cmap, xticklabels=False, yticklabels=False, ax=ax[0][0])
    sns.heatmap(sub.loc[:,sub.columns.get_level_values("Data")=="orf"], cmap=cmap, xticklabels=False, yticklabels=False, ax=ax[0][1])
    sns.heatmap(sub.loc[:,sub.columns.get_level_values("Data")=="cnv"], cmap=cmap, xticklabels=False, yticklabels=False, ax=ax[1][0])
    plt.savefig("Scripts/Data_Vis/Figure2a.pdf")
    plt.close()

    sub50 = fs_combined_50.loc[:,fs_combined_50.columns.get_level_values("Info")=="rank_per_bin"]
    sub50.replace({"(0.99, 1]":4, "(0.95, 0.99]":3, "(0.9, 0.95]":2, "(0, 0.90]":1}, inplace=True)
    sub50.fillna(0, inplace=True)
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 8.5), sharex=True)
    sns.heatmap(sub50.loc[:,sub50.columns.get_level_values("Data")=="snp"], cmap=cmap, xticklabels=False, yticklabels=False, ax=ax[0][0])
    sns.heatmap(sub50.loc[:,sub50.columns.get_level_values("Data")=="orf"], cmap=cmap, xticklabels=False, yticklabels=False, ax=ax[0][1])
    sns.heatmap(sub50.loc[:,sub50.columns.get_level_values("Data")=="cnv"], cmap=cmap, xticklabels=False, yticklabels=False, ax=ax[1][0])
    plt.savefig("Scripts/Data_Vis/Figure2a_top50.pdf")
    plt.close()

    ############################################################################
    # Figure 2
    ############################################################################
    pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # fitness data
    h2 = pd.read_csv("Results/Heritability_h2_H2_sommer.csv", index_col=0) # narrow-sense heritability
    
    # Single-environment models
    pcs = pd.read_csv("Results/RESULTS_RF_PCs_sorted.txt", sep="\t", index_col=0)
    snp = pd.read_csv("Results/RESULTS_RF_SNPs_FS.txt", sep="\t", index_col=0)
    orf = pd.read_csv("Results/RESULTS_RF_ORFs_FS.txt", sep="\t", index_col=0)
    cnv = pd.read_csv("Results/RESULTS_RF_CNVs_FS.txt", sep="\t", index_col=0)
    single = pd.concat([pcs.r2_test, snp.r2_test, orf.r2_test, cnv.r2_test], ignore_index=False, axis=1)
    single.columns = ["PCs", "SNPs", "PAVs", "CNVs"]
    single.to_csv("Results/RESULTS_RF_ALL_FS.txt", sep="\t")

    # Multi-environment models
    brr_res = pd.read_csv("/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/RESULTS_GxE.txt", sep="\t")
    brr_res["ID"] = brr_res.apply(lambda x: x.ID.split("_")[1], axis=1)
    brr_res = brr_res.loc[brr_res.ID==brr_res.Env,:]
    dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/multi-trait/cv1"
    mo_res =  [f"{dir}/{file}" for file in os.listdir(dir) if file.endswith("_preds_test_multitrait.txt")]
    multi = brr_res[["Alg", "Env", "r2_test"]]
    from sklearn.metrics import r2_score
    for file in mo_res:
        Env =  file.split("_")[-8]
        df = pd.read_csv(file, sep="\t", index_col=[0,1])
        means = df.groupby("Trait").mean().T
        r2 = r2_score(pheno.loc[means.index,Env], means[Env])
        toplot = toplot._append({"Alg":"Multi-output", "Env":Env, "r2_test":r2}, ignore_index=True)
    multi.to_csv("Results/RESULTS_ALL_ALG_GxE_PAVs.txt", sep="\t", index=False)

    #### Plot Single- and Multi-environment models results
    single = pd.read_csv("Results/RESULTS_RF_ALL_FS.txt", sep="\t", index_col=0) # single-environment models
    multi = pd.read_csv("Results/RESULTS_ALL_ALG_GxE_PAVs.txt", sep="\t") # multi-environment models
    meta = pd.concat([pheno.median(), h2["h2"]], ignore_index=False, axis=1) # median fitness and narrow-sense heritability
    media = pd.Series([v[0] for v in pheno.columns.str.split(r"(?<=YP(?!D))|(?<=YPD).+", regex=True)], index=pheno.columns) # growth media
    envs = [v[1] for v in pheno.columns.str.split(r"^YP(?!D)|YPD", regex=True)] # environments
    envs = {pheno.columns[i]:envs[i] for i in range(len(envs))}
    meta = pd.concat([media, meta], ignore_index=False, axis=1)
    meta.columns = ["media", "median", "h2"]
    
    # cluster single-environment models and reorder meta and multi accordingly
    dend = hierarchy.dendrogram(hierarchy.linkage(single, method="average"))
    single = single.iloc[dend["leaves"],:]
    multi = multi.pivot(index="Env", columns="Alg", values="r2_test") # reshape multi-env matrix
    multi = multi.loc[single.index,:]
    meta = meta.loc[single.index,:]

    # plot heatmaps
    single.rename(index=envs, inplace=True)
    multi.rename(index=envs, inplace=True)
    meta.rename(index=envs, inplace=True)
    meta.media.replace({"YPD":0, "YP":1}, inplace=True) # recode meta.media to binary
    fig, ax = plt.subplots(1, 4, figsize=(10, 11), sharey=True)
    sns.heatmap(meta.iloc[:,:2], cmap="plasma", ax=ax[0], cbar=True, cbar_kws={"location":"bottom"}, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    sns.heatmap(meta.iloc[:,2:], cmap="magma", ax=ax[1], cbar=True, cbar_kws={"location":"bottom"}, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    sns.heatmap(single, cmap="viridis", ax=ax[2], cbar=True, cbar_kws={"location":"bottom"}, center=0.3, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    sns.heatmap(multi, cmap="viridis", ax=ax[3], cbar=True, cbar_kws={"location":"bottom"}, center=0.3, xticklabels=True, yticklabels=True, annot=True, annot_kws={"fontsize": 6})
    plt.tight_layout()
    plt.savefig("Scripts/Data_Vis/Figure_2_model_performances.pdf")
    plt.close()

    ############################################################################
    # Figure 2
    ############################################################################
    fs_combined = pd.read_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types.tsv", sep="\t", index_col=0, header=[0,1,2])
    fs_combined.columns.names = ["Env", "Data", "Info"]
    fs_combined_50 = pd.read_csv("Scripts/Data_Vis/RF_FS_imp_all_data_types_top50.tsv", sep="\t", index_col=0, header=[0,1,2])
    fs_combined_50.columns.names = ["Env", "Data", "Info"]
    
    # Replace gene rank bin categories with 1s and missing values with 0s
    sub = fs_combined.loc[:,fs_combined.columns.get_level_values("Info")=="rank_per_bin"] # all the top genes in each environment based on feature selection
    sub = sub.notnull().astype("int")
    sub50 = fs_combined_50.loc[:,fs_combined_50.columns.get_level_values("Info")=="rank_per_bin"] # top 50 genes in each environmnet
    sub50 = sub50.notnull().astype("int")
    sub5p = {} # top 5% of genes in each environment
    for col in fs_combined.columns:
        if "rank_per" in col:
            idx = fs_combined.loc[fs_combined[col]<=0.05, col].index # top 5%
            sub5p[col] = fs_combined.loc[idx, (*col[:2], "rank_per_bin")]
    sub5p = pd.DataFrame(sub5p)
    sub5p.columns.names = ["Env", "Data", "Info"]
    sub5p = sub5p.notnull().astype("int")

    ## 2b. distribution of how many genes are shared in n number of environments
    fig1, ax1 = plt.subplots(nrows=3, ncols=3, figsize=(12, 11.25)) # logged x-axis
    fig2, ax2 = plt.subplots(nrows=3, ncols=3, figsize=(12, 11.25)) # env 1-35 x-axis
    fig3, ax3 = plt.subplots(nrows=3, ncols=3, figsize=(12, 11.25)) # logged x-axis; restricted to 1-20 envs
    fig4, ax4 = plt.subplots(nrows=3, ncols=3, figsize=(12, 11.25)) # env 1-35 x-axis; restricted to 1-20 envs
    ks_res_df1, ks_res_gps_df1, ks_rand_vs_gps_df1, ax1, ax2, ax3, ax4 = make_fig2a(sub, ("snp", "orf", "cnv"), ax1, ax2, ax3, ax4, 0, "all top genes") # all top genes based on feature selection
    ks_res_df2, ks_res_gps_df2, ks_rand_vs_gps_df2, ax1, ax2, ax3, ax4 = make_fig2a(sub50, ("snp", "orf", "cnv"), ax1, ax2, ax3, ax4, 1, "top 50 genes") # top 50 genes
    ks_res_df3, ks_res_gps_df3, ks_rand_vs_gps_df3, ax1, ax2, ax3, ax4 = make_fig2a(sub5p, ("snp", "orf", "cnv"), ax1, ax2, ax3, ax4, 2, "top 5% genes") # top 5% genes
    fig1.tight_layout(); fig2.tight_layout(); fig3.tight_layout(); fig4.tight_layout()
    fig1.savefig("Scripts/Data_Vis/Figure2b_logged_x-axis_v2.pdf") # _v2 added bc median line was plotted
    fig2.savefig("Scripts/Data_Vis/Figure2b_v2.pdf")
    fig3.savefig("Scripts/Data_Vis/Figure2b_logged_x-axis_20envs_v2.pdf")
    fig4.savefig("Scripts/Data_Vis/Figure2b_20envs_v2.pdf")
    plt.close(); plt.close(); plt.close(); plt.close()
    ks_res_df1.to_csv("Scripts/Data_Vis/Table_S_KS_res_fig2b_all_top_genes.csv", index=False)
    ks_res_df2.to_csv("Scripts/Data_Vis/Table_S_KS_res_grouped_fig2b_top50_genes.csv", index=False)
    ks_res_df3.to_csv("Scripts/Data_Vis/Table_S_KS_res_nongrouped_vs_grouped_fig2b_top5p_genes.csv", index=False)
    ks_res_gps_df1.to_csv("Scripts/Data_Vis/Table_S_KS_res_fig2b_all_top_genes.csv", index=False)
    ks_res_gps_df2.to_csv("Scripts/Data_Vis/Table_S_KS_res_grouped_fig2b_top50_genes.csv", index=False)
    ks_res_gps_df3.to_csv("Scripts/Data_Vis/Table_S_KS_res_nongrouped_vs_grouped_fig2b_top5p_genes.csv", index=False)
    ks_rand_vs_gps_df1.to_csv("Scripts/Data_Vis/Table_S_KS_res_fig2b_all_top_genes.csv", index=False)
    ks_rand_vs_gps_df2.to_csv("Scripts/Data_Vis/Table_S_KS_res_grouped_fig2b_top50_genes.csv", index=False)
    ks_rand_vs_gps_df3.to_csv("Scripts/Data_Vis/Table_S_KS_res_nongrouped_vs_grouped_fig2b_top5p_genes.csv", index=False)
    
    ## 2c. orf importance vs percent presence density plots and plots colored by gene copy number
    # Read presence/absence variation data and calculate percent orf presence in the population
    pav = pd.read_csv("Data/Peter_2018/ORFs_pres_abs.csv", index_col=0)
    percent_presence = pd.DataFrame((pav.sum()/pav.shape[0])*100, index=pav.columns).reset_index()
    percent_presence.columns = ["orf", "percent_presence"]
    percent_presence.set_index("orf", inplace=True)
    percent_presence.percent_presence = percent_presence.percent_presence.astype("float64")
    
    # Read in orf (gene) copy number data
    cnv = pd.read_csv("~/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv", index_col=0)
    cnv_stats = cnv.describe().T
    
    # Combine orf importance scores data (RF baseline, all ORF pres/abs features)
    dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline"
    pav_baseline_files = [os.path.join(dir, f"{x}_orf_baseline_imp") for x in mapping.keys()]
    imp = []
    for i in range(len(pav_baseline_files)):
        df = pd.read_csv(pav_baseline_files[i], sep="\t", index_col=0)
        imp.append(df)
    imp = pd.concat(imp, axis=1)
    imp.columns = [match.group() for s in pav_baseline_files for match in [re.search("[A-Z0-9]+(?=_orf)", s)] if match] # get column name that matches file name
    imp = imp.loc[percent_presence.index, :] # sort by percent presence order

    # Plot gene importance vs percent presence scatter plots colored by copy number
    fig, axes = plt.subplots(7, 5, figsize=(25, 28))
    for i, col in enumerate(imp.columns):
        row_idx, col_idx = divmod(i, 5)
        toplot = percent_presence.merge(imp[col], how="left", on="orf")
        toplot = toplot.merge(cnv_stats, how="left", left_on="orf", right_index=True)
        sns.scatterplot(data=toplot, x="percent_presence", y=col, c=cm.plasma(toplot["mean"]), alpha=0.3, ax=axes[row_idx, col_idx])
        norm = plt.Normalize(toplot["mean"].min(), toplot["mean"].max()) # normalize mean gene copy number values
        cbar = plt.cm.ScalarMappable(norm=norm, cmap="plasma") # colorbar
        axes[row_idx, col_idx].figure.colorbar(cbar, ax=axes[row_idx, col_idx]) # set colorbar
        axes[row_idx, col_idx].set_title(col)
        del toplot
    plt.tight_layout()
    plt.savefig(f'Scripts/Data_Vis/Figure2d.pdf')
    plt.close()

    # Plot density scatter plots
    res = []
    fig, axes = plt.subplots(7, 5, figsize=(15, 15)) # gene imp vs % presence density scatter plots
    for i, col in enumerate(imp.columns):
        # Density scatter plots
        spearman_corr, p = spearmanr(percent_presence["percent_presence"], imp[col]) # spearman's rho and p-value
        values = np.vstack([percent_presence["percent_presence"], imp[col]])
        kernel = stats.gaussian_kde(values)(values) # point kernel density
        row_idx, col_idx = divmod(i, 5)  # row and column index for subplots
        ax = axes[row_idx, col_idx] # set subplot
        sns.scatterplot(x=percent_presence["percent_presence"], y=imp[col],
                        c=kernel, edgecolor = 'none', size=-5, alpha=.3,
                        cmap="viridis", legend=False, ax=ax) # density scatter plot
        ax.set_title(f'{mapping[col]}')
        norm = plt.Normalize(kernel.min(), kernel.max()) # normalize kde values
        cbar = plt.cm.ScalarMappable(norm=norm, cmap="viridis") # colorbar
        ax.figure.colorbar(cbar, ax=ax) # set colorbar
        m, b, r, p_val, std_err = linregress(percent_presence["percent_presence"], imp[col]) # linear fit
        ax.plot(percent_presence["percent_presence"], m*percent_presence["percent_presence"] + b) # plot line
        res.append([col, spearman_corr, p, m, b, r, p_val, std_err]) # save stats results
        del spearman_corr, p, values, kernel, row_idx, col_idx, ax, norm, cbar, m, b, r, p_val, std_err
    plt.tight_layout()
    plt.savefig("Scripts/Data_Vis/Figure2c_baseline.pdf")
    plt.close()
    res = pd.DataFrame(res, columns=["Environment", "Spearman's rho", "p-value",
                                     "lm slope", "lm intercept", "lm pearson r",
                                     "lm p-value", "lm std_err"])
    res.to_csv("Scripts/Data_Vis/Table_S_baseline_importance_vs_percent_presence.csv", index=False)

    # Read in orf importance scores data for RF FS models
    dir = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs"
    pav_rf_res = pd.read_csv("Results/RESULTS_RF_ORFs_FS.txt", sep="\t")  # ORF pres/abs FS results
    pav_files = [os.path.join(dir, f"{x}_imp") for x in pav_rf_res['ID']]

    # Plot density scatter plots
    res = []
    fig, axes = plt.subplots(7, 5, figsize=(15, 15))
    fig2, axes2 = plt.subplots(7, 5, figsize=(15, 15))
    for i, col in enumerate(mapping.keys()):
        file = [f for f in pav_files if col in f] # Read feature importance score files
        imp = pd.read_csv(file[0], sep="\t", skiprows=1, names=["orf", "mean_imp"], index_col=0)
        imp = imp.merge(percent_presence, how="left", on="orf") # merge with percent_presence data
        imp = imp.merge(cnv_stats, how="left", left_on="orf", right_index=True) # merge with cnv_stats data
        spearman_corr, p = spearmanr(imp["percent_presence"], imp["mean_imp"]) # spearman's rho and p-value
        values = np.vstack([imp["percent_presence"], imp["mean_imp"]])
        kernel = stats.gaussian_kde(values)(values) # point kernel density
        row_idx, col_idx = divmod(i, 5)  # row and column index for subplots
        ax = axes[row_idx, col_idx] # set subplot
        ax2 = axes2[row_idx, col_idx] # set subplot
        # color by point density
        sns.scatterplot(x=imp["percent_presence"], y=imp["mean_imp"],
                        c=kernel, edgecolor = 'none', size=-5, alpha=.3,
                        cmap="viridis", legend=False, ax=ax) # density scatter plot
        ax.set_title(f"{mapping[col]}")
        norm = plt.Normalize(kernel.min(), kernel.max()) # normalize kde values
        cbar = plt.cm.ScalarMappable(norm=norm, cmap="viridis") # colorbar
        ax.figure.colorbar(cbar, ax=ax) # set colorbar
        m, b, r, p_val, std_err = linregress(imp["percent_presence"], imp["mean_imp"]) # linear fit
        ax.plot(imp["percent_presence"], m*imp["percent_presence"] + b) # plot line
        res.append([col, spearman_corr, p, m, b, r, p_val, std_err]) # save stats results
        # sizes by gene copy number
        sns.scatterplot(x=imp["percent_presence"], y=imp["mean_imp"], c=np.log10(imp["mean"]),
                        size=np.log10(imp["mean"]), edgecolor = 'gray', alpha=.15,
                        cmap="plasma", legend=True, ax=ax2) # density scatter plot
        ax2.set_title(f"{mapping[col]}")
        norm = plt.Normalize(np.log10(imp["mean"]).min(), np.log10(imp["mean"]).max()) # normalize gene copy number values
        cbar = plt.cm.ScalarMappable(norm=norm, cmap="plasma") # colorbar
        ax2.figure.colorbar(cbar, ax=ax2) # set colorbar
        ax2.plot(imp["percent_presence"], m*imp["percent_presence"] + b) # plot line
        del spearman_corr, p, values, kernel, row_idx, col_idx, ax, norm, cbar, m, b, r, p_val, std_err
    fig.tight_layout()
    fig.savefig("Scripts/Data_Vis/Figure2c_fs.pdf")
    fig2.tight_layout()
    fig2.savefig(f"Scripts/Data_Vis/Figure2d_fs.pdf")
    plt.close(); plt.close()
    res = pd.DataFrame(res, columns=["Environment", "Spearman's rho", "p-value",
                                     "lm slope", "lm intercept", "lm pearson r",
                                     "lm p-value", "lm std_err"])
    res.to_csv("Scripts/Data_Vis/Table_S_fs_importance_vs_percent_presence.csv", index=False)

    ############################################################################
    # Figure 3
    ############################################################################
    ## 3a. gene rank comparison between data types density plots
    baseline_combined = pd.read_csv("Scripts/Data_Vis/RF_baseline_imp_all_data_types.tsv", sep="\t", index_col=0, header=[0,1,2])
    baseline_combined.columns.names = ["Env", "Data", "Info"]
    res = []
    fig_svo, axes_svo = plt.subplots(7, 5, figsize=(20, 19))
    fig_svc, axes_svc = plt.subplots(7, 5, figsize=(20, 19))
    fig_ovc, axes_ovc = plt.subplots(7, 5, figsize=(20, 19))
    for i,env in tqdm(enumerate(mapping.keys()), total=len(mapping.keys())):
        # subset dataframe by data type and data type pairs
        snp_rank = baseline_combined.loc[:,(baseline_combined.columns.get_level_values("Env")==env) &
                                         (baseline_combined.columns.get_level_values("Data")=="snp")]
        pav_rank = baseline_combined.loc[:,(baseline_combined.columns.get_level_values("Env")==env) &
                                         (baseline_combined.columns.get_level_values("Data")=="orf")]
        cnv_rank = baseline_combined.loc[:,(baseline_combined.columns.get_level_values("Env")==env) &
                                         (baseline_combined.columns.get_level_values("Data")=="cnv")]
        # merge to ensure only common genes are used to calculate correlation
        toplot_svo = snp_rank.dropna(how="all").merge(pav_rank.dropna(how="all"), how="inner", on="gene")
        toplot_svc = snp_rank.dropna(how="all").merge(cnv_rank.dropna(how="all"), how="inner", on="gene")
        toplot_ovc = pav_rank.dropna(how="all").merge(cnv_rank.dropna(how="all"), how="inner", on="gene")
        toplot_svo.sort_values(by=(env, "snp", "rank_per"), inplace=True)
        toplot_svc.sort_values(by=(env, "snp", "rank_per"), inplace=True)
        toplot_ovc.sort_values(by=(env, "orf", "rank_per"), inplace=True)
        # calculate spearman's rho
        corr_svo, p_svo = spearmanr(toplot_svo[(env, "snp", "rank_per")], toplot_svo[(env, "orf", "rank_per")])
        corr_svc, p_svc = spearmanr(toplot_svc[(env, "snp", "rank_per")], toplot_svc[(env, "cnv", "rank_per")])
        corr_ovc, p_ovc = spearmanr(toplot_ovc[(env, "orf", "rank_per")], toplot_ovc[(env, "cnv", "rank_per")])
        # linear fit
        m1, b1, r1, p1, stde1 = linregress(toplot_svo[(env, "snp", "rank_per")], toplot_svo[(env, "orf", "rank_per")])
        m2, b2, r2, p2, stde2 = linregress(toplot_svc[(env, "snp", "rank_per")], toplot_svc[(env, "cnv", "rank_per")])
        m3, b3, r3, p3, stde3 = linregress(toplot_ovc[(env, "orf", "rank_per")], toplot_ovc[(env, "cnv", "rank_per")])
        res.append([env, "SNP vs ORF", corr_svo, p_svo, m1, b1, r1, p1, stde1])
        res.append([env, "SNP vs CNV", corr_svc, p_svc, m2, b2, r2, p2, stde2])
        res.append([env, "ORF vs CNV", corr_ovc, p_ovc, m3, b3, r3, p3, stde3])
        # density contour plot
        row_idx, col_idx = divmod(i, 5)  # row and column index for subplots
        sns.kdeplot(x=toplot_svo[env]["snp"]["rank_per"], y=toplot_svo[env]["orf"]["rank_per"],
                    cmap="viridis", fill=True, cbar=True, bw_adjust=.5, ax=axes_svo[row_idx, col_idx])
        axes_svo[row_idx, col_idx].set_title(f'{mapping[env]}')
        sns.kdeplot(x=toplot_svc[env]["snp"]["rank_per"], y=toplot_svc[env]["cnv"]["rank_per"],
                    cmap="viridis", fill=True, cbar=True, bw_adjust=.5, ax=axes_svc[row_idx, col_idx])
        axes_svc[row_idx, col_idx].set_title(f'{mapping[env]}')
        sns.kdeplot(x=toplot_ovc[(env, "orf", "rank_per")], y=toplot_ovc[(env, "cnv", "rank_per")],
                    cmap="viridis", fill=True, cbar=True, bw_adjust=.5, ax=axes_ovc[row_idx, col_idx])
        axes_ovc[row_idx, col_idx].set_title(f'{mapping[env]}')
    fig_svo.tight_layout()
    fig_svo.savefig(f"Scripts/Data_Vis/Figure3_snp_vs_pav.pdf")
    fig_svc.tight_layout()
    fig_svc.savefig(f"Scripts/Data_Vis/Figure3_snp_vs_cnv.pdf")
    fig_ovc.tight_layout()
    fig_ovc.savefig(f"Scripts/Data_Vis/Figure3_pav_vs_cnv.pdf")
    plt.close(); plt.close(); plt.close()
    res = pd.DataFrame(res, columns=["Environment", "Data Type Pair", "Spearman's rho",
                                     "p-value", "lm slope", "lm intercept", "lm pearson r",
                                     "lm p-value", "lm std_err"])
    res.sort_values(by="Spearman's rho", ascending=False, inplace=True)
    res.to_csv("Scripts/Data_Vis/Table_S_baseline_rank_corr_between_data_types.csv", index=False) # Table S5
    
    # 3b. Gene rank comparison between environments per data type


    ############################################################################
    # Figure 4: GO Enrichment
    ############################################################################
    # Read excel spreadsheets
    snp_go = pd.read_excel("Scripts/Data_Vis/SNP_Figures/SNP_GO_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
    snp_go.columns = ["Environment", "GO", "Gene_top_has_GO", "Gene_not_top_has_GO",
                      "Gene_top_no_GO", "Gene_not_top_no_GO", "direction",
                      "p.val", "odds ratio", "qvalues",  "BP", "CC", "MF"]
    snp_go["Data"] = "SNP"
    pav_go = pd.read_excel("Scripts/Data_Vis/ORF_CNV_Figures/ORF_GO_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
    pav_go.columns = ["Environment", "GO", "Gene_top_has_GO", "Gene_not_top_has_GO",
                      "Gene_top_no_GO", "Gene_not_top_no_GO", "direction",
                      "p.val", "odds ratio", "qvalues",  "BP", "CC", "MF"]
    pav_go["Data"] = "ORF"
    cnv_go = pd.read_excel("Scripts/Data_Vis/ORF_CNV_Figures/CNV_GO_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
    cnv_go.columns = ["Environment", "GO", "Gene_top_has_GO", "Gene_not_top_has_GO",
                      "Gene_top_no_GO", "Gene_not_top_no_GO", "direction",
                      "p.val", "odds ratio", "qvalues",  "BP", "CC", "MF"]
    cnv_go["Data"] = "CNV"
    
    # Combine
    go_combined = pd.concat([snp_go, pav_go, cnv_go], axis=0)
    go_combined = go_combined[go_combined["qvalues"]<=0.05] # filter for significant GO terms
    go_combined = go_combined[["Data", "Environment", "BP", "CC", "MF", "qvalues", "odds ratio", "direction"]]
    go_combined["log10(odds ratio)"] = np.log10(go_combined["odds ratio"]) # log10 transform odds ratio
    go_combined.loc[np.isinf(go_combined["odds ratio"]),"log10(odds ratio)"] = 1 # set infinite odds ratios to 1
    go_combined = go_combined.melt(id_vars=["Data", "Environment", "log10(odds ratio)"],
                                   value_vars=["BP", "CC", "MF"], var_name="GO Type",
                                   value_name="GO Description") # melt GO related columns
    go_combined.dropna(inplace=True) # drop rows with missing values
    go_combined.set_index(keys=["GO Type", "GO Description"], inplace=True) # create multi-index
    go_combined = go_combined.pivot(columns=["Data", "Environment"]) # unstack multi-index
    go_combined.to_csv("Scripts/Data_Vis/GO_enrichment_all.txt", sep="\t")
    
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
    plt.savefig("Scripts/Data_Vis/GO_enrichment_all_heatmap_v2.pdf")
    plt.close()

    ############################################################################
    # Figure 5: Pathway Enrichment
    ############################################################################

    # Read excel spreadsheets
    snp_pwy = pd.read_excel("Scripts/Data_Vis/SNP_Figures/SNP_pathway_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
    snp_pwy.columns = ["Environment", "PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY",
                      "Gene_top_no_PWY", "Gene_not_top_no_PWY", "direction",
                      "p.val", "odds ratio", "qvalues",  "PWY name", "PWY ID"]
    snp_pwy["Data"] = "SNP"
    pav_pwy = pd.read_excel("Scripts/Data_Vis/ORF_CNV_Figures/ORF_pathway_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
    pav_pwy.columns = ["Environment", "PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY",
                      "Gene_top_no_PWY", "Gene_not_top_no_PWY", "direction",
                      "p.val", "odds ratio", "qvalues",  "PWY name", "PWY ID"]
    pav_pwy["Data"] = "ORF"
    cnv_pwy = pd.read_excel("Scripts/Data_Vis/ORF_CNV_Figures/CNV_pathway_enrichment.xlsx", sheet_name="All").reset_index(drop=True)
    cnv_pwy.columns = ["Environment", "PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY",
                      "Gene_top_no_PWY", "Gene_not_top_no_PWY", "direction",
                      "p.val", "odds ratio", "qvalues",  "PWY name", "PWY ID"]
    cnv_pwy["Data"] = "CNV"
    
    # Combine
    pwy_combined = pd.concat([snp_pwy, pav_pwy, cnv_pwy], axis=0)
    pwy_combined = pwy_combined[pwy_combined["qvalues"]<=0.05] # filter for significant pathways, None were found

    ############################################################################
    # Figure ?: SHAP values of top genes VS genetic similarity to S288C
    ############################################################################
    test = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", header=None)
    geno = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv").to_pandas()
    geno.set_index("ID", inplace=True)
    
    # Add S288C genotypes as a row of -1s
    geno = geno.T
    geno["S288C"] = -1
    geno = geno.T
    geno.to_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_with_S288C.csv")

    # Calculate the euclidean distance and the jaccard distance
    from scipy.spatial.distance import pdist, squareform
    geno_train = geno.loc[~geno.index.isin(test[0]),:]
    eu_dist = pdist(geno_train.values, metric="euclidean")
    eu_dist = pd.DataFrame(squareform(eu_dist), columns=geno_train.index, index=geno_train.index)
    ja_dist = pdist(geno_train.values, metric="jaccard") # proportion of dissimilarity
    ja_dist = pd.DataFrame(squareform(ja_dist), columns=geno_train.index, index=geno_train.index)

    # Get the most similar strain to S288C per clade
    clades = pd.read_excel("Data/Peter_2018/Peter_2018_Supplementary_Tables.xls", sheet_name="Table S1", skiprows=3, nrows=1011) # isolate clades
    clades = clades[["Standardized name", "Clades"]] # subset relevant columns
    clades = clades.loc[clades["Standardized name"].isin(geno_train.index),:] # only training instances
    clades.fillna("Unknown", inplace=True) # some instances did not fit into a clade

    ja_dist.loc["S288C"].iloc[:625].idxmin() # SACE_YAF most similar to S288C based on jaccard distance
    eu_dist.loc["S288C"].iloc[:625].idxmin() # CHF most similar to S288C based on euclidean distance
    ja_dist.loc["S288C"].loc[['SACE_YAF','CHF']] # looking at the dissimilarity
    eu_dist.loc["S288C"].loc[['SACE_YAF','CHF']]

    res = []
    for clade in clades.Clades.unique():
        print(clade)
        clades_sub = clades.loc[clades.Clades==clade, :]
        res.append([clade,\
        ja_dist.loc[ja_dist.index.isin(clades_sub["Standardized name"]), "S288C"].idxmin(),\
        ja_dist.loc[ja_dist.index.isin(clades_sub["Standardized name"]), "S288C"].min(),\
        eu_dist.loc[eu_dist.index.isin(clades_sub["Standardized name"]), "S288C"].idxmin(),\
        eu_dist.loc[eu_dist.index.isin(clades_sub["Standardized name"]), "S288C"].min()])
    
    res = pd.DataFrame(res, columns=["clade", "most similar by Jaccard",\
        "jaccard distance", "most similar by euclidean", "euclidean distance"])

    # Assign colors to each clade for plotting
    cmap = cm.get_cmap('hsv', clades.Clades.nunique()) # sample hsv color palette
    color_list = {clades.Clades.unique()[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
    clades["colors"] = clades['Clades'].map(color_list) # assign colors to isolates
    clades.set_index("Standardized name", inplace=True)

    # Plot the top feature shap value and similarity value for the representative isolate for each clade
    top_features = pd.read_csv("Scripts/Data_Vis/Section_5/shap_value_top_features_per_data_type.csv")
    target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
                "YPDSODIUMMETAARSENITE"]
    d = "/mnt/gs21/scratch/seguraab/yeast_project/SHAP_Interaction"
    for env in target_envs:
        # read in the shap values
        snp_shap_all = pd.read_csv(glob.glob(f"{d}/SNP/SHAP_values_sorted_snp_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
        pav_shap_all = pd.read_csv(glob.glob(f"{d}/PAV/SHAP_values_sorted_pav_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
        cnv_shap_all = pd.read_csv(glob.glob(f"{d}/CNV/SHAP_values_sorted_cnv_{env}_top_*_plus_comb_lit_genes_training.txt")[0], sep="\t", index_col=0)
        pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix orf IDs
        pav_shap_all.columns = pav_shap_all.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
        cnv_shap_all.columns = cnv_shap_all.apply(lambda x: re.sub("^X", "", x.name), axis=0)
        cnv_shap_all.columns = cnv_shap_all.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
        ### Violin plot
        sns.violinplot(data=pd.concat([clades.colors, snp_shap_all], ignore_index=False, axis=1),
            y=top_features.loc[top_features["Unnamed: 0"]==env, "snp"].values[0],
            hue="colors")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_5/{env}_top_gene_shap_vs_clade_violin.pdf")

        ### Plot all isolates, regardless of clade, against the SHAP values of the top feature
        fig, ax = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(10,3))
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "snp"].values[0]
        y = snp_shap_all.loc[snp_shap_all.index.isin(res["most similar by Jaccard"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "snp"].values[0]]
        data = pd.concat([res.set_index("most similar by Jaccard").loc[:,"jaccard distance"], y], axis=1)
        sns.regplot(x=data["jaccard distance"], y=data[feat], ci=None, ax=ax[0])
        m, b, r, p, sd = stats.linregress(x=data["jaccard distance"], y=data[feat])
        ax[0].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[0].set_xlabel("Jaccard Distance with S288C (one isolate per clade)")
        ax[0].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'snp'].values[0]}")
        # now PAVs
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "pav"].values[0]
        y = pav_shap_all.loc[pav_shap_all.index.isin(res["most similar by Jaccard"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "pav"].values[0]]
        data = pd.concat([res.set_index("most similar by Jaccard").loc[:,"jaccard distance"], y], axis=1)
        sns.regplot(x=data["jaccard distance"], y=data[feat], ci=None, ax=ax[1])
        m, b, r, p, sd = stats.linregress(x=data["jaccard distance"], y=data[feat])
        ax[1].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[1].set_xlabel("Jaccard Distance with S288C (one isolate per clade)")
        ax[1].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'pav'].values[0]}")
        # now CNVs
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "cnv"].values[0]
        y = cnv_shap_all.loc[cnv_shap_all.index.isin(res["most similar by Jaccard"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "cnv"].values[0]]
        data = pd.concat([res.set_index("most similar by Jaccard").loc[:,"jaccard distance"], y], axis=1)
        sns.regplot(x=data["jaccard distance"], y=data[feat], ci=None, ax=ax[2])
        m, b, r, p, sd = stats.linregress(x=data["jaccard distance"], y=data[feat])
        ax[2].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[2].set_xlabel("Jaccard Distance with S288C (one isolate per clade)")
        ax[2].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'cnv'].values[0]}")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_5/{env}_top_gene_shap_vs_clade_jaccard_distance_to_S288C.pdf")
        plt.close()
        #
        ### Plot the linear regression plot for jaccard distance first for SNP-based SHAP values; one isolate per clade
        fig, ax = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(10,3))
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "snp"].values[0]
        y = snp_shap_all.loc[snp_shap_all.index.isin(res["most similar by Jaccard"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "snp"].values[0]]
        data = pd.concat([res.set_index("most similar by Jaccard").loc[:,"jaccard distance"], y], axis=1)
        sns.regplot(x=data["jaccard distance"], y=data[feat], ci=None, ax=ax[0])
        m, b, r, p, sd = stats.linregress(x=data["jaccard distance"], y=data[feat])
        ax[0].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[0].set_xlabel("Jaccard Distance with S288C (one isolate per clade)")
        ax[0].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'snp'].values[0]}")
        # now PAVs
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "pav"].values[0]
        y = pav_shap_all.loc[pav_shap_all.index.isin(res["most similar by Jaccard"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "pav"].values[0]]
        data = pd.concat([res.set_index("most similar by Jaccard").loc[:,"jaccard distance"], y], axis=1)
        sns.regplot(x=data["jaccard distance"], y=data[feat], ci=None, ax=ax[1])
        m, b, r, p, sd = stats.linregress(x=data["jaccard distance"], y=data[feat])
        ax[1].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[1].set_xlabel("Jaccard Distance with S288C (one isolate per clade)")
        ax[1].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'pav'].values[0]}")
        # now CNVs
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "cnv"].values[0]
        y = cnv_shap_all.loc[cnv_shap_all.index.isin(res["most similar by Jaccard"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "cnv"].values[0]]
        data = pd.concat([res.set_index("most similar by Jaccard").loc[:,"jaccard distance"], y], axis=1)
        sns.regplot(x=data["jaccard distance"], y=data[feat], ci=None, ax=ax[2])
        m, b, r, p, sd = stats.linregress(x=data["jaccard distance"], y=data[feat])
        ax[2].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[2].set_xlabel("Jaccard Distance with S288C (one isolate per clade)")
        ax[2].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'cnv'].values[0]}")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_5/{env}_top_gene_shap_vs_clade_jaccard_distance_to_S288C.pdf")
        plt.close()
        #
        ### Plot the linear regression plot for euclidean distance first for SNP-based SHAP values; one isolate per clade
        fig, ax = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(10,3))
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "snp"].values[0]
        y = snp_shap_all.loc[snp_shap_all.index.isin(res["most similar by euclidean"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "snp"].values[0]]
        data = pd.concat([res.set_index("most similar by euclidean").loc[:,"euclidean distance"], y], axis=1)
        sns.regplot(x=data["euclidean distance"], y=data[feat], ci=None, ax=ax[0])
        m, b, r, p, sd = stats.linregress(x=data["euclidean distance"], y=data[feat])
        ax[0].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[0].set_xlabel("Euclidean Distance with S288C (one isolate per clade)")
        ax[0].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'snp'].values[0]}")
        # now PAVs
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "pav"].values[0]
        y = pav_shap_all.loc[pav_shap_all.index.isin(res["most similar by euclidean"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "pav"].values[0]]
        data = pd.concat([res.set_index("most similar by euclidean").loc[:,"euclidean distance"], y], axis=1)
        sns.regplot(x=data["euclidean distance"], y=data[feat], ci=None, ax=ax[1])
        m, b, r, p, sd = stats.linregress(x=data["euclidean distance"], y=data[feat])
        ax[1].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[1].set_xlabel("Euclidean Distance with S288C (one isolate per clade)")
        ax[1].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'pav'].values[0]}")
        # now CNVs
        feat = top_features.loc[top_features["Unnamed: 0"]==env, "cnv"].values[0]
        y = cnv_shap_all.loc[cnv_shap_all.index.isin(res["most similar by euclidean"]),\
                             top_features.loc[top_features["Unnamed: 0"]==env, "cnv"].values[0]]
        data = pd.concat([res.set_index("most similar by euclidean").loc[:,"euclidean distance"], y], axis=1)
        sns.regplot(x=data["euclidean distance"], y=data[feat], ci=None, ax=ax[2])
        m, b, r, p, sd = stats.linregress(x=data["euclidean distance"], y=data[feat])
        ax[2].annotate(f"Slope: {m}\nIntercept: {b}\nPCC: {r}\nP-value: {p}\nSE: {sd}",
            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, va='top')
        ax[2].set_xlabel("Euclidean Distance with S288C (one isolate per clade)")
        ax[2].set_ylabel(f"SHAP value of top gene {top_features.loc[top_features['Unnamed: 0']==env, 'cnv'].values[0]}")
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/Section_5/{env}_top_gene_shap_vs_clade_euclidean_distance_to_S288C.pdf")
        plt.close()
