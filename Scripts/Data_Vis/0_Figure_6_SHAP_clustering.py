#!/usr/bin/env python3
"""
Cluster SHAP values and plot the cluster dendrogram, heatmap and clustered
fitness label values. Conduct a linear regression to compare median fitness
between clusters and a Mann-Whitney U test to compare fitness distributions
between clusters.
"""

__author__ = "Kenia Segura Abá"

import os, re
import pandas as pd
import datatable as dt
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import rgb2hex
from scipy.cluster import hierarchy
from scipy.stats import linregress, mannwhitneyu, pearsonr

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/")

def sort_clusters(dend, df_top, fitness, iso_col):
    """
    Sort the dendrogram clusters by median fitness value
    """
    
    # Get the cluster median fitness
    cluster_fitness = fitness.groupby(iso_col["Dend_Cluster"]).median().sort_values()
    
    # Sort the clusters by median fitness
    iso_col.insert(4, "idx", range(iso_col.shape[0])) # add index column
    fitness_sorted = {}
    df_top_sorted = {}
    dend_sorted = {k:[] for k in dend.keys()}
    # dict_keys(['icoord', 'dcoord', 'ivl', 'leaves', 'color_list', 'leaves_color_list'])
    
    for c in cluster_fitness.index:
        fitness_sorted[c] = fitness[iso_col["Dend_Cluster"]==c]
        df_top_sorted[c] = df_top.loc[:,iso_col["Dend_Cluster"]==c]
        
        idx_lst = iso_col[iso_col["Dend_Cluster"]==c].idx
        dend_sorted["ivl"].extend([dend["ivl"][i] for i in idx_lst])
        dend_sorted["leaves"].extend([dend["leaves"][i] for i in idx_lst])
        dend_sorted["leaves_color_list"].extend([dend["leaves_color_list"][i] for i in idx_lst])
        
        idx_lst_trunc = [i for i in range(len(dend["color_list"])) if dend["color_list"][i]==c]
        dend_sorted["icoord"].extend([dend["icoord"][i] for i in idx_lst_trunc])
        dend_sorted["dcoord"].extend([dend["dcoord"][i] for i in idx_lst_trunc])
        dend_sorted["color_list"].extend([dend["color_list"][i] for i in idx_lst_trunc])
    
    fitness_tmp = pd.DataFrame()
    df_top_tmp = pd.DataFrame()
    for c in fitness_sorted.keys():
        fitness_tmp = pd.concat([fitness_tmp, fitness_sorted[c]])
        df_top_tmp = pd.concat([df_top_tmp, df_top_sorted[c]], axis=1)
    
    iso_col_tmp = iso_col.loc[fitness_tmp.index,:]
    
    return(fitness_tmp, iso_col_tmp, df_top_tmp, dend_sorted)


def plotTree(D_dendro, nrow, ncol, idx_ax_dend, save:bool, save_dir, name=""):
    """ 
    Plot a dendrogram
    Modified from: https://stackoverflow.com/questions/36538090/bigger-color-palette-in-matplotlib-for-scipys-dendrogram-python?noredirect=1#comment60681856_36538090
    """
    
    fig,ax = plt.subplots(nrows=nrow, ncols=ncol, figsize=(25, 5*nrow))
    icoord = np.array( D_dendro['icoord'] )
    dcoord = np.array( D_dendro['dcoord'] )
    color_list = np.array( D_dendro['color_list'] )
    x_min, x_max = icoord.min(), icoord.max()
    y_min, y_max = dcoord.min(), dcoord.max()
    
    for xs, ys, color in zip(icoord, dcoord, color_list):
        ax[idx_ax_dend].plot(xs, ys,  color)
    
    ax[idx_ax_dend].set_xlim(x_min-20, x_max + 0.05*abs(x_max))
    ax[idx_ax_dend].set_ylim(y_min, y_max + 0.1*abs(y_max))
    ax[idx_ax_dend].set_title("Dendrogram", fontsize=20)
    ax[idx_ax_dend].set_xlabel("Clusters", fontsize=15)
    ax[idx_ax_dend].set_ylabel("Distance", fontsize=15)
    ax[idx_ax_dend].tick_params(axis="y", labelsize = 10)
    
    if save == True:
        print(f"Plotting Scripts/Data_Vis/Section_6/{name}.png")
        # Save just the portion _inside_ the second axis's boundaries
        extent = ax[idx_ax_dend].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(f"{save_dir}/{name}.pdf", bbox_inches=extent.expanded(1.1, 1.15)) # pad the area by 15% and 10% in x and y directions
    
    return(fig,ax)


# Need to add the masking and cluster coloring to this function.
def cluster_shap_indiv(f, target_envs, pheno, iso_colors, dtype, plot:bool,
                       save_out:bool, save_dir="", dend_thrsh={'YPDCAFEIN40':0.01},
                       cbar_lim={'snp':{'YPDCAFEIN40':[-0.002, 0.002]}}):
    """ 
    Cluster SHAP values and color the dendrogram by clades for the 
    target environments
    """
    
    env = [e for e in f.split("/")[-1].split("_") if e in target_envs][0]
    if env in target_envs:
        linreg_res = {} # store the linear regression results
        mwu_res = {} # store the mann-whitney U test results
        
        df = pd.read_csv(f"{f}", sep="\t", index_col=0).T # SHAP values dataframe
        
        # select the top 20 features based on median absolute shap value across isolates
        top_feat = df.abs().median(axis=1).sort_values(ascending=False)[:20].index
        df_top = df.loc[top_feat,:]
        
        fitness = pheno.loc[df_top.columns, env]
        iso_col = iso_colors.loc[df_top.columns]
        
        if plot==True:
            print(f"Plotting {f}")
            
            # Cluster isolates and reorder all dataframes
            Z = hierarchy.linkage(df_top.T, method="ward", metric="euclidean")
            dend = hierarchy.dendrogram(Z, no_plot=True, color_threshold=dend_thrsh[env])
            print(f"Number of clusters: {len(np.unique(dend['leaves_color_list']))}")
            df_top = df_top.iloc[:,dend["leaves"]]
            fitness = fitness.loc[df_top.columns]
            iso_col = iso_colors.loc[df_top.columns,:]
            
            # Assign each dendrogram cluster a color
            cmap = mpl.colormaps.get_cmap("hsv").resampled(len(np.unique(dend["leaves_color_list"]))) # sample hsv color palette
            color_list = {np.unique(dend["leaves_color_list"])[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
            iso_col["Dend_Cluster"] = dend["leaves_color_list"]
            iso_col["Dend_Cluster_Color"] = iso_col["Dend_Cluster"].map(color_list) # assign colors to isolate
            
            # Sort the dendrogram clusters by median fitness value and plot (don't save to file yet)
            fitness_sorted, iso_col_sorted, df_top_sorted, dend_sorted = sort_clusters(dend, df_top, fitness, iso_col)
            
            if save_out==True:
                print(f"Saving sorted shap values, fitness values, and isolates to files...")
                df_top_sorted.to_csv(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_sorted_v2.tsv", sep="\t")
                iso_col_sorted.to_csv(f"{save_dir}/isolates_clustered_{env}_top20_{dtype}_sorted_v2.tsv", sep="\t")
            
            plotTree(dend, 2, 1, 0, True, save_dir, f"SHAP_clustered_{env}_top20_{dtype}_dendrogram_unsorted") ; plt.close() # save unsorted dendrogram
            fig, ax_dend = plotTree(dend_sorted, 3, 1, 0, False, save_dir) # save the sorted dendrogram
            
            # Fit a line to the median fitness values in each cluster to get the pearson r and p-value
            grouped_data = pd.DataFrame({
                'Dend_Cluster': iso_col_sorted["Dend_Cluster"],
                'Fitness': fitness_sorted[0]})
            grouped_data_med = grouped_data.groupby('Dend_Cluster').median().reset_index() # median fitness values per cluster
            
            m, b, r_value, p_value, std_err = \
            linregress(grouped_data_med.index, grouped_data_med['Fitness']) # calculate R2 and p-value
            linreg_res[env] = {'m':m, 'b':b, 'r':r_value, 'pval':p_value, 'SE':std_err}
            
            # Conduct a mann-whitney U test to compare fitness distributions between clusters
            mannwhitney_res = {}
            clusters = grouped_data['Dend_Cluster'].unique()
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    group1 = grouped_data[grouped_data['Dend_Cluster']==clusters[i]]['Fitness']
                    group2 = grouped_data[grouped_data['Dend_Cluster']==clusters[j]]['Fitness']
                    stat, p_value = mannwhitneyu(group1, group2)
                    mannwhitney_res[(clusters[i], clusters[j])] = {'stat': stat, 'pval': p_value}
            mwu_res[env] = mannwhitney_res
            
            # Add masked values to dataframe to visualize clusters in heatmap
            mask_top = {}
            mask_iso_col = {}
            for i in range(iso_col_sorted.shape[0]):
                mask_top[df_top_sorted.columns[i]] = df_top_sorted.iloc[:,i]
                mask_iso_col[iso_col_sorted.index[i]] = iso_col_sorted.iloc[i,:]
                
                if i==iso_col_sorted.shape[0]-1:
                    break
                elif iso_col_sorted.iloc[i,2]!=iso_col_sorted.iloc[i+1,2]: # cluster divider mask
                    mask_top[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=df_top_sorted.index)
                    mask_iso_col[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series("#ffffff", index=iso_col_sorted.columns)
            
            mask_top = pd.DataFrame(mask_top)
            mask_iso_col = pd.DataFrame(mask_iso_col).T
            
            # Adjust SHAP values based on given limits passed to cbar_lim
            mask_top.clip(lower=cbar_lim[dtype][env][0], upper=cbar_lim[dtype][env][1], inplace=True)
            
            # Plot
            sns.heatmap(mask_top, cmap="RdBu_r", center=0,
                        cbar_kws={"orientation": "vertical", "location":"right"},
                        yticklabels=True, xticklabels=False, ax=ax_dend[1])
            
            # Create an additional axis for the column colors
            ax_dend[1].tick_params(axis="x", which="major", pad=20, length=0) # extra padding for column colors
            for i in range(mask_iso_col.shape[0]):
                ax_dend[1].add_patch(plt.Rectangle(xy=(i,-.4), width=0.016, height=0.5, color=mask_iso_col["Clade_Color"][i]))
            
            sns.violinplot(x=iso_col_sorted["Dend_Cluster"], y=fitness_sorted[0], hue=iso_col_sorted["Dend_Cluster"], palette=color_list, inner="quart", ax=ax_dend[2])
            fig.savefig(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_sorted_v2.pdf")
            plt.close()
            print(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_sorted_v2.pdf")
            
        return df_top, linreg_res, mwu_res


def cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno, iso_colors, dend_thrsh, cbar_lim, save_dir):
    '''Cluster SHAP values based on the SNP model and re-order the PAV and CNV models'''
    linreg_res = {}
    mwu_res = {}
    for env in target_envs:
        print("Clustering SHAP values for", env)
        
        # Get file and read in SHAP values
        f_snp = [f"{f}" for f in snp_files if [f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        f_pav = [f"{f}" for f in pav_files if [f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        f_cnv = [f"{f}" for f in cnv_files if [f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        top_snp, lin_reg_snp, mwu_snp = cluster_shap_indiv(f_snp, target_envs, pheno, iso_colors, dtype="snp", plot=False, save_out=False, save_dir=save_dir, dend_thrsh=dend_thrsh)
        top_pav, lin_reg_pav, mwu_pav = cluster_shap_indiv(f_pav, target_envs, pheno, iso_colors, dtype="pav", plot=False, save_out=False, save_dir=save_dir, dend_thrsh=dend_thrsh)
        top_cnv, lin_reg_cnv, mwu_cnv = cluster_shap_indiv(f_cnv, target_envs, pheno, iso_colors, dtype="cnv", plot=False, save_out=False, save_dir=save_dir, dend_thrsh=dend_thrsh)
        linreg_res[env] = {'snp':lin_reg_snp, 'pav':lin_reg_pav, 'cnv':lin_reg_cnv}
        mwu_res[env] = {'snp':mwu_snp, 'pav':mwu_pav, 'cnv':mwu_cnv}
        
        # Cluster isolates by SNPs
        Z = hierarchy.linkage(top_snp.T, method="ward", metric="euclidean")
        dend_snp = hierarchy.dendrogram(Z, no_plot=True, color_threshold=dend_thrsh[env])
        top_snp = top_snp.iloc[:,dend_snp["leaves"]]
        fitness = pheno.loc[top_snp.columns, env]
        iso_col = iso_colors.loc[top_snp.columns]
        
        # Assign each dendrogram cluster a color
        cmap = mpl.colormaps.get_cmap("hsv").resampled(len(np.unique(dend_snp["leaves_color_list"]))) # sample hsv color palette
        color_list = {np.unique(dend_snp["leaves_color_list"])[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
        iso_col["Dend_Cluster"] = dend_snp["leaves_color_list"]
        iso_col["Dend_Cluster_Color"] = iso_col["Dend_Cluster"].map(color_list) # assign colors to isolate
        
        # Sort the dendrogram clusters by median fitness value and reorder all dataframes
        fitness_sorted, iso_col_sorted, top_snp_sorted, dend_snp_sorted = sort_clusters(dend_snp, top_snp, fitness, iso_col)
        plotTree(dend_snp, 2, 1, 0, True, save_dir, f"SHAP_clustered_by_snp_{env}_top20_dendrogram_v2") ; plt.close()
        top_pav = top_pav.iloc[:,dend_snp_sorted["leaves"]]
        top_cnv = top_cnv.iloc[:,dend_snp_sorted["leaves"]]
        
        # Add masked values to dataframe to visualize clusters in heatmap
        mask_top_snp = {}
        mask_top_pav = {}
        mask_top_cnv = {}
        mask_iso_col = {}
        for i in range(iso_col_sorted.shape[0]):
            mask_top_snp[top_snp_sorted.columns[i]] = top_snp_sorted.iloc[:,i]
            mask_top_pav[top_pav.columns[i]] = top_pav.iloc[:,i]
            mask_top_cnv[top_cnv.columns[i]] = top_cnv.iloc[:,i]
            mask_iso_col[iso_col_sorted.index[i]] = iso_col_sorted.iloc[i,:]
            
            if i==624:
                break
            elif iso_col_sorted.iloc[i,2]!=iso_col_sorted.iloc[i+1,2]: # cluster division
                mask_top_snp[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=top_snp_sorted.index)
                mask_top_cnv[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=top_cnv.index)
                mask_top_pav[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=top_pav.index)
                mask_iso_col[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series("#ffffff", index=iso_col_sorted.columns)
        
        mask_top_snp = pd.DataFrame(mask_top_snp)
        mask_top_pav = pd.DataFrame(mask_top_pav)
        mask_top_cnv = pd.DataFrame(mask_top_cnv)
        mask_iso_col = pd.DataFrame(mask_iso_col).T
        
        # Adjust SHAP values based on given limits passed to cbar_lim
        mask_top_snp.clip(lower=cbar_lim["snp"][env][0], upper=cbar_lim["snp"][env][1], inplace=True)
        mask_top_pav.clip(lower=cbar_lim["pav"][env][0], upper=cbar_lim["pav"][env][1], inplace=True)
        mask_top_cnv.clip(lower=cbar_lim["cnv"][env][0], upper=cbar_lim["cnv"][env][1], inplace=True)
        
        # Plot
        fig, ax = plt.subplots(4,1, figsize=(13,8.5))
        s = sns.heatmap(mask_top_snp, cmap="RdBu_r", center=0,
                        cbar_kws={"orientation": "vertical", "location":"right"},
                        yticklabels=True, xticklabels=False, ax=ax[0])
        
        # Create an additional axis for the column colors
        ax[0].tick_params(axis="x", which="major", pad=20, length=0) # extra padding for column colors
        for i in range(iso_col_sorted.shape[0]):
            ax[0].add_patch(plt.Rectangle(xy=(i,-.4), width=0.0032, height=0.5, color=mask_iso_col["Clade_Color"][i]))
            # ax[0].add_patch(plt.Rectangle(xy=(i,-.2), width=0.0064, height=0.6, color=mask_iso_col["Dend_Cluster_Color"][i]))
        
        # Plot the other two data types
        p = sns.heatmap(mask_top_pav, cmap="RdBu_r", center=0,
                        cbar_kws={"orientation": "vertical", "location":"right"},
                        yticklabels=True, xticklabels=False, ax=ax[1])
        c = sns.heatmap(mask_top_cnv, cmap="RdBu_r", center=0,
                        cbar_kws={"orientation": "vertical", "location":"right"},
                        yticklabels=True, xticklabels=False, ax=ax[2])
        
        # Group fitness by dendrogram cluster and plot violin
        v = sns.violinplot(x=iso_col_sorted["Dend_Cluster"],
                           y=fitness_sorted[0],
                           hue=iso_col_sorted["Dend_Cluster"],
                           palette=color_list, inner="quart", ax=ax[3])
        plt.tight_layout()
        plt.savefig(f"{save_dir}/SHAP_clustered_by_snp_{env}_top20_sorted_v2.pdf")
        plt.close()
        print(f"Saved as {save_dir}/SHAP_clustered_by_snp_{env}_top20_sorted_v2.pdf")
        
        del top_snp, top_pav, top_cnv, fitness, iso_col, mask_iso_col
        
    return linreg_res, mwu_res


def save_results(linreg_res, mwu_res, linreg_save, mwu_save):
    # Convert linear regression results to a dataframe
    linreg_df = pd.DataFrame.from_dict({(i,j,k): linreg_res[i][j][k]
                                         for i in linreg_res.keys()
                                         for j in linreg_res[i].keys()
                                         for k in linreg_res[i][j].keys()},
                                        orient='index')
    linreg_df = linreg_df.droplevel(0)
    linreg_df.sort_index(inplace=True)
    linreg_df.to_csv(linreg_save)
    
    # Convert Mann-Whitney U test results to a dataframe
    mwu_df = pd.DataFrame.from_dict({(i,j,k,l): mwu_res[i][j][k][l]
                                      for i in mwu_res.keys()
                                      for j in mwu_res[i].keys()
                                      for k in mwu_res[i][j].keys()
                                      for l in mwu_res[i][j][k].keys()},
                                     orient='index')
    mwu_df = mwu_df.droplevel(0)
    mwu_df.sort_index(inplace=True)
    mwu_df.reset_index(level=[2], inplace=True)
    mwu_df.insert(0, "Cluster1", [i[0] for i in mwu_df.level_2])
    mwu_df.insert(1, "Cluster2", [i[1] for i in mwu_df.level_2])
    mwu_df.drop(columns="level_2").to_csv(mwu_save)


if __name__ == "__main__":
    # Clade & fitness information
    clades = pd.read_excel("Data/Peter_2018/0_raw_data/Peter_2018_Supplementary_Tables.xls", sheet_name="Table S1", skiprows=3, nrows=1011) # isolate clades
    pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # isolate fitness

    # Map clades to isolates and create a colormap
    clades = clades[["Standardized name", "Clades"]] # subset relevant columns
    clades = clades.loc[clades["Standardized name"].isin(pheno.index)] # diploid isolates
    clades.set_index("Standardized name", inplace=True)
    clades.loc[clades.Clades.isna(),"Clades"] = "Unknown" # replace NaN with Unknown

    cmap = mpl.colormaps.get_cmap("hsv").resampled(clades.Clades.nunique()) # sample hsv color palette
    color_list = {clades.Clades.unique()[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
    iso_colors = clades['Clades'].map(color_list) # assign colors to isolates
    iso_colors = pd.concat([clades["Clades"], iso_colors], axis=1)
    iso_colors.columns = ["Clades", "Clade_Color"]

    map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t") # map of SNPs to genes
    map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t") # map of ORFs to genes
    
    ############### Read in SHAP value files for baseline models ###############
    # Note: plots saved to /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_6/baseline_top20/
    
    snp_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/SNP/baseline"
    files = os.listdir(snp_path)
    snp_files = [f"{snp_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    pav_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/PAV/baseline"
    files = os.listdir(pav_path)
    pav_files = [f"{pav_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    cnv_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/CNV/baseline"
    files = os.listdir(cnv_path)
    cnv_files = [f"{cnv_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]

    target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    snp_files = [f for f in snp_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    pav_files = [f for f in pav_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    cnv_files = [f for f in cnv_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    
    # Cluster SHAP values
    dend_thrsh_snp={"YPDCAFEIN40":0.002, "YPDCAFEIN50":0.002, "YPDSODIUMMETAARSENITE":0.01, "YPDBENOMYL500":0.002, "YPDCUSO410MM":0.04}
    dend_thrsh_pav={"YPDCAFEIN40":0.008, "YPDCAFEIN50":0.008, "YPDSODIUMMETAARSENITE":0.01, "YPDBENOMYL500":0.01, "YPDCUSO410MM":0.04}
    dend_thrsh_cnv={"YPDCAFEIN40":0.01, "YPDCAFEIN50":0.01, "YPDSODIUMMETAARSENITE":0.01, "YPDBENOMYL500":0.01, "YPDCUSO410MM":0.02}
    cbar_lim={"snp":{"YPDCAFEIN40":[-0.002, 0.002], "YPDCAFEIN50":[-0.002, 0.002], "YPDSODIUMMETAARSENITE":[-0.02, 0.02], "YPDBENOMYL500":[-0.002,0.002], "YPDCUSO410MM":[-0.02, 0.02]}, # SHAP value limits for colorbar
              "pav":{"YPDCAFEIN40":[-0.005, 0.005], "YPDCAFEIN50":[-0.005, 0.005], "YPDSODIUMMETAARSENITE":[-0.01, 0.01], "YPDBENOMYL500":[-0.005, 0.005], "YPDCUSO410MM":[-0.03, 0.03]},
              "cnv":{"YPDCAFEIN40":[-0.005, 0.005], "YPDCAFEIN50":[-0.005, 0.005], "YPDSODIUMMETAARSENITE":[-0.01, 0.01], "YPDBENOMYL500":[-0.007, 0.007], "YPDCUSO410MM":[-0.005, 0.005]}}
    
    cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno,
                          iso_colors, dend_thrsh_snp, cbar_lim=cbar_lim,
                          save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
    
    linreg_res = {}
    mwu_res = {}
    for i in range(len(snp_files)):
        _, snp_lin, snp_mwu = cluster_shap_indiv(f"{snp_files[i]}",
            target_envs, pheno, iso_colors, "snp", plot=True, save_out=False,
            dend_thrsh=dend_thrsh_snp, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
        _, pav_lin, pav_mwu = cluster_shap_indiv(f"{pav_files[i]}",
           target_envs, pheno, iso_colors, "pav", plot=True, save_out=False,
            dend_thrsh=dend_thrsh_pav, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
        _, cnv_lin, cnv_mwu = cluster_shap_indiv(f"{cnv_files[i]}",
            target_envs, pheno, iso_colors, "cnv", plot=True, save_out=False,
            dend_thrsh=dend_thrsh_cnv, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
        linreg_res[i] = {"snp":snp_lin, "pav":pav_lin, "cnv":cnv_lin}
        mwu_res[i] = {"snp":snp_mwu, "pav":pav_mwu, "cnv":cnv_mwu}

    linreg_save = "Scripts/Data_Vis/Section_6/baseline_top20/SHAP_cluster_median_fitness_linreg_RF_baseline_models.csv"
    mwu_save = "Scripts/Data_Vis/Section_6/baseline_top20/SHAP_cluster_fitness_mwu_RF_baseline_models.csv"
    save_results(linreg_res, mwu_res, linreg_save, mwu_save)

    ############## Read in SHAP value files for optimized models ###############
    # Note: plots saved to /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_6/fs_top20/
    
    snp_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/SNP/fs"
    files = os.listdir(snp_path)
    snp_files = [f"{snp_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    pav_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/PAV/fs"
    files = os.listdir(pav_path)
    pav_files = [f"{pav_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    cnv_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/CNV/fs"
    files = os.listdir(cnv_path)
    cnv_files = [f"{cnv_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]

    # Cluster SHAP values
    target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    snp_files = [f for f in snp_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    pav_files = [f for f in pav_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    cnv_files = [f for f in cnv_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]

    dend_thrsh_snp={"YPDCAFEIN40":0.008, "YPDCAFEIN50":0.02, "YPDSODIUMMETAARSENITE":0.025, "YPDBENOMYL500":0.01, "YPDCUSO410MM":0.05}
    dend_thrsh_pav={"YPDCAFEIN40":0.04, "YPDCAFEIN50":0.03, "YPDSODIUMMETAARSENITE":0.05, "YPDBENOMYL500":0.05, "YPDCUSO410MM":0.1}
    dend_thrsh_cnv={"YPDCAFEIN40":0.05, "YPDCAFEIN50":0.03, "YPDSODIUMMETAARSENITE":0.125, "YPDBENOMYL500":0.05, "YPDCUSO410MM":0.5}
    cbar_lim={"snp":{"YPDCAFEIN40":[-0.008, 0.008], "YPDCAFEIN50":[-0.02, 0.02], "YPDSODIUMMETAARSENITE":[-0.025, 0.025], "YPDBENOMYL500":[-0.01,0.01], "YPDCUSO410MM":[-0.05, 0.05]}, # SHAP value limits for colorbar
                "pav":{"YPDCAFEIN40":[-0.04, 0.04], "YPDCAFEIN50":[-0.03, 0.03], "YPDSODIUMMETAARSENITE":[-0.05, 0.05], "YPDBENOMYL500":[-0.05, 0.05], "YPDCUSO410MM":[-0.1, 0.1]},
                "cnv":{"YPDCAFEIN40":[-0.05, 0.05], "YPDCAFEIN50":[-0.03, 0.03], "YPDSODIUMMETAARSENITE":[-0.125, 0.125], "YPDBENOMYL500":[-0.05, 0.05], "YPDCUSO410MM":[-0.5, 0.5]}}

    cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno,
                            iso_colors, dend_thrsh_snp, cbar_lim=cbar_lim,
                            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20")

    linreg_res = {}
    mwu_res = {}
    for i in range(len(snp_files)):
        _, snp_lin, snp_mwu = cluster_shap_indiv(f"{snp_files[i]}",
            target_envs, pheno, iso_colors, "snp", plot=True, save_out=True,
            dend_thrsh=dend_thrsh_snp, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20")
        _, pav_lin, pav_mwu = cluster_shap_indiv(f"{pav_files[i]}",
            target_envs, pheno, iso_colors, "pav", plot=True, save_out=True,
            dend_thrsh=dend_thrsh_pav, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20")
        _, cnv_lin, cnv_mwu = cluster_shap_indiv(f"{cnv_files[i]}",
            target_envs, pheno, iso_colors, "cnv", plot=True, save_out=True,
            dend_thrsh=dend_thrsh_cnv, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20")
        linreg_res[i] = {"snp":snp_lin, "pav":pav_lin, "cnv":cnv_lin}
        mwu_res[i] = {"snp":snp_mwu, "pav":pav_mwu, "cnv":cnv_mwu}


    linreg_save = "Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_cluster_median_fitness_linreg_RF_FS_models.csv"
    mwu_save = "Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_cluster_fitness_mwu_RF_FS_models.csv"
    save_results(linreg_res, mwu_res, linreg_save, mwu_save)

    # Plot SHAP values for the top 20 features in each SNP, PAV, CNV model
    ## label which ones are the benchmark genes
    ## Fitness on x axis
    ## feature value on y axis
    ## color by SHAP value
    target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # isolate fitness
    pav = dt.fread("Data/Peter_2018/ORFs_pres_abs.csv").to_pandas() # ORF presence/absence
    cnv = dt.fread("Data/Peter_2018/ORFs_no_NA.csv").to_pandas() # CNV data
    snp = dt.fread("Data/Peter_2018/geno.csv").to_pandas() # SNP data
    pav.set_index("ID", inplace=True)
    cnv.set_index("ID", inplace=True)
    snp.set_index("ID", inplace=True)
    pav.columns = pav.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix ORF IDs to match map_orfs
    pav.columns = pav.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
    cnv.columns = cnv.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix ORF IDs to match map_orfs
    cnv.columns = cnv.apply(lambda x: re.sub("\.", "-", x.name), axis=0)

    res = []
    for env in target_envs:
        for data_type in ["snp", "pav", "cnv"]:
            df = pd.read_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v2.tsv", sep="\t", index_col=0)
            
            if data_type == "snp": # Add benchmark gene information to the dataframe
                df_info = pd.DataFrame(df.index).merge(
                    map_snps.set_index("snp").drop(columns=["chr", "pos"]),
                    left_on=0, right_index=True, how="left")
            else:
                df.index = df.apply(lambda x: re.sub("^X", "", x.name), axis=1) # fix ORF IDs to match map_orfs
                df.index = df.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
                df_info = pd.DataFrame(df.index).merge(
                    map_orfs.set_index("orf").drop(columns="organism"),
                    left_on=0, right_index=True, how="left")
            
            # plots for individual features
            fig, ax = plt.subplots(10, 2, figsize=(9, 26))
            for i in range(20):
                # Center the SHAP values legend at 0
                vmin = df.iloc[i, :].min()
                vmax = df.iloc[i, :].max()
                norm = plt.Normalize(vmin=vmin, vmax=vmax)
                
                if data_type == "pav":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=pav.loc[df.columns, df.index[i]].astype(int),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor=None, alpha=0.7, legend=False)
                    ax[i//2][i%2].set_ylabel("PAV value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=pav.loc[df.columns, df.index[i]].astype(int))
                    res.append([env, data_type, df.index[i], r, p])
                if data_type == "cnv":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=cnv.loc[df.columns, df.index[i]].astype(float),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor=None, alpha=0.7, legend=False)
                    ax[i//2][i%2].set_ylabel("CNV value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=cnv.loc[df.columns, df.index[i]].astype(float))
                    res.append([env, data_type, df.index[i], r, p])
                if data_type == "snp":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=snp.loc[df.columns, df.index[i]].astype(int),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor=None, alpha=0.7, legend=False)
                    ax[i//2][i%2].set_ylabel("SNP value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=snp.loc[df.columns, df.index[i]].astype(int))
                    res.append([env, data_type, df.index[i], r, p])
                
                ax[i//2][i%2].set_title(f"Feature {'//'.join(df_info.iloc[i].astype(str))}",
                                        fontsize=8)
                ax[i//2][i%2].set_xlabel("Fitness")
                sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
                sm.set_array([])
                plt.colorbar(sm, ax=ax[i//2][i%2], orientation="vertical", label="SHAP value")
            
            for axis in ax.flat:
                axis.set_box_aspect(1)
            
            plt.tight_layout()
            plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v2_indiv_feat_trends.pdf",
                        bbox_inches="tight", dpi=300)
            plt.close()
            
            # plots for pairs of features
            fig, ax = plt.subplots(20, 20, figsize=(80, 80))
            for i in range(20):
                for j in range(20):
                    if (i != j) & (i < j):
                        vmax = pheno.loc[df.columns, env].max()
                        norm = plt.Normalize(vmin=0, vmax=vmax)
                        
                        if data_type == "pav":
                            size_column = pav.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="viridis", edgecolor=None, alpha=0.7, hue_norm=norm,
                                size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]}")
                            ax[i][j].set_xlabel(f"{df.index[i]}")
                        if data_type == "cnv":
                            size_column = cnv.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="viridis", edgecolor=None, alpha=0.7, hue_norm=norm,
                                size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]}")
                            ax[i][j].set_xlabel(f"{df.index[i]}")
                        if data_type == "snp":
                            size_column = snp.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="viridis", edgecolor=None, alpha=0.7, hue_norm=norm,
                                size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]} SHAP")
                            ax[i][j].set_xlabel(f"{df.index[i]} SHAP")
                        
                        ax[i][j].set_title(f"{df.index[i]} vs {df.index[j]}",
                                            fontsize=8)
                        ax[i][j].axvline(x=0, color="black", linestyle="--")
                        ax[i][j].axhline(y=0, color="black", linestyle="--")
                        
                        # Create separate legends for hue y size
                        handles, labels = scatter.get_legend_handles_labels()
                        size_legend = ax[i][j].legend(
                            handles[-len(size_column.unique()):],
                            labels[-len(size_column.unique()):],
                            bbox_to_anchor=(1.05, 0.5), loc='center left',
                            title="Feature value", prop={'size': 6})
                        hue_legend = ax[i][j].legend(
                            handles[1:-len(size_column.unique())-1],
                            labels[1:-len(size_column.unique())-1],
                            bbox_to_anchor=(1.05, 1), loc='upper left',
                            title="Fitness", prop={'size': 6})
                        ax[i][j].add_artist(size_legend)
            
            for axis in ax.flat:
                axis.set_box_aspect(1)
            
            plt.tight_layout()
            plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v2_feat_pair_trends.pdf",
                        bbox_inches="tight", dpi=300)
            plt.close()


    res_df = pd.DataFrame(res, columns=["Environment", "Data_Type", "Feature", "R", "P-value"])
    res_df.to_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_FS_top20_indiv_feat_trends_pearsonr.tsv")
    
    ######## Read in SHAP value files for max_gini_from_RF_baseline optimized FS models ########
    path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction"
    files = os.listdir(f"{path}/SNP/best_rf_fs")
    snp_files = [f"{path}/SNP/best_rf_fs/{f}" for f in files if re.search(r"SHAP_values_sorted_snp_[A-Za-z0-9_]+_[RFa-z0-9-_]+_training\.txt", f)]
    files = os.listdir(f"{path}/PAV/best_rf_fs")
    pav_files = [f"{path}/PAV/best_rf_fs/{f}" for f in files if re.search(r"SHAP_values_sorted_pav_[A-Za-z0-9_]+_[RFa-z0-9-_]+_training\.txt", f)]
    files = os.listdir(f"{path}/CNV/best_rf_fs")
    cnv_files = [f"{path}/CNV/best_rf_fs/{f}" for f in files if re.search(r"SHAP_values_sorted_cnv_[A-Za-z0-9_]+_[RFa-z0-9-_]+_training\.txt", f)]

    # Cluster SHAP values
    target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    dend_thrsh_snp={"YPDCAFEIN40":0.03, "YPDCAFEIN50":0.025, "YPDSODIUMMETAARSENITE":0.04, "YPDBENOMYL500":0.03, "YPDCUSO410MM":0.2}
    dend_thrsh_pav={"YPDCAFEIN40":0.045, "YPDCAFEIN50":0.04, "YPDSODIUMMETAARSENITE":0.15, "YPDBENOMYL500":0.06, "YPDCUSO410MM":0.2}
    dend_thrsh_cnv={"YPDCAFEIN40":0.045, "YPDCAFEIN50":0.05, "YPDSODIUMMETAARSENITE":0.3, "YPDBENOMYL500":0.04, "YPDCUSO410MM":1}
    cbar_lim={"snp":{"YPDCAFEIN40":[-0.03, 0.03], "YPDCAFEIN50":[-0.025, 0.025], "YPDSODIUMMETAARSENITE":[-0.04, 0.04], "YPDBENOMYL500":[-0.03, 0.03], "YPDCUSO410MM":[-0.2, 0.2]}, # SHAP value limits for colorbar
              "pav":{"YPDCAFEIN40":[-0.045, 0.045], "YPDCAFEIN50":[-0.05, 0.05], "YPDSODIUMMETAARSENITE":[-0.15, 0.15], "YPDBENOMYL500":[-0.06, 0.06], "YPDCUSO410MM":[-0.2, 0.2]},
              "cnv":{"YPDCAFEIN40":[-0.045, 0.045], "YPDCAFEIN50":[-0.05, 0.05], "YPDSODIUMMETAARSENITE":[-0.3, 0.3], "YPDBENOMYL500":[-0.04, 0.04], "YPDCUSO410MM":[-1, 1]}}
    
    cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno,
                          iso_colors, dend_thrsh_snp, cbar_lim=cbar_lim,
                          save_dir="Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20") # cluster PAV and CNV by SNPs
    
    linreg_res = {}
    mwu_res = {}
    for i in range(len(snp_files)):
        _, snp_lin, snp_mwu = cluster_shap_indiv(f"{snp_files[i]}",
            target_envs, pheno, iso_colors, "snp", plot=True, save_out=True,
            dend_thrsh=dend_thrsh_snp, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20")
        _, pav_lin, pav_mwu = cluster_shap_indiv(f"{pav_files[i]}",
           target_envs, pheno, iso_colors, "pav", plot=True, save_out=True,
            dend_thrsh=dend_thrsh_pav, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20")
        _, cnv_lin, cnv_mwu = cluster_shap_indiv(f"{cnv_files[i]}",
            target_envs, pheno, iso_colors, "cnv", plot=True, save_out=True,
            dend_thrsh=dend_thrsh_cnv, cbar_lim=cbar_lim,
            save_dir="Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20")
        linreg_res[i] = {"snp":snp_lin, "pav":pav_lin, "cnv":cnv_lin}
        mwu_res[i] = {"snp":snp_mwu, "pav":pav_mwu, "cnv":cnv_mwu}
    
    linreg_save = "Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20/SHAP_cluster_median_fitness_linreg_shap_int_best_rf_fs_models.csv"
    mwu_save = "Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20/SHAP_cluster_fitness_mwu_shap_int_best_rf_fs_models.csv"
    save_results(linreg_res, mwu_res, linreg_save, mwu_save)
    
    # Plot SHAP values for the top 20 features in each SNP, PAV, CNV model
    ## label which ones are the benchmark genes
    ## Fitness on x axis
    ## feature value on y axis
    ## color by SHAP value
    target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0) # isolate fitness
    pav = dt.fread("Data/Peter_2018/ORFs_pres_abs.csv").to_pandas() # ORF presence/absence
    cnv = dt.fread("Data/Peter_2018/ORFs_no_NA.csv").to_pandas() # CNV data
    snp = dt.fread("Data/Peter_2018/geno.csv").to_pandas() # SNP data
    pav.set_index("ID", inplace=True)
    cnv.set_index("ID", inplace=True)
    snp.set_index("ID", inplace=True)
    pav.columns = pav.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix ORF IDs to match map_orfs
    pav.columns = pav.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
    cnv.columns = cnv.apply(lambda x: re.sub("^X", "", x.name), axis=0) # fix ORF IDs to match map_orfs
    cnv.columns = cnv.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
    
    res = []
    for env in target_envs:
        for data_type in ["snp", "pav", "cnv"]:
            df = pd.read_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v2.tsv", sep="\t", index_col=0)
            
            if data_type == "snp": # Add benchmark gene information to the dataframe
                df_info = pd.DataFrame(df.index).merge(
                    map_snps.set_index("snp").drop(columns=["chr", "pos"]),
                    left_on=0, right_index=True, how="left")
            else:
                df.index = df.apply(lambda x: re.sub("^X", "", x.name), axis=1) # fix ORF IDs to match map_orfs
                df.index = df.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
                df_info = pd.DataFrame(df.index).merge(
                    map_orfs.set_index("orf").drop(columns="organism"),
                    left_on=0, right_index=True, how="left")
            
            # plots for individual features
            fig, ax = plt.subplots(10, 2, figsize=(9, 26))
            for i in range(20):
                # Center the SHAP values legend at 0
                vmin = df.iloc[i, :].min()
                vmax = df.iloc[i, :].max()
                norm = plt.Normalize(vmin=vmin, vmax=vmax)
                
                if data_type == "pav":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=pav.loc[df.columns, df.index[i]].astype(int),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor=None, alpha=0.7, legend=False)
                    ax[i//2][i%2].set_ylabel("PAV value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=pav.loc[df.columns, df.index[i]].astype(int))
                    res.append([env, data_type, df.index[i], r, p])
                if data_type == "cnv":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=cnv.loc[df.columns, df.index[i]].astype(float),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor=None, alpha=0.7, legend=False)
                    ax[i//2][i%2].set_ylabel("CNV value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=cnv.loc[df.columns, df.index[i]].astype(float))
                    res.append([env, data_type, df.index[i], r, p])
                if data_type == "snp":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=snp.loc[df.columns, df.index[i]].astype(int),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor=None, alpha=0.7, legend=False)
                    ax[i//2][i%2].set_ylabel("SNP value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=snp.loc[df.columns, df.index[i]].astype(int))
                    res.append([env, data_type, df.index[i], r, p])
                
                ax[i//2][i%2].set_title(f"Feature {'//'.join(df_info.iloc[i].astype(str))}",
                                        fontsize=8)
                ax[i//2][i%2].set_xlabel("Fitness")
                sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
                sm.set_array([])
                plt.colorbar(sm, ax=ax[i//2][i%2], orientation="vertical", label="SHAP value")
            
            for axis in ax.flat:
                axis.set_box_aspect(1)
            
            plt.tight_layout()
            plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v2_indiv_feat_trends.pdf",
                        bbox_inches="tight", dpi=300)
            plt.close()
            
            # plots for pairs of features
            fig, ax = plt.subplots(20, 20, figsize=(80, 80))
            for i in range(20):
                for j in range(20):
                    if (i != j) & (i < j):
                        vmax = pheno.loc[df.columns, env].max()
                        norm = plt.Normalize(vmin=0, vmax=vmax)
                        # cmap = sns.diverging_palette(240, 10, n=9, as_cmap=True)
                        
                        if data_type == "pav":
                            size_column = pav.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="viridis", edgecolor=None, alpha=0.7, hue_norm=norm,
                                size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]}")
                            ax[i][j].set_xlabel(f"{df.index[i]}")
                        if data_type == "cnv":
                            size_column = cnv.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="viridis", edgecolor=None, alpha=0.7, hue_norm=norm,
                                size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]}")
                            ax[i][j].set_xlabel(f"{df.index[i]}")
                        if data_type == "snp":
                            size_column = snp.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="viridis", edgecolor=None, alpha=0.7, hue_norm=norm,
                                size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]} SHAP")
                            ax[i][j].set_xlabel(f"{df.index[i]} SHAP")
                        
                        ax[i][j].set_title(f"{df.index[i]} vs {df.index[j]}",
                                            fontsize=8)
                        ax[i][j].axvline(x=0, color="black", linestyle="--")
                        ax[i][j].axhline(y=0, color="black", linestyle="--")
                        
                        # Create separate legends for hue y size
                        handles, labels = scatter.get_legend_handles_labels()
                        size_legend = ax[i][j].legend(
                            handles[-len(size_column.unique()):],
                            labels[-len(size_column.unique()):],
                            bbox_to_anchor=(1.05, 0.5), loc='center left',
                            title="Feature value", prop={'size': 6})
                        hue_legend = ax[i][j].legend(
                            handles[1:-len(size_column.unique())-1],
                            labels[1:-len(size_column.unique())-1],
                            bbox_to_anchor=(1.05, 1), loc='upper left',
                            title="Fitness", prop={'size': 6})
                        ax[i][j].add_artist(size_legend)
            
            for axis in ax.flat:
                axis.set_box_aspect(1)
            
            plt.tight_layout()
            plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v2_feat_pair_trends.pdf",
                        bbox_inches="tight", dpi=300)
            plt.close()


    res_df = pd.DataFrame(res, columns=["Environment", "Data_Type", "Feature", "R", "P-value"])
    res_df.to_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/shap_int_best_rf_fs_top20/SHAP_clustered_FS_top20_indiv_feat_trends_pearsonr.tsv")

