"""
Cluster SHAP values and plot the cluster dendrogram, heatmap and clustered fitness label values
"""
__author__ = "Kenia Segura Abá"
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import rgb2hex
from scipy.cluster import hierarchy

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


def plotTree(D_dendro, nrow, ncol, idx_ax_dend, save:bool, name=""):
    """ 
    Plot a dendrogram
    Modified from: https://stackoverflow.com/questions/36538090/bigger-color-palette-in-matplotlib-for-scipys-dendrogram-python?noredirect=1#comment60681856_36538090
    """
    print(f"Plotting Scripts/Data_Vis/{save}.png")
    fig,ax = plt.subplots(nrows=nrow, ncols=ncol, figsize=(25, 5*nrow))
    icoord = np.array( D_dendro['icoord'] )
    dcoord = np.array( D_dendro['dcoord'] )
    color_list = np.array( D_dendro['color_list'] )
    x_min, x_max = icoord.min(), icoord.max()
    y_min, y_max = dcoord.min(), dcoord.max()
    for xs, ys, color in zip(icoord, dcoord, color_list):
        ax[idx_ax_dend].plot(xs, ys,  color)
    ax[idx_ax_dend].set_xlim( x_min-10, x_max + 0.1*abs(x_max) )
    ax[idx_ax_dend].set_ylim( y_min, y_max + 0.1*abs(y_max) )
    ax[idx_ax_dend].set_title("Dendrogram", fontsize=20)
    ax[idx_ax_dend].set_xlabel("Clusters", fontsize=15)
    ax[idx_ax_dend].set_ylabel("Distance", fontsize=15)
    ax[idx_ax_dend].tick_params(axis="y", labelsize = 10)
    if save==True:
        # Save just the portion _inside_ the second axis's boundaries
        extent = ax[idx_ax_dend].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(f"Scripts/Data_Vis/{name}.pdf", bbox_inches=extent.expanded(1.1, 1.15)) # pad the area by 15% and 10% in x and y directions
    return(fig,ax)


# Need to add the masking and cluster coloring to this function.
def cluster_shap_indiv(f, target_envs, pheno, iso_colors, dtype, plot:bool,
                       dend_thrsh={'YPDCAFEIN40':0.01}, cbar_lim={'snp':{'YPDCAFEIN40':[-0.002, 0.002]}}):
    """ 
    Cluster SHAP values and color the dendrogram by clades for the 
    target environments
    """
    env = f.split("/")[-1].split("_")[3]
    if env in target_envs:
        df = pd.read_csv(f"{f}", sep="\t", index_col=0).T
        # select the top 10 features based on median shap value across isolates
        top_feat = df.median(axis=1).sort_values(ascending=False)[:10].index
        df_top = df.loc[top_feat,:]
        fitness = pheno.loc[df_top.columns, env]
        iso_col = iso_colors.loc[df_top.columns]
        if plot==True:
            print(f"Plotting {f}")
            # Cluster isolates and reorder all dataframes
            Z = hierarchy.linkage(df_top.T, method="ward", metric="euclidean")
            dend = hierarchy.dendrogram(Z, no_plot=True, color_threshold=dend_thrsh[env])
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
            fig, ax_dend = plotTree(dend_sorted, 3, 1, 0, False)
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
            fig.savefig(f"Scripts/Data_Vis/SHAP_clustered_{env}_top10_{dtype}_sorted.pdf")
            plt.close()
            print(f"Saved as Scripts/Data_Vis/SHAP_clustered_{env}_top10_{dtype}_sorted.pdf")   
        return df_top


def cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno, iso_colors, dend_thrsh, cbar_lim):
    for env in target_envs:
        print("Clustering SHAP values for", env)
        # Get file and read inS
        f_snp = [f"Scripts/Data_Vis/SNP_Figures/SHAP/{f}" for f in snp_files if f.split("_")[3] == env][0]
        f_pav = [f"Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_presence_absence/{f}" for f in pav_files if f.split("_")[3] == env][0]
        f_cnv = [f"Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_copy_number/{f}" for f in cnv_files if f.split("_")[3] == env][0]
        top_snp = cluster_shap_indiv(f_snp, target_envs, pheno, iso_colors, dtype="snp", plot=False, dend_thrsh=dend_thrsh)
        top_pav = cluster_shap_indiv(f_pav, target_envs, pheno, iso_colors, dtype="pav", plot=False, dend_thrsh=dend_thrsh)
        top_cnv = cluster_shap_indiv(f_cnv, target_envs, pheno, iso_colors, dtype="cnv", plot=False, dend_thrsh=dend_thrsh)
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
        plotTree(dend_snp, 2, 1, 0, True, f"SHAP_clustered_by_snp_{env}_top10_dendrogram") ; plt.close()
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
        sns.violinplot(x=iso_col_sorted["Dend_Cluster"], y=fitness_sorted[0], hue=iso_col_sorted["Dend_Cluster"], palette=color_list, inner="quart", ax=ax[3])
        plt.tight_layout()
        plt.savefig(f"Scripts/Data_Vis/SHAP_clustered_by_snp_{env}_top10_sorted.pdf")
        plt.close()
        print(f"Saved as Scripts/Data_Vis/SHAP_clustered_by_snp_{env}_top10_sorted.pdf")
        del top_snp, top_pav, top_cnv, fitness, iso_col, mask_iso_col


if __name__ == "__main__":
    # Read in SHAP value files
    files = os.listdir("Scripts/Data_Vis/SNP_Figures/SHAP/")
    snp_files = [f for f in files if f.startswith("SHAP_values_sorted_Y")]
    files = os.listdir("Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_presence_absence/")
    pav_files = [f for f in files if f.startswith("SHAP_values_sorted_Y")]
    files = os.listdir("Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_copy_number/")
    cnv_files = [f for f in files if f.startswith("SHAP_values_sorted_Y")]

    # Clade & fitness information
    clades = pd.read_excel("Data/Peter_2018/Peter_2018_Supplementary_Tables.xls", sheet_name="Table S1", skiprows=3, nrows=1011) # isolate clades
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

    # Cluster SHAP values
    target_envs = ["YPDCAFEIN40", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    dend_thrsh_snp={"YPDCAFEIN40":0.007, "YPDSODIUMMETAARSENITE":0.02, "YPDBENOMYL500":0.008, "YPDCUSO410MM":0.07}
    dend_thrsh_pav={"YPDCAFEIN40":0.008, "YPDSODIUMMETAARSENITE":0.02, "YPDBENOMYL500":0.01, "YPDCUSO410MM":0.09}
    dend_thrsh_cnv={"YPDCAFEIN40":0.0125, "YPDSODIUMMETAARSENITE":0.02, "YPDBENOMYL500":0.01, "YPDCUSO410MM":0.08}
    cbar_lim={"snp":{"YPDCAFEIN40":[-0.002, 0.002], "YPDSODIUMMETAARSENITE":[-0.02, 0.02], "YPDBENOMYL500":[-0.002,0.002], "YPDCUSO410MM":[-0.02, 0.02]}, # SHAP value limits for colorbar
              "pav":{"YPDCAFEIN40":[-0.005, 0.005], "YPDSODIUMMETAARSENITE":[-0.01, 0.01], "YPDBENOMYL500":[-0.005, 0.005], "YPDCUSO410MM":[-0.02, 0.02]},
              "cnv":{"YPDCAFEIN40":[-0.005, 0.005], "YPDSODIUMMETAARSENITE":[-0.01, 0.01], "YPDBENOMYL500":[-0.007, 0.007], "YPDCUSO410MM":[-0.004, 0.004]}}
    cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno, iso_colors, dend_thrsh_snp, cbar_lim=cbar_lim)
    for i in range(len(snp_files)):
        cluster_shap_indiv(f"Scripts/Data_Vis/SNP_Figures/SHAP/{snp_files[i]}",
                           target_envs, pheno, iso_colors, "snp", plot=True, dend_thrsh=dend_thrsh_snp, cbar_lim=cbar_lim)
        cluster_shap_indiv(f"Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_presence_absence/{pav_files[i]}",
                           target_envs, pheno, iso_colors, "pav", plot=True, dend_thrsh=dend_thrsh_pav, cbar_lim=cbar_lim)
        cluster_shap_indiv(f"Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_copy_number/{cnv_files[i]}",
                           target_envs, pheno, iso_colors, "cnv", plot=True, dend_thrsh=dend_thrsh_cnv, cbar_lim=cbar_lim)
