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

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/")

def sort_clusters(fitness, dend): # df_top, iso_col
    """
    Sort the dendrogram clusters by median fitness value
    """
    
    # Get the cluster median fitness
    # cluster_fitness = fitness.groupby(iso_col["Dend_Cluster"]).median().sort_values()
    # assert fitness.index.tolist() == dend["ivl"], "Fitness index does not match dendrogram leaves order. Sorting will not work correctly."
    cluster_fitness = fitness.groupby(dend["leaves_color_list"]).median().sort_values()
    
    # Sort the clusters by median fitness
    # iso_col_tmp = iso_col.copy(deep=True)
    # iso_col_tmp.insert(4, "idx", range(iso_col_tmp.shape[0])) # add index column
    fitness_sorted = {}
    # df_top_sorted = {}
    # dend_sorted = {k:[] for k in dend.keys()}
    # dict_keys(['icoord', 'dcoord', 'ivl', 'leaves', 'color_list', 'leaves_color_list'])
    map_df = pd.DataFrame({"ID": fitness.index.tolist(), "cluster": dend["leaves_color_list"]})
    
    for c in cluster_fitness.index:
        # fitness_sorted[c] = fitness[iso_col_tmp["Dend_Cluster"]==c]
        # df_top_sorted[c] = df_top.loc[:,iso_col_tmp["Dend_Cluster"]==c]
        fitness_sorted[c] = fitness[map_df.set_index("ID")["cluster"]==c]
        # df_top_sorted[c] = df_top.loc[:, map_df.set_index("ID")["cluster"]==c]
        
        # idx_lst = iso_col_tmp[iso_col_tmp["Dend_Cluster"]==c].idx
        # dend_sorted["ivl"].extend([dend["ivl"][i] for i in idx_lst])
        # dend_sorted["leaves"].extend([dend["leaves"][i] for i in idx_lst])
        # dend_sorted["leaves_color_list"].extend([dend["leaves_color_list"][i] for i in idx_lst])
        
        # idx_lst_trunc = [i for i in range(len(dend["color_list"])) if dend["color_list"][i]==c]
        # dend_sorted["icoord"].extend([dend["icoord"][i] for i in idx_lst_trunc])
        # dend_sorted["dcoord"].extend([dend["dcoord"][i] for i in idx_lst_trunc])
        # dend_sorted["color_list"].extend([dend["color_list"][i] for i in idx_lst_trunc])
    
    fitness_tmp = pd.DataFrame()
    # df_top_tmp = pd.DataFrame()
    sorted_colnames = []
    for c in fitness_sorted.keys():
        fitness_tmp = pd.concat([fitness_tmp, fitness_sorted[c]])
        # insert an extra 2 columns between clusters, for better visualization
        sorted_colnames.extend(fitness_sorted[c].index.tolist())
        sorted_colnames.extend([np.nan, np.nan]) # add two NaN columns for spacing
        # if df_top_tmp.empty:
        #     df_top_tmp = df_top_sorted[c]
        # else:
        #     df_top_tmp = pd.concat([df_top_tmp, pd.DataFrame(index=df_top_sorted[c].index, columns=[np.nan]*2)], axis=1)
        #     df_top_tmp = pd.concat([df_top_tmp, df_top_sorted[c]], axis=1)
            
    # iso_col_tmp = iso_col.loc[fitness_tmp.index,:]
    
    return(fitness_tmp, sorted_colnames) #df_top_tmp, iso_col_tmp, dend_sorted)


def plotTree(D_dendro, color_list, nrow, ncol, idx_ax_dend, save:bool, save_dir, name=""):
    """ 
    Plot a dendrogram
    Modified from: https://stackoverflow.com/questions/36538090/bigger-color-palette-in-matplotlib-for-scipys-dendrogram-python?noredirect=1#comment60681856_36538090
    """
    
    fig,ax = plt.subplots(nrows=nrow, ncols=ncol, figsize=(25, 5*nrow))
    icoord = np.array(D_dendro['icoord'])
    dcoord = np.array(D_dendro['dcoord'])
    leaves_color_list = np.array(D_dendro['leaves_color_list'])
    x_min, x_max = icoord.min(), icoord.max()
    y_min, y_max = dcoord.min(), dcoord.max()
    
    for xs, ys, color in zip(icoord, dcoord, leaves_color_list):
        ax[idx_ax_dend].plot(xs, ys,  color_list[color])
    
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
    
    return(fig,ax,color_list)


# Need to add the masking and cluster coloring to this function.
def cluster_shap_indiv(f, target_envs, pheno, iso_colors, plot:bool,
                       save_out:bool, dtype="", save_dir="",
                       dend_thrsh={'YPDCAFEIN40':0.01},
                       cbar_lim={'snp':{'YPDCAFEIN40':[-0.002, 0.002]}}):
    """ 
    Cluster SHAP values and color the dendrogram by clades for the 
    target environments
    """
    
    env = [e for e in f.split("/")[-1].split("_") if e in target_envs][0]
    if env in target_envs:
        # linreg_res = {} # store the linear regression results
        mwu_res = {"fitness":{}, "shap":{}} # store the mann-whitney U test results
        shap_v_fitness = {} # cluster shap vs median fitness linreg
        
        df = pd.read_csv(f"{f}", sep="\t", index_col=0).T # SHAP values dataframe
        
        # select the top 20 features based on median absolute shap value across isolates
        top_feat = df.abs().median(axis=1).sort_values(ascending=False)[:20].index
        df_top = df.loc[top_feat,:] #; print("Before clustering:", df_top.head())
        # print(env, df_top.T.describe().T) # stats across isolates for each feature
        print(dtype, env, df_top.abs().median(axis=1).min(), df_top.abs().median(axis=1).max())
        
        fitness = pheno.loc[df_top.columns, env] # order isolates
        # iso_col = iso_colors.loc[df_top.columns]
        
        if plot==True:
            print(f"Plotting {f}")
            
            # Cluster the isolates in the SHAP values dataframe
            Z_iso = hierarchy.linkage(df_top.T, method="ward", metric="euclidean")
            # og_leaves_list = hierarchy.leaves_list(Z_iso)
            Z_iso_ordered = hierarchy.optimal_leaf_ordering(Z_iso, df_top.T, metric="euclidean")
            optimal_leaves_list = hierarchy.leaves_list(Z_iso_ordered)
            dend = hierarchy.dendrogram(Z_iso_ordered, no_plot=True, color_threshold=dend_thrsh[env])
            print(f"Number of clusters: {len(np.unique(dend['leaves_color_list']))}")
            
            # Now, cluster the features in the SHAP values dataframe
            Z_feat = hierarchy.linkage(df_top, method="ward", metric="euclidean")
            Z_feat_ordered = hierarchy.optimal_leaf_ordering(Z_feat, df_top, metric="euclidean")
            dend_feat = hierarchy.dendrogram(Z_feat_ordered, no_plot=True)
            
            # Assign each dendrogram cluster a color
            cmap = mpl.colormaps.get_cmap("hsv").resampled(len(np.unique(dend["leaves_color_list"]))) # sample hsv color palette
            color_list = {np.unique(dend["leaves_color_list"])[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
            # iso_col["Dend_Cluster"] = dend["leaves_color_list"]
            # iso_col["Dend_Cluster_Color"] = iso_col["Dend_Cluster"].map(color_list) # assign colors to isolate
            
            # Plot the unsorted dendrogram with colored clusters
            hierarchy.set_link_color_palette(list(color_list.values())) # set the color palette for the dendrogram
            fig, ax = plt.subplots(3, 1, figsize=(20, 15))
            dend_iso = hierarchy.dendrogram(Z_iso_ordered, color_threshold=dend_thrsh[env],
                                 labels=df_top.iloc[:,optimal_leaves_list].columns, ax=ax[0],
                                 leaf_font_size=2.5)
            '''
            # dend_iso["leaves"]==optimal_leaves_list # sanity check: True!
            print("BEFORE", color_list)
            new_color_list = {np.unique(dend_iso["leaves_color_list"])[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)}
            print("AFTER", new_color_list)
            color_list = color_list.update(new_color_list)
            
            # Plot the SHAP values heatmap
            # adjust SHAP values based on given limits passed to cbar_lim
            df_top.clip(lower=cbar_lim[dtype][env][0], upper=cbar_lim[dtype][env][1], inplace=True)
            sns.heatmap(df_top.iloc[dend_feat["leaves"], dend_iso["leaves"]],
                        cmap="RdBu_r", center=0,
                        cbar_kws={"orientation": "vertical", "location":"right"},
                        yticklabels=True, xticklabels=False, ax=ax[1])
            
            # Plot the violin plot of fitness values by dendrogram cluster
            fitness_toplot = pd.DataFrame({
                "Fitness": fitness.loc[df_top.iloc[:,optimal_leaves_list].columns],
                "Cluster_Color": dend_iso["leaves_color_list"]}).reset_index()
            cluster_medians = fitness_toplot.groupby("Cluster_Color").\
                median("Fitness").sort_values(by="Fitness", ascending=True) # median fitness
            fitness_toplot["Cluster_Color"] = pd.Categorical(
                fitness_toplot["Cluster_Color"], categories=cluster_medians.index,
                ordered=True) # convert cluster color into an ordered categorical variable
            fitness_toplot = fitness_toplot.sort_values("Cluster_Color") # sort by median cluster fitness
            
            sns.violinplot(data=fitness_toplot, x="Cluster_Color", y="Fitness",
                hue="Cluster_Color", palette=color_list, inner="quart", ax=ax[2])
            
            plt.savefig(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_dendrogram_unsorted_v4.pdf")
            plt.close()
            '''
            # Re-order the SHAP values and fitness values based on the dendrogram
            df_top = df_top.iloc[dend_feat["leaves"], dend_iso["leaves"]] #; print("After clustering:", df_top.head())
            df_top.to_csv(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_sorted_v4.tsv", sep="\t")
            # fitness = fitness.loc[df_top.columns]
            # iso_col = iso_colors.loc[df_top.columns,:]
            
            '''
            # Conduct a mann-whitney U test to compare fitness distributions and shap values between clusters
            mannwhitney_res = {} # save fitness distribution comparisons
            mwu_res_shap = {} # save shap value comparisons
            clusters = cluster_medians.index.to_list() ### FIX THISSSS
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    # compare fitness distributions
                    group1 = fitness_toplot[fitness_toplot['Cluster_Color']==clusters[i]]['Fitness']
                    group2 = fitness_toplot[fitness_toplot['Cluster_Color']==clusters[j]]['Fitness']
                    stat, p_value = mannwhitneyu(group1, group2)
                    mannwhitney_res[(clusters[i], clusters[j])] = {'stat': stat, 'pval': p_value}
                    
                    # compare absolute shap values
                    # group1_shap = df_top[df_top.columns[fitness_toplot['Cluster_Color']==clusters[i]]]
                    group1_shap = df_top.loc[:,fitness_toplot.loc[fitness_toplot['Cluster_Color']==clusters[i], "ID"]]
                    group1_shap = group1_shap.abs().to_numpy().flatten()
                    # group2_shap = df_top[df_top.columns[fitness_toplot['Cluster_Color']==clusters[j]]]
                    group2_shap = df_top.loc[:,fitness_toplot.loc[fitness_toplot['Cluster_Color']==clusters[j], "ID"]]
                    group2_shap = group2_shap.abs().to_numpy().flatten()
                    stat_shap, p_value_shap = mannwhitneyu(group1_shap, group2_shap)
                    mwu_res_shap[(clusters[i], clusters[j])] = {'stat': stat_shap, 'pval': p_value_shap}
            mwu_res["fitness"][env] = mannwhitney_res
            mwu_res["shap"][env] = mwu_res_shap
            '''
            #
            # Determine the relationship between cluster fitness and absolute shap values
            # tmp = df_top.T.abs().merge(fitness_toplot.set_index("ID")["Cluster_Color"], 
            #              left_index=True, right_index=True, how="inner")
            # tmp.groupby("Cluster_Color").median().median(axis=1) # for sanity check purposes. The line below works as expected.
            shap_med_per_cluster = df_top.T.abs().groupby(fitness_toplot.set_index("ID")["Cluster_Color"]).median().median(axis=1)
            m, b, r_value, p_value, std_err = \
                linregress(shap_med_per_cluster.values, cluster_medians.values.flatten())
            shap_v_fitness[env] = {'m':m, 'b':b, 'r':r_value, 'pval':p_value, 'SE':std_err}
            del m, b, r_value, p_value, std_err
            
            # # Trend in median fitness of clusters
            # print(cluster_medians.reset_index().index)
            # m, b, r_value, p_value, std_err = \
            # linregress(cluster_medians.reset_index().index, cluster_medians['Fitness']) # calculate R2 and p-value
            # linreg_res[env] = {'m':m, 'b':b, 'r':r_value, 'pval':p_value, 'SE':std_err}
            # del m, b, r_value, p_value, std_err
            
            # # Sort the dendrogram clusters by median fitness value and plot (don't save to file yet)
            # fitness_sorted, iso_col_sorted, df_top_sorted, dend_sorted = sort_clusters(dend, df_top, fitness, iso_col) # dend,  #, dend_sorted
            
            # if save_out==True:
            #     print(f"Saving sorted shap values, fitness values, and isolates to files...")
            #     df_top_sorted.to_csv(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_sorted_v3.tsv", sep="\t")
            #     iso_col_sorted.to_csv(f"{save_dir}/isolates_clustered_{env}_top20_{dtype}_sorted_v3.tsv", sep="\t")
            
            # # plotTree(dend, color_list, 2, 1, 0, True, save_dir, f"SHAP_clustered_{env}_top20_{dtype}_dendrogram_unsorted_v3") ; plt.close() # save unsorted dendrogram
            # fig, ax_dend, fig_clist = plotTree(dend_sorted, color_list, 3, 1, 0, True, save_dir, f"SHAP_clustered_{env}_top20_{dtype}_dendrogram_sorted_v3") # save the sorted dendrogram (since icoord and dcoord are not recalculated, it's the same as the unsorted dendrogram)
            # print("FIG COLOR LIST:", fig_clist)
            
            # # Fit a line to the median fitness values in each cluster to get the pearson r and p-value
            # grouped_data = pd.concat([iso_col_sorted["Dend_Cluster"], fitness_sorted[0]], axis=1, ignore_index=False)
            # grouped_data.columns = ['Dend_Cluster', 'Fitness']
            # grouped_data_med = grouped_data.groupby('Dend_Cluster').median().reset_index() # median fitness values per cluster
            
            # # Trend in median fitness of clusters
            # m, b, r_value, p_value, std_err = \
            # linregress(grouped_data_med.index, grouped_data_med['Fitness']) # calculate R2 and p-value
            # linreg_res[env] = {'m':m, 'b':b, 'r':r_value, 'pval':p_value, 'SE':std_err}
            # del m, b, r_value, p_value, std_err
            
            # # Conduct a mann-whitney U test to compare fitness distributions and shap values between clusters
            # mannwhitney_res = {} # save fitness distribution comparisons
            # mwu_res_shap = {} # save shap value comparisons
            # clusters = grouped_data['Dend_Cluster'].unique()
            # for i in range(len(clusters)):
            #     for j in range(i + 1, len(clusters)):
            #         # compare fitness distributions
            #         group1 = grouped_data[grouped_data['Dend_Cluster']==clusters[i]]['Fitness']
            #         group2 = grouped_data[grouped_data['Dend_Cluster']==clusters[j]]['Fitness']
            #         stat, p_value = mannwhitneyu(group1, group2)
            #         mannwhitney_res[(clusters[i], clusters[j])] = {'stat': stat, 'pval': p_value}
                    
            #         # compare absolute shap values
            #         group1_shap = df_top_sorted[df_top_sorted.columns[iso_col_sorted['Dend_Cluster']==clusters[i]]]
            #         group1_shap = group1_shap.abs().to_numpy().flatten()
            #         group2_shap = df_top_sorted[df_top_sorted.columns[iso_col_sorted['Dend_Cluster']==clusters[j]]]
            #         group2_shap = group2_shap.abs().to_numpy().flatten()
            #         stat_shap, p_value_shap = mannwhitneyu(group1_shap, group2_shap)
            #         mwu_res_shap[(clusters[i], clusters[j])] = {'stat': stat_shap, 'pval': p_value_shap}
            # mwu_res["fitness"][env] = mannwhitney_res
            # mwu_res["shap"][env] = mwu_res_shap
            
            # # Determine the relationship between cluster fitness and absolute shap values
            # shap_med_per_cluster = df_top_sorted.T.abs().groupby(iso_col_sorted["Dend_Cluster"]).median().median(axis=1)
            # fitness_med_per_cluster = fitness_sorted.groupby(iso_col["Dend_Cluster"]).median()
            # m, b, r_value, p_value, std_err = \
            #     linregress(shap_med_per_cluster.values, fitness_med_per_cluster.values.flatten())
            # shap_v_fitness[env] = {'m':m, 'b':b, 'r':r_value, 'pval':p_value, 'SE':std_err}
            
            # # Add masked values to dataframe to visualize clusters in heatmap
            # mask_top = {}
            # mask_iso_col = {}
            # for i in range(iso_col_sorted.shape[0]):
            #     mask_top[df_top_sorted.columns[i]] = df_top_sorted.iloc[:,i]
            #     mask_iso_col[iso_col_sorted.index[i]] = iso_col_sorted.iloc[i,:]
                
            #     if i==iso_col_sorted.shape[0]-1:
            #         break
            #     elif iso_col_sorted.iloc[i,2]!=iso_col_sorted.iloc[i+1,2]: # cluster divider mask
            #         mask_top[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=df_top_sorted.index)
            #         mask_iso_col[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series("#ffffff", index=iso_col_sorted.columns)
            
            # mask_top = pd.DataFrame(mask_top)
            # mask_iso_col = pd.DataFrame(mask_iso_col).T
            
            # # Adjust SHAP values based on given limits passed to cbar_lim
            # mask_top.clip(lower=cbar_lim[dtype][env][0], upper=cbar_lim[dtype][env][1], inplace=True)
            
            # # Plot
            # sns.heatmap(mask_top, cmap="RdBu_r", center=0,
            #             cbar_kws={"orientation": "vertical", "location":"right"},
            #             yticklabels=True, xticklabels=False, ax=ax_dend[1])
            
            # # Create an additional axis for the column colors
            # ax_dend[1].tick_params(axis="x", which="major", pad=20, length=0) # extra padding for column colors
            # for i in range(mask_iso_col.shape[0]):
            #     ax_dend[1].add_patch(plt.Rectangle(xy=(i,-.4), width=0.016, height=0.5, color=mask_iso_col["Clade_Color"][i]))
            
            # sns.violinplot(x=iso_col_sorted["Dend_Cluster"], y=fitness_sorted[0], hue=iso_col_sorted["Dend_Cluster"], palette=color_list, inner="quart", ax=ax_dend[2])
            # fig.savefig(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_sorted_v3.pdf")
            # plt.close()
            # print(f"{save_dir}/SHAP_clustered_{env}_top20_{dtype}_sorted_v3.pdf")
            
        return df_top, mwu_res, shap_v_fitness # linreg_res


def cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno, iso_colors, dend_thrsh, cbar_lim, save_dir, by="snp"):
    '''Cluster SHAP values based on the SNP model and re-order the PAV and CNV models'''
    correlations = {} # store the shap vs fitness linear regression results
    for env in target_envs:
        # if env == "YPDBENOMYL500":
        print("Clustering SHAP values for", env)
        
        # Get file and read in SHAP values; get the top 20 features
        f_snp = [f"{f}" for f in snp_files if [f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        f_pav = [f"{f}" for f in pav_files if [f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        f_cnv = [f"{f}" for f in cnv_files if [f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        top_snp, _, _ = cluster_shap_indiv(f_snp, target_envs, pheno, iso_colors, plot=False, save_out=False)
        top_pav, _, _ = cluster_shap_indiv(f_pav, target_envs, pheno, iso_colors, plot=False, save_out=False)
        top_cnv, _, _ = cluster_shap_indiv(f_cnv, target_envs, pheno, iso_colors, plot=False, save_out=False)
        
        if by=="snp":
            # Cluster isolates by SNP SHAP values
            Z_iso = hierarchy.linkage(top_snp.T, method="ward", metric="euclidean")
            Z_iso_ordered = hierarchy.optimal_leaf_ordering(Z_iso, top_snp.T, metric="euclidean")
        if by=="pav":
            # Cluster isolates by PAV SHAP values
            Z_iso = hierarchy.linkage(top_pav.T, method="ward", metric="euclidean")
            Z_iso_ordered = hierarchy.optimal_leaf_ordering(Z_iso, top_pav.T, metric="euclidean")
        if by=="cnv":
            # Cluster isolates by CNV SHAP values
            Z_iso = hierarchy.linkage(top_cnv.T, method="ward", metric="euclidean")
            Z_iso_ordered = hierarchy.optimal_leaf_ordering(Z_iso, top_cnv.T, metric="euclidean")
        
        optimal_leaves_list = hierarchy.leaves_list(Z_iso_ordered)
        dend = hierarchy.dendrogram(Z_iso_ordered, no_plot=True, color_threshold=dend_thrsh[env])
        
        # Cluster features in the SNP, PAV, and CNV SHAP value dataframe
        Z_feat_snp = hierarchy.linkage(top_snp, method="ward", metric="euclidean")
        Z_feat_snp_ordered = hierarchy.optimal_leaf_ordering(Z_feat_snp, top_snp, metric="euclidean")
        dend_feat_snp = hierarchy.dendrogram(Z_feat_snp_ordered, no_plot=True)
        
        Z_feat_pav = hierarchy.linkage(top_pav, method="ward", metric="euclidean")
        Z_feat_pav_ordered = hierarchy.optimal_leaf_ordering(Z_feat_pav, top_pav, metric="euclidean")
        dend_feat_pav = hierarchy.dendrogram(Z_feat_pav_ordered, no_plot=True)
        
        Z_feat_cnv = hierarchy.linkage(top_cnv, method="ward", metric="euclidean")
        Z_feat_cnv_ordered = hierarchy.optimal_leaf_ordering(Z_feat_cnv, top_cnv, metric="euclidean")
        dend_feat_cnv = hierarchy.dendrogram(Z_feat_cnv_ordered, no_plot=True)
        
        # Assign each dendrogram cluster a color
        cmap = mpl.colormaps.get_cmap("hsv").resampled(len(np.unique(dend["leaves_color_list"]))) # sample hsv color palette
        color_list = {np.unique(dend["leaves_color_list"])[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)} # assign colors to clades
        
        # Plot the unsorted dendrogram with colored clusters
        hierarchy.set_link_color_palette(list(color_list.values())) # set the color palette for the dendrogram
        fig, ax = plt.subplots(5, 1, figsize=(20, 20))
        if by=="snp":
            dend_iso = hierarchy.dendrogram(Z_iso_ordered, ax=ax[0],
                color_threshold=dend_thrsh[env], leaf_font_size=2.5,
                labels=top_snp.iloc[:, optimal_leaves_list].columns)
        if by=="pav":
            dend_iso = hierarchy.dendrogram(Z_iso_ordered, ax=ax[0],
                color_threshold=dend_thrsh[env], leaf_font_size=2.5,
                labels=top_pav.iloc[:, optimal_leaves_list].columns)
        if by=="cnv":
            dend_iso = hierarchy.dendrogram(Z_iso_ordered, ax=ax[0],
                color_threshold=dend_thrsh[env], leaf_font_size=2.5,
                labels=top_cnv.iloc[:, optimal_leaves_list].columns)
        print("BEFORE", color_list)
        new_color_list = {np.unique(dend_iso["leaves_color_list"])[i]:rgb2hex(cmap(i)[:3]) for i in range(cmap.N)}
        print("AFTER", new_color_list)
        color_list = color_list.update(new_color_list)
        iso_colors["Clade_Color"] = iso_colors["Clade_Color"].map(new_color_list) # assign colors to isolate
        
        # Re-order the SNP SHAP values dataframe and fitness values
        top_snp = top_snp.iloc[dend_feat_snp["leaves"], dend_iso["leaves"]]
        top_pav = top_pav.iloc[dend_feat_pav["leaves"], dend_iso["leaves"]]
        top_cnv = top_cnv.iloc[dend_feat_cnv["leaves"], dend_iso["leaves"]]
        fitness = pheno.loc[top_snp.columns, env]
        
        # Sort the clusters by median fitness value and reorder all shap dataframes
        fitness_sorted, sorted_colnames = sort_clusters(fitness, dend_iso)
        no_na_cols = [col for col in sorted_colnames if not pd.isna(col)] # get the indices of the non-NaN columns
        top_snp_sorted = top_snp.loc[:, no_na_cols] # reorder the SNP dataframe
        top_pav_sorted = top_pav.loc[:, no_na_cols]
        top_cnv_sorted = top_cnv.loc[:, no_na_cols]
        na_cols = [i for i in range(len(sorted_colnames)) if pd.isna(sorted_colnames[i])] # get the indices of the extra columns
        na_cols.pop(-1)
        na_cols.pop(-1) # remove the last two NaN columns, which are not needed for the heatmap
        
        # Add the masked values to the remaining 2 dataframes to visualize clusters in heatmap
        for i in na_cols:
            top_snp_sorted.insert(i, f"NaN{i}", pd.Series(np.nan, index=top_snp_sorted.index))
            top_pav_sorted.insert(i, f"NaN{i}", pd.Series(np.nan, index=top_pav_sorted.index))
            top_cnv_sorted.insert(i, f"NaN{i}", pd.Series(np.nan, index=top_cnv_sorted.index))
        
        # iso_col = iso_colors.loc[top_snp.columns]
        
        top_snp_sorted.clip(lower=cbar_lim["snp"][env][0], upper=cbar_lim["snp"][env][1], inplace=True)
        sns.heatmap(top_snp_sorted, cmap="RdBu_r", center=0, ax=ax[1],
                    yticklabels=True, xticklabels=False,
                    cbar_kws={"orientation": "vertical", "location":"right"})
        
        top_pav_sorted.clip(lower=cbar_lim["pav"][env][0], upper=cbar_lim["pav"][env][1], inplace=True)
        sns.heatmap(top_pav_sorted, cmap="RdBu_r", center=0, ax=ax[2],
                    yticklabels=True, xticklabels=False,
                    cbar_kws={"orientation": "vertical", "location":"right"})
        
        top_cnv_sorted.clip(lower=cbar_lim["cnv"][env][0], upper=cbar_lim["cnv"][env][1], inplace=True)
        sns.heatmap(top_cnv_sorted, cmap="RdBu_r", center=0, ax=ax[3],
                    yticklabels=True, xticklabels=False,
                    cbar_kws={"orientation": "vertical", "location":"right"})
        
        # Plot the violin plot of fitness values by dendrogram cluster
        fitness_toplot = pd.DataFrame({
            "Fitness": fitness,
            "Cluster_Color": dend_iso["leaves_color_list"]}).reset_index()
        cluster_medians = fitness_toplot.groupby("Cluster_Color").\
            median("Fitness").sort_values(by="Fitness", ascending=True) # median fitness
        fitness_toplot["Cluster_Color"] = pd.Categorical(
            fitness_toplot["Cluster_Color"], categories=cluster_medians.index,
            ordered=True) # convert cluster color into an ordered categorical variable
        fitness_toplot = fitness_toplot.sort_values("Cluster_Color") # sort by median cluster fitness
        
        # try:
        sns.violinplot(data=fitness_toplot, x="Cluster_Color", y="Fitness",
            hue="Cluster_Color", palette=color_list, inner="quart", ax=ax[4])
        # except ValueError:
        #     for color in fitness_toplot["Cluster_Color"].unique():
        #         if (color not in color_list.keys()) and (color not in color_list.values()):
        #             color_list[color] = color
        #         elif (color not in color_list.keys()) and (color in color_list.values()):
        #             new_color = 
        #     sns.violinplot(data=fitness_toplot, x="Cluster_Color", y="Fitness",
        #         hue="Cluster_Color", palette=color_list, inner="quart", ax=ax[4])
        
        plt.savefig(f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_sorted_v6.pdf")
        plt.close()
        
        # Save the sorted SHAP values and fitness values to files
        top_snp_sorted.to_csv(f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_snp_sorted_v6.tsv", sep="\t")
        top_pav_sorted.to_csv(f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_pav_sorted_v6.tsv", sep="\t")
        top_cnv_sorted.to_csv(f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_cnv_sorted_v6.tsv", sep="\t")
        fitness_toplot.to_csv(f"{save_dir}/SHAP_clustered_by_{by}_{env}_fitness_sorted_v6.tsv", sep="\t")
        
        # Determine the relationship between cluster fitness and absolute shap values
        if by=="snp":
            top = top_snp_sorted.copy(deep=True)
        elif by=="pav":
            top = top_pav_sorted.copy(deep=True)
        elif by=="cnv":
            top = top_cnv_sorted.copy(deep=True)
        
        for i in range(top.shape[0]): # relationship between each feature and fitness
            shap_row = top.iloc[i, :].loc[no_na_cols] # 625 SHAP values for feature i
            print(shap_row, fitness_sorted)
            m, b, r, p, se = linregress(shap_row, fitness_sorted.values.flatten())
            correlations[f'sorted_by_{by}:{env}:{top.index[i]}'] = \
                {'m':m, 'b':b, 'r':r, 'p-value':p, 'se':se}
            del m, b, r, p, se
        # relationship between median shap values and fitness
        print(np.median(np.abs(top.loc[:, no_na_cols]), axis=0))
        m, b, r, p, se = linregress(np.median(np.abs(top.loc[:, no_na_cols]), axis=0), 
                                    fitness_sorted.values.flatten())
        correlations[f'sorted_by_{by}:{env}:median_SHAP'] = \
            {'m':m, 'b':b, 'r':r, 'p-value':p, 'se':se}
        shap_v_fitness_df = pd.DataFrame.from_dict(correlations, orient='index')
        # shap_med_per_cluster = df_top.T.abs().groupby(fitness_toplot.set_index("ID")["Cluster_Color"]).median().median(axis=1)
        # m, b, r_value, p_value, std_err = \
        #     linregress(shap_med_per_cluster.values, cluster_medians.values.flatten())
        # shap_v_fitness[env] = {'m':m, 'b':b, 'r':r_value, 'pval':p_value, 'SE':std_err}
        
        
        # # Sort the dendrogram clusters by median fitness value and reorder all dataframes by SNP-based isolate order
        # fitness_sorted, iso_col_sorted, top_snp_sorted, dend_snp_sorted = sort_clusters(dend_snp, top_snp, fitness, iso_col) # dend_snp,  # , dend_snp_sorted 
        # plotTree(dend_snp, color_list, 2, 1, 0, True, save_dir, f"SHAP_clustered_by_snp_{env}_top20_dendrogram_v3") ; plt.close()
        # top_pav = top_pav.iloc[:,dend_snp_sorted["leaves"]] # dend_snp_sorted
        # top_cnv = top_cnv.iloc[:,dend_snp_sorted["leaves"]] # dend_snp_sorted
        
        # # Add masked values to dataframe to visualize clusters in heatmap
        # mask_top_snp = {}
        # mask_top_pav = {}
        # mask_top_cnv = {}
        # mask_iso_col = {}
        # for i in range(iso_col_sorted.shape[0]):
        #     mask_top_snp[top_snp_sorted.columns[i]] = top_snp_sorted.iloc[:,i]
        #     mask_top_pav[top_pav.columns[i]] = top_pav.iloc[:,i]
        #     mask_top_cnv[top_cnv.columns[i]] = top_cnv.iloc[:,i]
        #     mask_iso_col[iso_col_sorted.index[i]] = iso_col_sorted.iloc[i,:]
            
        #     if i==624:
        #         break
        #     elif iso_col_sorted.iloc[i,2]!=iso_col_sorted.iloc[i+1,2]: # cluster division
        #         mask_top_snp[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=top_snp_sorted.index)
        #         mask_top_cnv[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=top_cnv.index)
        #         mask_top_pav[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series(np.nan, index=top_pav.index)
        #         mask_iso_col[f"{iso_col_sorted.iloc[i,2]}_{iso_col_sorted.iloc[i+1,2]}_mask"] = pd.Series("#ffffff", index=iso_col_sorted.columns)
        
        # mask_top_snp = pd.DataFrame(mask_top_snp)
        # mask_top_pav = pd.DataFrame(mask_top_pav)
        # mask_top_cnv = pd.DataFrame(mask_top_cnv)
        # mask_iso_col = pd.DataFrame(mask_iso_col).T
        
        # # Adjust SHAP values based on given limits passed to cbar_lim
        # mask_top_snp.clip(lower=cbar_lim["snp"][env][0], upper=cbar_lim["snp"][env][1], inplace=True)
        # mask_top_pav.clip(lower=cbar_lim["pav"][env][0], upper=cbar_lim["pav"][env][1], inplace=True)
        # mask_top_cnv.clip(lower=cbar_lim["cnv"][env][0], upper=cbar_lim["cnv"][env][1], inplace=True)
        
        # # Plot the SNP heatmap
        # fig, ax = plt.subplots(4,1, figsize=(13,8.5))
        # s = sns.heatmap(mask_top_snp, cmap="RdBu_r", center=0,
        #                 cbar_kws={"orientation": "vertical", "location":"right"},
        #                 yticklabels=True, xticklabels=False, ax=ax[0])
        
        # # # Create an additional axis for the column colors
        # # ax[0].tick_params(axis="x", which="major", pad=20, length=0) # extra padding for column colors
        # # for i in range(iso_col_sorted.shape[0]):
        # #     ax[0].add_patch(plt.Rectangle(xy=(i,-.4), width=0.0032, height=0.5, color=mask_iso_col["Clade_Color"][i]))
        # #     # ax[0].add_patch(plt.Rectangle(xy=(i,-.2), width=0.0064, height=0.6, color=mask_iso_col["Dend_Cluster_Color"][i]))
        
        # # Plot the other two data types
        # p = sns.heatmap(mask_top_pav, cmap="RdBu_r", center=0,
        #                 cbar_kws={"orientation": "vertical", "location":"right"},
        #                 yticklabels=True, xticklabels=False, ax=ax[1])
        # c = sns.heatmap(mask_top_cnv, cmap="RdBu_r", center=0,
        #                 cbar_kws={"orientation": "vertical", "location":"right"},
        #                 yticklabels=True, xticklabels=False, ax=ax[2])
        
        # # Group fitness by dendrogram cluster and plot violin
        # v = sns.violinplot(x=iso_col_sorted["Dend_Cluster"],
        #                 y=fitness_sorted[0],
        #                 hue=iso_col_sorted["Dend_Cluster"],
        #                 palette=color_list, inner="quart", ax=ax[3])
        # plt.tight_layout()
        # plt.savefig(f"{save_dir}/SHAP_clustered_by_snp_{env}_top20_sorted_v3.pdf")
        # plt.close()
        # print(f"Saved as {save_dir}/SHAP_clustered_by_snp_{env}_top20_sorted_v3.pdf")
        
        del top_snp, top_pav, top_cnv, fitness #, iso_col, mask_iso_col
    return shap_v_fitness_df


def save_results(mwu_res, shap_v_fitness_res, mwu_save, shap_v_fitness_save):
    # # Convert linear regression results to a dataframe
    # linreg_df = pd.DataFrame.from_dict({(i,j,k): linreg_res[i][j][k]
    #                                      for i in linreg_res.keys()
    #                                      for j in linreg_res[i].keys()
    #                                      for k in linreg_res[i][j].keys()},
    #                                     orient='index')
    # linreg_df = linreg_df.droplevel(0)
    # linreg_df.sort_index(inplace=True)
    # linreg_df.to_csv(linreg_save)
    
    # Convert Mann-Whitney U test results to a dataframe
    mwu_df = pd.DataFrame.from_dict({(i,j,k,l,m): mwu_res[i][j][k][l][m]
                                      for i in mwu_res.keys()
                                      for j in mwu_res[i].keys()
                                      for k in mwu_res[i][j].keys()
                                      for l in mwu_res[i][j][k].keys()
                                      for m in mwu_res[i][j][k][l].keys()},
                                     orient='index')
    mwu_df = mwu_df.droplevel(0)
    mwu_df.sort_index(inplace=True)
    mwu_df.reset_index(level=[3], inplace=True)
    mwu_df.insert(0, "Cluster1", [i[0] for i in mwu_df.level_3])
    mwu_df.insert(1, "Cluster2", [i[1] for i in mwu_df.level_3])
    mwu_df.drop(columns="level_3").to_csv(mwu_save)
    
    # Median absolute SHAP values vs median fitness of clusters linear regression
    shap_v_fitness_df = pd.DataFrame.from_dict({(i,j,k): shap_v_fitness_res[i][j][k]
                                         for i in shap_v_fitness_res.keys()
                                         for j in shap_v_fitness_res[i].keys()
                                         for k in shap_v_fitness_res[i][j].keys()},
                                        orient='index')
    shap_v_fitness_df = shap_v_fitness_df.droplevel(0)
    shap_v_fitness_df.sort_index(inplace=True)
    shap_v_fitness_df.to_csv(shap_v_fitness_save)


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

    # ############### Read in SHAP value files for baseline models ###############
    # # Note: plots saved to /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_6/baseline_top20/
    
    # snp_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/SNP/baseline"
    # files = os.listdir(snp_path)
    # snp_files = [f"{snp_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    # pav_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/PAV/baseline"
    # files = os.listdir(pav_path)
    # pav_files = [f"{pav_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    # cnv_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/CNV/baseline"
    # files = os.listdir(cnv_path)
    # cnv_files = [f"{cnv_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]

    # target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    # snp_files = [f for f in snp_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    # pav_files = [f for f in pav_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    # cnv_files = [f for f in cnv_files if len([f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    
    # # Cluster SHAP values
    # dend_thrsh_snp={"YPDCAFEIN40":0.002, "YPDCAFEIN50":0.002, "YPDSODIUMMETAARSENITE":0.01, "YPDBENOMYL500":0.002, "YPDCUSO410MM":0.04}
    # dend_thrsh_pav={"YPDCAFEIN40":0.008, "YPDCAFEIN50":0.008, "YPDSODIUMMETAARSENITE":0.01, "YPDBENOMYL500":0.01, "YPDCUSO410MM":0.04}
    # dend_thrsh_cnv={"YPDCAFEIN40":0.01, "YPDCAFEIN50":0.01, "YPDSODIUMMETAARSENITE":0.01, "YPDBENOMYL500":0.01, "YPDCUSO410MM":0.02}
    # cbar_lim={"snp":{"YPDCAFEIN40":[-0.002, 0.002], "YPDCAFEIN50":[-0.002, 0.002], "YPDSODIUMMETAARSENITE":[-0.02, 0.02], "YPDBENOMYL500":[-0.002,0.002], "YPDCUSO410MM":[-0.02, 0.02]}, # SHAP value limits for colorbar
    #           "pav":{"YPDCAFEIN40":[-0.005, 0.005], "YPDCAFEIN50":[-0.005, 0.005], "YPDSODIUMMETAARSENITE":[-0.01, 0.01], "YPDBENOMYL500":[-0.005, 0.005], "YPDCUSO410MM":[-0.03, 0.03]},
    #           "cnv":{"YPDCAFEIN40":[-0.005, 0.005], "YPDCAFEIN50":[-0.005, 0.005], "YPDSODIUMMETAARSENITE":[-0.01, 0.01], "YPDBENOMYL500":[-0.007, 0.007], "YPDCUSO410MM":[-0.005, 0.005]}}
    
    # cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno,
    #                       iso_colors, dend_thrsh_snp, cbar_lim=cbar_lim,
    #                       save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
    
    # linreg_res = {}
    # mwu_res = {}
    # for i in range(len(snp_files)):
    #     _, snp_lin, snp_mwu, snp_shap_v_fitness = cluster_shap_indiv(f"{snp_files[i]}",
    #         target_envs, pheno, iso_colors, "snp", plot=True, save_out=False,
    #         dend_thrsh=dend_thrsh_snp, cbar_lim=cbar_lim,
    #         save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
    #     _, pav_lin, pav_mwu, pav_shap_v_fitness = cluster_shap_indiv(f"{pav_files[i]}",
    #         target_envs, pheno, iso_colors, "pav", plot=True, save_out=False,
    #         dend_thrsh=dend_thrsh_pav, cbar_lim=cbar_lim,
    #         save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
    #     _, cnv_lin, cnv_mwu, cnv_shap_v_fitness = cluster_shap_indiv(f"{cnv_files[i]}",
    #         target_envs, pheno, iso_colors, "cnv", plot=True, save_out=False,
    #         dend_thrsh=dend_thrsh_cnv, cbar_lim=cbar_lim,
    #         save_dir="Scripts/Data_Vis/Section_6/shap_clusters/baseline_top20")
        
    #     linreg_res[i] = {"snp":snp_lin, "pav":pav_lin, "cnv":cnv_lin}
    #     mwu_res[i] = {"snp":snp_mwu, "pav":pav_mwu, "cnv":cnv_mwu}

    # linreg_save = "Scripts/Data_Vis/Section_6/baseline_top20/SHAP_cluster_median_fitness_linreg_RF_baseline_models.csv"
    # mwu_save = "Scripts/Data_Vis/Section_6/baseline_top20/SHAP_cluster_fitness_mwu_RF_baseline_models.csv"
    # save_results(linreg_res, mwu_res, linreg_save, mwu_save)

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

    dend_thrsh_snp={"YPDCAFEIN40":0.014, "YPDCAFEIN50":0.035, "YPDSODIUMMETAARSENITE":0.05, "YPDBENOMYL500":0.012, "YPDCUSO410MM":0.09}
    dend_thrsh_pav={"YPDCAFEIN40":0.08, "YPDCAFEIN50":0.055, "YPDSODIUMMETAARSENITE":0.08, "YPDBENOMYL500":0.09, "YPDCUSO410MM":0.2}
    dend_thrsh_cnv={"YPDCAFEIN40":0.1, "YPDCAFEIN50":0.07, "YPDSODIUMMETAARSENITE":0.25, "YPDBENOMYL500":0.08, "YPDCUSO410MM":0.8}
    cbar_lim={"snp":{"YPDCAFEIN40":[-0.002, 0.002], "YPDCAFEIN50":[-0.005, 0.005], "YPDSODIUMMETAARSENITE":[-0.02, 0.02], "YPDBENOMYL500":[-0.001,0.001], "YPDCUSO410MM":[-0.03, 0.03]}, # SHAP value limits for colorbar
                "pav":{"YPDCAFEIN40":[-0.04, 0.04], "YPDCAFEIN50":[-0.01, 0.01], "YPDSODIUMMETAARSENITE":[-0.04, 0.04], "YPDBENOMYL500":[-0.01, 0.01], "YPDCUSO410MM":[-0.04, 0.04]},
                "cnv":{"YPDCAFEIN40":[-0.04, 0.04], "YPDCAFEIN50":[-0.02, 0.02], "YPDSODIUMMETAARSENITE":[-0.06, 0.06], "YPDBENOMYL500":[-0.01, 0.01], "YPDCUSO410MM":[-0.5, 0.5]}}

    shap_v_fit_snp = cluster_shap_combined(snp_files, pav_files, cnv_files,
        target_envs, pheno, iso_colors, dend_thrsh_snp, cbar_lim=cbar_lim,
        save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20", by="snp")
    shap_v_fit_pav = cluster_shap_combined(snp_files, pav_files, cnv_files,
        target_envs, pheno, iso_colors, dend_thrsh_pav, cbar_lim=cbar_lim,
        save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20", by="pav")
    shap_v_fit_cnv = cluster_shap_combined(snp_files, pav_files, cnv_files,
        target_envs, pheno, iso_colors, dend_thrsh_cnv, cbar_lim=cbar_lim,
        save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20", by="cnv")

    # Save the SHAP values vs fitness linear regression results
    shap_v_fitness = pd.concat([shap_v_fit_snp, shap_v_fit_pav, shap_v_fit_cnv], axis=0)
    shap_v_fitness.index = shap_v_fitness.index.str.split(":", expand=True)
    shap_v_fitness.index = shap_v_fitness.index.set_names(["Sorting Method", "Environment", "Feature"])
    shap_v_fitness.to_csv(
        "Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_cluster_values_vs_fitness_linreg_RF_FS_models_v5.csv")
    
    # mwu_res = {}
    # shap_v_fitness = {}
    # for i in range(len(snp_files)):
    #     _, snp_mwu, snp_shap_v_fitness= cluster_shap_indiv(f"{snp_files[i]}",
    #         target_envs, pheno, iso_colors, "snp", plot=True, save_out=True,
    #         dend_thrsh=dend_thrsh_snp, cbar_lim=cbar_lim,
    #         save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20")
    #     _, pav_mwu, pav_shap_v_fitness = cluster_shap_indiv(f"{pav_files[i]}",
    #         target_envs, pheno, iso_colors, "pav", plot=True, save_out=True,
    #         dend_thrsh=dend_thrsh_pav, cbar_lim=cbar_lim,
    #         save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20")
    #     _, cnv_mwu, cnv_shap_v_fitness = cluster_shap_indiv(f"{cnv_files[i]}",
    #         target_envs, pheno, iso_colors, "cnv", plot=True, save_out=True,
    #         dend_thrsh=dend_thrsh_cnv, cbar_lim=cbar_lim,
    #         save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20")
    #     mwu_res[i] = {"snp":snp_mwu, "pav":pav_mwu, "cnv":cnv_mwu}
    #     shap_v_fitness[i] = {"snp": snp_shap_v_fitness, "pav":pav_shap_v_fitness, "cnv":cnv_shap_v_fitness}

    '''Min and max of the median absolute SHAP values of the top 20 features
    snp YPDBENOMYL500 0.0002354298444995 0.0007214783908874
    pav YPDCAFEIN40  0.0019436788700108 0.0058764563217105
    cnv YPDCAFEIN50 0.0013561125454816 0.0040035949655501
    snp YPDCAFEIN40 0.0003227890148337 0.00056593105846
    pav YPDCUSO410MM 0.0059695306881159 0.0200749856723022
    cnv YPDBENOMYL500 0.0011791918268477 0.0046210646092313
    snp YPDCAFEIN50 0.0005349312344456 0.0015409023144522
    pav YPDSODIUMMETAARSENITE 0.0018329834955045 0.0099934814158312
    cnv YPDSODIUMMETAARSENITE 0.0022798527122261 0.017888163432244
    snp YPDSODIUMMETAARSENITE 0.0017063827902945 0.0049816312935072
    pav YPDBENOMYL500 0.0013243758280785 0.0046946133274949
    cnv YPDCAFEIN40 0.0017196560204147 0.0051414742889729
    snp YPDCUSO410MM 0.0030447193434712 0.0113774728699881
    pav YPDCAFEIN50 0.0012298897373209 0.0040465473392858
    cnv YPDCUSO410MM 0.0023580262972865 0.3308934626544931'''

    # linreg_save = "Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_cluster_median_fitness_linreg_RF_FS_models_v4.csv"
    # mwu_save = "Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_cluster_fitness_mwu_RF_FS_models_v4.csv"
    # shap_v_fitness_save = "Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_cluster_values_vs_fitness_linreg_RF_FS_models_v4.csv"
    # save_results(mwu_res, shap_v_fitness, mwu_save, shap_v_fitness_save)

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
            df = pd.read_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v4.tsv", sep="\t", index_col=0)
            df.T.corr(method="pearson").to_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v4_pcc_r_values.tsv", sep="\t")
            df.T.corr(method=lambda x, y: pearsonr(x, y)[1]).to_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v4_pcc_p_values.tsv", sep="\t")
            
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
                                    palette="RdBu_r", edgecolor="black", alpha=0.3, legend=False)
                    ax[i//2][i%2].set_ylabel("PAV value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=pav.loc[df.columns, df.index[i]].astype(int))
                    res.append([env, data_type, df.index[i], r, p])
                if data_type == "cnv":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=cnv.loc[df.columns, df.index[i]].astype(float),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor="black", alpha=0.3, legend=False)
                    ax[i//2][i%2].set_ylabel("CNV value")
                    #
                    r, p = pearsonr(x=pheno.loc[df.columns, env],
                                y=cnv.loc[df.columns, df.index[i]].astype(float))
                    res.append([env, data_type, df.index[i], r, p])
                if data_type == "snp":
                    sns.scatterplot(x=pheno.loc[df.columns, env],
                                    y=snp.loc[df.columns, df.index[i]].astype(int),
                                    hue=df.iloc[i,:], ax=ax[i//2][i%2], hue_norm=norm,
                                    palette="RdBu_r", edgecolor="black", alpha=0.3, legend=False)
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
            plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v4_indiv_feat_trends.pdf",
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
                                palette="GnBu", edgecolor="black", alpha=0.3, hue_norm=norm)
                                # size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]}")
                            ax[i][j].set_xlabel(f"{df.index[i]}")
                        if data_type == "cnv":
                            size_column = cnv.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="GnBu", edgecolor="black", alpha=0.3, hue_norm=norm),
                                # size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]}")
                            ax[i][j].set_xlabel(f"{df.index[i]}")
                        if data_type == "snp":
                            size_column = snp.loc[df.columns, df.index[i]].astype(int)
                            scatter = sns.scatterplot(x=df.iloc[i, :], y=df.iloc[j, :],
                                hue=pheno.loc[df.columns, env], ax=ax[i][j],
                                palette="GnBu", edgecolor="black", alpha=0.3, hue_norm=norm)
                                # size=size_column)
                            ax[i][j].set_ylabel(f"{df.index[j]} SHAP")
                            ax[i][j].set_xlabel(f"{df.index[i]} SHAP")
                        
                        ax[i][j].set_title(f"{df.index[i]} vs {df.index[j]}",
                                            fontsize=8)
                        ax[i][j].axvline(x=0, color="black", linestyle="--")
                        ax[i][j].axhline(y=0, color="black", linestyle="--")
                        
                        # # Create separate legends for hue y size
                        # handles, labels = scatter.get_legend_handles_labels()
                        # size_legend = ax[i][j].legend(
                        #     handles[-len(size_column.unique()):],
                        #     labels[-len(size_column.unique()):],
                        #     bbox_to_anchor=(1.05, 0.5), loc='center left',
                        #     title="Feature value", prop={'size': 6})
                        # hue_legend = ax[i][j].legend(
                        #     handles[1:-len(size_column.unique())-1],
                        #     labels[1:-len(size_column.unique())-1],
                        #     bbox_to_anchor=(1.05, 1), loc='upper left',
                        #     title="Fitness", prop={'size': 6})
                        # ax[i][j].add_artist(size_legend)
            
            for axis in ax.flat:
                axis.set_box_aspect(1)
            
            plt.tight_layout()
            plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_{env}_top20_{data_type}_sorted_v4_feat_pair_trends.pdf",
                        bbox_inches="tight", dpi=300)
            plt.close()


    res_df = pd.DataFrame(res, columns=["Environment", "Data_Type", "Feature", "R", "P-value"])
    res_df.to_csv(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_FS_top20_indiv_feat_trends_pearsonr_v4.tsv")
    
    # # Plot YBR008C (FLR1) SNP value vs Fitness, colored by SHAP (chromosome2_252680)
    # flr1_snps = map_snps.loc[map_snps.gene == 'YBR008C',:].snp.to_list()
    # df = pd.read_csv("/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/SNP/fs/SHAP_values_sorted_YPDBENOMYL500_snp_rf_fs_6000_top_6000_training.txt", sep="\t", index_col=0).T # SHAP values dataframe
    # df.loc[df.index.isin(flr1_snps),:] # only one snp (chromosome2_252680)
    # df_info = pd.DataFrame(df.index).merge(
    #                 map_snps.set_index("snp").drop(columns=["chr", "pos"]),
    #                 left_on=0, right_index=True, how="left")
    
    # vmin = df.loc["chromosome2_252680", :].min()
    # vmax = df.loc["chromosome2_252680", :].max()
    # norm = plt.Normalize(vmin=vmin, vmax=vmax) # center the SHAP values legend at 0
    # sns.scatterplot(x=pheno.loc[df.columns, "YPDBENOMYL500"],
    #                 y=snp.loc[df.columns, "chromosome2_252680"].astype(int),
    #                 hue=df.loc["chromosome2_252680",:], hue_norm=norm,
    #                 palette="RdBu_r", edgecolor=None, alpha=0.7, legend=False)
    # plt.title(f"Feature {'//'.join(df_info.loc[df_info[0]=='chromosome2_252680',:].astype(str).values.flatten())}",
    #           fontsize=8)
    # plt.xlabel("Fitness")
    # sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
    # sm.set_array([])
    # plt.colorbar(sm, ax=plt.gca(), orientation="vertical", label="SHAP value")

    # plt.gca().set_box_aspect(1)

    # plt.tight_layout()
    # plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_YPDBENOMYL500_flr1_snp_sorted_v2_indiv_feat_trend.pdf",
    #             bbox_inches="tight", dpi=300)
    # plt.close()
    
    # # Plot YBR008C (FLR1) vs YBR180W (DTR1, chromosome2_591388)
    # vmax = pheno.loc[df.columns, "YPDBENOMYL500"].max()
    # norm = plt.Normalize(vmin=0, vmax=vmax)
    # scatter = sns.scatterplot(x=df.loc["chromosome2_252680",:],
    #                           y=df.loc["chromosome2_591388",:],
    #                           hue=pheno.loc[df.columns, "YPDBENOMYL500"],
    #                           palette="viridis", edgecolor=None, alpha=0.7,
    #                           hue_norm=norm)
    # plt.ylabel("chromosome2_252680 (FLR1) SHAP")
    # plt.xlabel("chromosome2_591388 (DTR1) SHAP")
    # plt.title(f"chromosome2_252680 vs chromosome2_591388", fontsize=8)
    # plt.gca().axvline(x=0, color="black", linestyle="--")
    # plt.gca().axhline(y=0, color="black", linestyle="--")
    # plt.legend(title="Fitness")
    # plt.gca().set_box_aspect(1)
    # plt.tight_layout()
    # plt.savefig(f"Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_clustered_YPDBENOMYL500_flr1_v_dtr1_snp_sorted_v2_feat_pair_trends.pdf",
    #             bbox_inches="tight", dpi=300)
    # plt.close()
    
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

