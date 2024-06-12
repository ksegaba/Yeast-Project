#!/usr/bin/env python3
################################################################################
# Figure S1
################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

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