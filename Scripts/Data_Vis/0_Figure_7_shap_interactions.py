#! /usr/bin/env python3
"""Analysis for the section on genetic interactions underlying benomyl stress in
yeast """

import joblib, os
import datatable as dt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score
from statannotations.Annotator import Annotator
from itertools import combinations
from scipy.stats import mannwhitneyu

#### Combine the optimized RF SNP, PAV, and CNV features into a single dataframe
# SNP, PAV, CNV feature tables
data = "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018"
snp = dt.fread(f"{data}/geno.csv").to_pandas().set_index("ID")
pav = dt.fread(f"{data}/ORFs_pres_abs.csv").to_pandas().set_index("ID")
cnv = dt.fread(f"{data}/ORFs_no_NA.csv").to_pandas().set_index("ID")

# list of best model 
save_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/"
snp_mod_list = pd.read_csv(f"{save_path}/SNP_yeast_RF_results/fs/snp_best_fs_models.txt", header=None)
pav_mod_list = pd.read_csv(f"{save_path}/ORF_yeast_RF_results/fs/pav_best_fs_models.txt", header=None)
cnv_mod_list = pd.read_csv(f"{save_path}/ORF_yeast_RF_results/fs/cnv_best_fs_models.txt", header=None)

for env in snp_mod_list[0].str.split("_").apply(lambda x: x[0]).unique():
	snp_mod_name = snp_mod_list.loc[snp_mod_list[0].str.contains(env, regex=False), 0].values[0]
	pav_mod_name = pav_mod_list.loc[pav_mod_list[0].str.contains(env, regex=False), 0].values[0]
	cnv_mod_name = cnv_mod_list.loc[cnv_mod_list[0].str.contains(env, regex=False), 0].values[0]
	
	# load the models
	snp_mod = joblib.load(f"{save_path}/SNP_yeast_RF_results/fs/{snp_mod_name}")
	pav_mod = joblib.load(f"{save_path}/ORF_yeast_RF_results/fs/{pav_mod_name}")
	cnv_mod = joblib.load(f"{save_path}/ORF_yeast_RF_results/fs/{cnv_mod_name}")
	
	# get the feature names
	snp_feat = snp_mod.feature_names_in_
	pav_feat = pav_mod.feature_names_in_
	cnv_feat = cnv_mod.feature_names_in_
	
	# integrate the features into a single dataframe
	integrated = pd.concat([snp[snp_feat], pav[pav_feat], cnv[cnv_feat]], axis=1, ignore_index=False)
	
	# save the integrated feature tables to files
	integrated.to_csv(
		f"{save_path}/SHAP_Interaction/optimized_single_env_RF/integrated_models/Integrated_all_variants_{env}.csv",
		index=True) # snp + pav + cnv features
	pd.concat([snp[snp_feat], pav[pav_feat]], axis=1, ignore_index=False).to_csv(
		f"{save_path}/SHAP_Interaction/optimized_single_env_RF/integrated_models/Integrated_snp_pav_{env}.csv",
		index=True) # snp + pav features
	pd.concat([snp[snp_feat], cnv[pav_feat]], axis=1, ignore_index=False).to_csv(
		f"{save_path}/SHAP_Interaction/optimized_single_env_RF/integrated_models/Integrated_snp_cnv_{env}.csv",
		index=True) # snp + cnv features
	pd.concat([pav[pav_feat], cnv[cnv_feat]], axis=1, ignore_index=False).to_csv(
		f"{save_path}/SHAP_Interaction/optimized_single_env_RF/integrated_models/Integrated_pav_cnv_{env}.csv",
		index=True) # pav + cnv features
	
	del integrated, snp_mod, pav_mod, cnv_mod, snp_feat, pav_feat, cnv_feat


#### Combine the benomyl benchmark gene SNP, PAV, and CNV features into a single dataframe
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"

# features where all variants are kept per gene
# these were generated in 0_Table_S10-11_make_benchmark_gene_models.py
snp_feat = pd.read_csv(f"{path}/Features_all_variants_per_gene_benomyl_500ugml_snp.txt", sep="\t", header=None)
pav_feat = pd.read_csv(f"{path}/Features_all_variants_per_gene_benomyl_500ugml_pav.txt", sep="\t", header=None)
cnv_feat = pd.read_csv(f"{path}/Features_all_variants_per_gene_benomyl_500ugml_cnv.txt", sep="\t", header=None)

# features where only the best variant is kept per gene
snp_1feat = pd.read_csv(f"{path}/Features_max_gini_from_RF_baseline_imp_snp_YPDBENOMYL500.txt", sep="\t", header=None)
pav_1feat = pd.read_csv(f"{path}/Features_max_gini_from_RF_baseline_imp_pav_YPDBENOMYL500.txt", sep="\t", header=None)
cnv_1feat = pd.read_csv(f"{path}/Features_max_gini_from_RF_baseline_imp_cnv_YPDBENOMYL500.txt", sep="\t", header=None)

# feature tables and feature to gene maps
data = "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018"
snp = dt.fread(f"{data}/geno.csv").to_pandas().set_index("ID")
pav = dt.fread(f"{data}/ORFs_pres_abs.csv").to_pandas().set_index("ID")
cnv = dt.fread(f"{data}/ORFs_no_NA.csv").to_pandas().set_index("ID")

map_snps = pd.read_csv(f"{data}/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t", index_col=0) # needed for *_feat1 files to extract benomyl genes
map_orfs = pd.read_csv(f"{data}/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t", index_col=0)
map_orfs.index = map_orfs.index.str.replace("-", ".", regex=False) # replace . with - in the index names
map_orfs.index = map_orfs.apply(lambda x: "X" + x.name, axis=1) # remove X from the beginning of the index names

# get one variant per gene
snp_gini = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_4/RF_baseline_imp_snp.tsv").to_pandas()
pav_gini = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_4/RF_baseline_imp_pav.tsv").to_pandas()
cnv_gini = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/Section_4/RF_baseline_imp_cnv.tsv").to_pandas()

snp_gini.set_index(["snp", "gene"], inplace=True)
pav_gini.set_index(["orf", "gene"], inplace=True)
cnv_gini.set_index(["orf", "gene"], inplace=True)

pav_gini.index = pav_gini.index.set_levels(pav_gini.index.levels[0].str.replace("-", ".", regex=False), level=0) # replace . with - in the index names
cnv_gini.index = cnv_gini.index.set_levels(cnv_gini.index.levels[0].str.replace("-", ".", regex=False), level=0)
pav_gini.index = pav_gini.index.set_levels("X" + pav_gini.index.levels[0].astype(str), level=0) # add X to the beginning of the index names
cnv_gini.index = cnv_gini.index.set_levels("X" + cnv_gini.index.levels[0].astype(str), level=0) # add X to the beginning of the index names

snp_gini = snp_gini.loc[:,"YPDBENOMYL500"].sort_values(ascending=False)
pav_gini = pav_gini.loc[:,"YPDBENOMYL500"].sort_values(ascending=False)
cnv_gini = cnv_gini.loc[:,"YPDBENOMYL500"].sort_values(ascending=False)
# snp_gini = snp_gini.loc[snp_gini != 0.0] # remove features with 0 gini importance
# pav_gini = pav_gini.loc[pav_gini != 0.0]
# cnv_gini = cnv_gini.loc[cnv_gini != 0.0]

''' The following are duplicated:
tmp = pav_gini.loc[pav_gini.index.get_level_values(0).isin(ben_pav_1feat),:].reset_index('gene')
tmp.loc[tmp.gene.duplicated(keep=False)]
	X793.augustus_masked.CPI_4.15676  YHL007C   9.048511e-07
	X4200.YHL007C                     YHL007C   4.738338e-07
	X572.snap_masked.12142.CPK_4      YDL225W   2.702119e-07
	X2670.YDL225W                     YDL225W   0.000000e+00
'''

ben_snp_1feat = list(map_snps.loc[snp_gini.index.levels[0].values].loc[map_snps.Benomyl==1].index)
ben_pav_1feat = list(map_orfs.loc[map_orfs.index.isin(pav_gini.index.levels[0].values),:].loc[map_orfs.Benomyl==1].index)
ben_cnv_1feat = list(map_orfs.loc[map_orfs.index.isin(cnv_gini.index.levels[0].values),:].loc[map_orfs.Benomyl==1].index)

final_ben_snp_1feat = snp_gini.loc[snp_gini.index.get_level_values(0).isin(ben_snp_1feat),:].groupby("gene").idxmax() # variants with the largest gini importance
final_ben_pav_1feat = pav_gini.loc[pav_gini.index.get_level_values(0).isin(ben_pav_1feat),:].groupby("gene").idxmax()
final_ben_cnv_1feat = cnv_gini.loc[cnv_gini.index.get_level_values(0).isin(ben_cnv_1feat),:].groupby("gene").idxmax()

pd.Series(pd.DataFrame(final_ben_snp_1feat).apply(lambda x: x[0][0], axis=1).values).to_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp.txt", sep="\t", header=None, index=False)
pd.Series(pd.DataFrame(final_ben_pav_1feat).apply(lambda x: x[0][0], axis=1).values).to_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_pav.txt", sep="\t", header=None, index=False)
pd.Series(pd.DataFrame(final_ben_cnv_1feat).apply(lambda x: x[0][0], axis=1).values).to_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_cnv.txt", sep="\t", header=None, index=False)

# these are the unique PAV and CNV ORFs (this info is important later for determining shap interaction type: pav-cnv, snp-pav-cnv)
final_ben_pav_1feat.reset_index()[pd.Series(pd.DataFrame(final_ben_pav_1feat).apply(lambda x: x[0][0], axis=1).values)!=pd.Series(pd.DataFrame(final_ben_cnv_1feat).apply(lambda x: x[0][0], axis=1).values)]
#       gene                                YPDBENOMYL500
# 77   YDL225W      (X572.snap_masked.12142.CPK_4, YDL225W) # unique to PAVs
# 157  YHL007C  (X793.augustus_masked.CPI_4.15676, YHL007C) # unique to PAVs
final_ben_cnv_1feat.reset_index()[pd.Series(pd.DataFrame(final_ben_cnv_1feat).apply(lambda x: x[0][0], axis=1).values)!=pd.Series(pd.DataFrame(final_ben_pav_1feat).apply(lambda x: x[0][0], axis=1).values)]
#         gene             YPDBENOMYL500
# 77   YDL225W  (X2670.YDL225W, YDL225W) # unique to CNVs
# 157  YHL007C  (X4200.YHL007C, YHL007C) # unique to CNVs

# integrate the features into single dataframes
integrated = pd.concat([snp[list(snp_feat.T.values[0])], pav[list(pav_feat.T.values[0])],
           cnv[list(cnv_feat.T.values[0])]], axis=1, ignore_index=False)

integrated_1feat = pd.concat([
	snp[list(pd.DataFrame(final_ben_snp_1feat).apply(lambda x: x[0][0], axis=1).values)],
	pav[list(pd.DataFrame(final_ben_pav_1feat).apply(lambda x: x[0][0], axis=1).values)],
	cnv[list(pd.DataFrame(final_ben_cnv_1feat).apply(lambda x: x[0][0], axis=1).values)]],
	axis=1, ignore_index=False)

# save the integrated feature tables to files
integrated.to_csv(f"{path}/Features_all_variants_per_gene_benomyl_500ugml_snp_pav_cnv.txt",
    sep=",", index=True) # snp + pav + cnv features
pd.concat([snp[list(snp_feat.T.values[0])], pav[list(pav_feat.T.values[0])]],
    axis=1, ignore_index=False).to_csv(
    f"{path}/Features_all_variants_per_gene_benomyl_500ugml_snp_pav.txt",
    sep=",", index=True) # snp + pav features
pd.concat([snp[list(snp_feat.T.values[0])], cnv[list(cnv_feat.T.values[0])]],
    axis=1, ignore_index=False).to_csv(
    f"{path}/Features_all_variants_per_gene_benomyl_500ugml_snp_cnv.txt",
    sep=",", index=True) # snp + cnv features
pd.concat([pav[list(pav_feat.T.values[0])], cnv[list(cnv_feat.T.values[0])]],
    axis=1, ignore_index=False).to_csv(
    f"{path}/Features_all_variants_per_gene_benomyl_500ugml_pav_cnv.txt",
    sep=",", index=True) # pav + cnv features

integrated_1feat.to_csv(f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp_pav_cnv.txt",
                        index=True) # snp + pav + cnv features
pd.concat([
	snp[list(pd.DataFrame(final_ben_snp_1feat).apply(lambda x: x[0][0], axis=1).values)],
	pav[list(pd.DataFrame(final_ben_pav_1feat).apply(lambda x: x[0][0], axis=1).values)]],
	axis=1, ignore_index=False).to_csv(
	f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp_pav.txt", index=True) # snp + pav features
pd.concat([
	snp[list(pd.DataFrame(final_ben_snp_1feat).apply(lambda x: x[0][0], axis=1).values)],
	cnv[list(pd.DataFrame(final_ben_cnv_1feat).apply(lambda x: x[0][0], axis=1).values)]],
	axis=1, ignore_index=False).to_csv(
	f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp_cnv.txt", index=True) # snp + cnv features
pd.concat([
	pav[list(pd.DataFrame(final_ben_pav_1feat).apply(lambda x: x[0][0], axis=1).values)],
	cnv[list(pd.DataFrame(final_ben_cnv_1feat).apply(lambda x: x[0][0], axis=1).values)]],
	axis=1, ignore_index=False).to_csv(
f"{path}/Features_one_variant_per_gene_benomyl_500ugml_pav_cnv.txt", index=True) # pav + cnv features

##### Get the best benomyl model to use to estimate shap interaction scores ####
test = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", sep="\t", header=None)
snp_scores = pd.read_csv(
    "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/snp_one_variant_per_gene_benomyl_500ugml_scores.txt",
    sep="\t", index_col=0)
pav_scores = pd.read_csv(
    "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/pav_one_variant_per_gene_benomyl_500ugml_scores.txt",
	sep="\t", index_col=0)
cnv_scores = pd.read_csv(
    "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/cnv_one_variant_per_gene_benomyl_500ugml_scores.txt",
	sep="\t", index_col=0)
snp_pav_cnv_scores = pd.read_csv(
	"/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/integrated_one_variant_per_gene_benomyl_500ugml_scores.txt",
	sep="\t", index_col=0)
snp_pav_scores = pd.read_csv(
	"/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/integrated_snp_pav_one_variant_per_gene_benomyl_500ugml_scores.txt",
	sep="\t", index_col=0)
snp_cnv_scores = pd.read_csv(
	"/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/integrated_snp_cnv_one_variant_per_gene_benomyl_500ugml_scores.txt",
	sep="\t", index_col=0)
pav_cnv_scores = pd.read_csv(
	"/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/integrated_pav_cnv_one_variant_per_gene_benomyl_500ugml_scores.txt",
	sep="\t", index_col=0)

save_names = ["snp_one_variant_per_gene_benomyl_500ugml",
			  "pav_one_variant_per_gene_benomyl_500ugml",
			  "cnv_one_variant_per_gene_benomyl_500ugml",
			  "integrated_one_variant_per_gene_benomyl_500ugml",
			  "integrated_snp_pav_one_variant_per_gene_benomyl_500ugml",
			  "integrated_snp_cnv_one_variant_per_gene_benomyl_500ugml",
			  "integrated_pav_cnv_one_variant_per_gene_benomyl_500ugml"]

with open("/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models/best_benomyl_models.txt", "w") as f:
	for i,scores in enumerate([snp_scores, pav_scores, cnv_scores,
		snp_pav_cnv_scores, snp_pav_scores, snp_cnv_scores, pav_cnv_scores]):
		train_scores = scores.loc[~scores.index.isin(test[0]), :]
		best_rep = train_scores.drop(columns=["Mean", "stdev"]).apply(
			lambda c: r2_score(train_scores["Y"], train_scores[c.name] \
			if c.name != "Y" else c), axis=0).iloc[1:].idxmax()
		f.write(f"{save_names[i]}_models_rep_{int(best_rep.split('_')[1]) - 1}.pkl\n")

### Plot A. Model performances ------------------------------------------------#
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models"
rf = pd.read_csv(f"{d}/RESULTS_reg.txt", sep="\t")
rrblup = pd.read_csv(f"{d}/rrBLUP_benomyl_models/RESULTS_rrblup.txt", sep="\t")

rf = rf.loc[rf["DateTime"].str.contains("2025-05-06")] # relevant models
rrblup = rrblup.iloc[7:10,:]

# for plotting purposes, set x-axis labels
rf["Model"] = ["SNP", "PAV", "CNV", "SNP + PAV + CNV", "SNP + PAV", "SNP + CNV", "PAV + CNV"]
rrblup["Model"] = ["SNP", "PAV", "CNV"]

fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(5, 3.5))
sns.barplot(data=rf, x="Model", y="r2_test", ax=ax[0])
ax[0].errorbar(data=rf, x="Model", y="r2_test", yerr="r2_test_sd", fmt="o", color="black", capsize=5)
ax[0].set_ylim(0, 1) ; ax[0].set_ylabel("RF R2 test")
ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=45, ha="right")
sns.barplot(data=rrblup, x="Model", y="R2_test", ax=ax[1])
ax[1].errorbar(data=rrblup, x="Model", y="R2_test", yerr="R2_test_sd", fmt="o", color="black", capsize=5)
ax[1].set_ylabel("rrBLUP R2 test")
# for axis in ax.flat:
#     axis.set_aspect('equal', adjustable='box')  # Makes the axis square

plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7A_benomyl_model_performances_RF_rrBLUP.pdf")
plt.close()
### ---------------------------------------------------------------------------#

################ Exploratory analysis of SHAP interaction scores ###############
os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

def get_unique_gp(df, Gene1="Gene1", Gene2="Gene2"):
	'''
	Get the unique gene pairs from the dataframe
	'''
	df_gp = df.apply(lambda x: set([x[Gene1], x[Gene2]]), axis=1).values # gene pairs
	df_gp = {frozenset(sorted(set)) for set in df_gp} # get unique interactions
	return df_gp


# SNP and ORF gene maps
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")
map_snps_dict = map_snps[["snp", "gene"]].set_index("snp").to_dict(orient="dict")
map_orfs_dict = map_orfs[["orf", "gene"]].set_index("orf").to_dict(orient="dict")

#### SHAP interaction scores from SNP, PAV, CNV, and Integrated models
path="/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/benomyl_models"
snp_shap_int = pd.read_csv(f"{path}/SNP/shap_interaction_scores_snp_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
pav_shap_int = pd.read_csv(f"{path}/PAV/shap_interaction_scores_pav_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
cnv_shap_int = pd.read_csv(f"{path}/CNV/shap_interaction_scores_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
spc_shap_int = pd.read_csv(f"{path}/INTEGRATED/snp_pav_cnv/shap_interaction_scores_integrated_snp_pav_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
sp_shap_int = pd.read_csv(f"{path}/INTEGRATED/snp_pav/shap_interaction_scores_integrated_snp_pav_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
sc_shap_int = pd.read_csv(f"{path}/INTEGRATED/snp_cnv/shap_interaction_scores_integrated_snp_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
pc_shap_int = pd.read_csv(f"{path}/INTEGRATED/pav_cnv/shap_interaction_scores_integrated_pav_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")

### Plot ?. SHAP interaction score violin plots -------------------------------#
snp_shap_int.insert(3, "Model", "SNP")
pav_shap_int.insert(3, "Model", "PAV")
cnv_shap_int.insert(3, "Model", "CNV")
spc_shap_int.insert(3, "Model", "SNP + PAV + CNV")
sp_shap_int.insert(3, "Model", "SNP + PAV")
sc_shap_int.insert(3, "Model", "SNP + CNV")
pc_shap_int.insert(3, "Model", "PAV + CNV")

violin_data = pd.concat([snp_shap_int, pav_shap_int, cnv_shap_int,
	spc_shap_int, sp_shap_int, sc_shap_int, pc_shap_int], axis=0,
	ignore_index=True)

from statannotations.Annotator import Annotator
from itertools import combinations
from scipy.stats import mannwhitneyu
models = violin_data["Model"].unique()
pairs = list(combinations(models, 2))  # all pairwise combinations

ax = sns.violinplot(data=violin_data, inner="box", fill=False, x="Model", log_scale=10,
	y="Interaction", hue="Model")
# sns.stripplot(data=violin_data, x="Model", y="Interaction", hue="Model", dodge=True,
# 	alpha=0.5, size=1)
plt.axhline(violin_data["Interaction"].quantile(0.95), color="red", linestyle="--")
plt.axhline(violin_data["Interaction"].quantile(0.99), color="blue", linestyle="--")
annotator = Annotator(ax, pairs, data=violin_data,
                      x="Model", y="Interaction", order=models)
annotator.configure(test="Mann-Whitney", text_format="star", loc="outside", comparisons_correction="bonferroni")
annotator.apply_and_annotate()
'''Annotator results:
p-value annotation legend:
      ns: 5.00e-02 < p <= 1.00e+00
       *: 1.00e-02 < p <= 5.00e-02
      **: 1.00e-03 < p <= 1.00e-02
     ***: 1.00e-04 < p <= 1.00e-03
    ****: p <= 1.00e-04

SNP vs. PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:3.527e-05 U_stat=9.511e+06
PAV vs. CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:3.311e-15 U_stat=4.107e+06
CNV vs. SNP + PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.269e-02 U_stat=9.109e+08
SNP + PAV + CNV vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=2.716e+09
SNP + PAV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:3.163e-33 U_stat=1.796e+09
SNP + CNV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=5.451e+08
SNP vs. CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=3.878e+08
PAV vs. SNP + PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:3.543e-18 U_stat=3.488e+07
CNV vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=3.330e+08
SNP + PAV + CNV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=2.075e+09
SNP + PAV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=8.318e+08
SNP vs. SNP + PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=3.286e+09
PAV vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.000e+00 U_stat=1.578e+07
CNV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=2.541e+08
SNP + PAV + CNV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.065e-296 U_stat=1.107e+09
SNP vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=1.501e+09
PAV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.000e+00 U_stat=1.112e+07
CNV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:4.409e-198 U_stat=1.344e+08
SNP vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=1.043e+09
PAV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.911e-41 U_stat=4.817e+06
SNP vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=4.493e+08
'''
plt.ylabel("SHAP interaction score")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7B_violin_log10_summed_shap_interaction.pdf")
plt.close()

rows_above_95th = violin_data[violin_data["Interaction"] > violin_data["Interaction"].quantile(0.95)]
rows_above_99th = violin_data[violin_data["Interaction"] > violin_data["Interaction"].quantile(0.99)]
rows_above_95th["Model"].value_counts()
# Model
# SNP + CNV          4853
# SNP                4496
# SNP + PAV + CNV    3243
# SNP + PAV          1744
# CNV                 649
# PAV + CNV           380
# PAV                  89
rows_above_99th["Model"].value_counts()
# SNP + CNV          1293
# SNP                 963
# SNP + PAV + CNV     395
# CNV                 163
# SNP + PAV           138
# PAV + CNV            96
# PAV                  43

### ---------------------------------------------------------------------------#
"""To determine the number of unique genetic interactions, there are different
levels to consider: variant-to-variant (might have duplicate interactions at the
gene level) and gene-to-gene (this is the one we care about most)."""

# Determine from which dataset Feature1 and Feature2 come from in the integrated data
spc_shap_int.insert(4, "Feature1_Data", spc_shap_int.apply(lambda x: "CNV" if x["Feature1"].endswith(".0") else \
		("PAV" if x["Feature1"].startswith("X") else "SNP"), axis=1))
spc_shap_int.insert(5, "Feature2_Data", spc_shap_int.apply(lambda x: "CNV" if x["Feature2"].endswith(".0") else \
		("PAV" if x["Feature2"].startswith("X") else "SNP"), axis=1))
sp_shap_int.insert(4, "Feature1_Data", sp_shap_int.apply(lambda x: "CNV" if x["Feature1"].endswith(".0") else \
		("PAV" if x["Feature1"].startswith("X") else "SNP"), axis=1))
sp_shap_int.insert(5, "Feature2_Data", sp_shap_int.apply(lambda x: "CNV" if x["Feature2"].endswith(".0") else \
		("PAV" if x["Feature2"].startswith("X") else "SNP"), axis=1))
sc_shap_int.insert(4, "Feature1_Data", sc_shap_int.apply(lambda x: "CNV" if x["Feature1"].endswith(".0") else \
		("PAV" if x["Feature1"].startswith("X") else "SNP"), axis=1))
sc_shap_int.insert(5, "Feature2_Data", sc_shap_int.apply(lambda x: "CNV" if x["Feature2"].endswith(".0") else \
		("PAV" if x["Feature2"].startswith("X") else "SNP"), axis=1))
pc_shap_int.insert(4, "Feature1_Data", pc_shap_int.apply(lambda x: "CNV" if x["Feature1"].endswith(".0") else \
		("PAV" if x["Feature1"].startswith("X") else "SNP"), axis=1))
pc_shap_int.insert(5, "Feature2_Data", pc_shap_int.apply(lambda x: "CNV" if x["Feature2"].endswith(".0") else \
		("PAV" if x["Feature2"].startswith("X") else "SNP"), axis=1))

# Prepare ORF identifiers for mapping to genes
pav_shap_int.Feature1 = pav_shap_int.Feature1.str.replace("^X", "", regex=True) # remove X from the beginning of the column names
pav_shap_int.Feature2 = pav_shap_int.Feature2.str.replace("^X", "", regex=True)
pav_shap_int.Feature1 = pav_shap_int.Feature1.str.replace("\.", "-", regex=True) # replace . with - in the column names
pav_shap_int.Feature2 = pav_shap_int.Feature2.str.replace("\.", "-", regex=True)

cnv_shap_int.Feature1 = cnv_shap_int.Feature1.str.replace("^X", "", regex=True) # remove X from the beginning of the column names
cnv_shap_int.Feature2 = cnv_shap_int.Feature2.str.replace("^X", "", regex=True)
cnv_shap_int.Feature1 = cnv_shap_int.Feature1.str.replace("\.", "-", regex=True) # replace . with - in the column names
cnv_shap_int.Feature2 = cnv_shap_int.Feature2.str.replace("\.", "-", regex=True)

spc_shap_int.Feature1 = spc_shap_int.Feature1.str.replace("\.0$", "", regex=True) # remove the .0 from the end of the column names
sp_shap_int.Feature1 = sp_shap_int.Feature1.str.replace("\.0$", "", regex=True)
sc_shap_int.Feature1 = sc_shap_int.Feature1.str.replace("\.0$", "", regex=True)
pc_shap_int.Feature1 = pc_shap_int.Feature1.str.replace("\.0$", "", regex=True)
spc_shap_int.Feature2 = spc_shap_int.Feature2.str.replace("\.0$", "", regex=True)
sp_shap_int.Feature2 = sp_shap_int.Feature2.str.replace("\.0$", "", regex=True)
sc_shap_int.Feature2 = sc_shap_int.Feature2.str.replace("\.0$", "", regex=True)
pc_shap_int.Feature2 = pc_shap_int.Feature2.str.replace("\.0$", "", regex=True)

spc_shap_int.Feature1 = spc_shap_int.Feature1.str.replace("^X", "", regex=True) # remove X from the beginning of the column names
sp_shap_int.Feature1 = sp_shap_int.Feature1.str.replace("^X", "", regex=True)
sc_shap_int.Feature1 = sc_shap_int.Feature1.str.replace("^X", "", regex=True)
pc_shap_int.Feature1 = pc_shap_int.Feature1.str.replace("^X", "", regex=True)
spc_shap_int.Feature2 = spc_shap_int.Feature2.str.replace("^X", "", regex=True)
sp_shap_int.Feature2 = sp_shap_int.Feature2.str.replace("^X", "", regex=True)
sc_shap_int.Feature2 = sc_shap_int.Feature2.str.replace("^X", "", regex=True)
pc_shap_int.Feature2 = pc_shap_int.Feature2.str.replace("^X", "", regex=True)

spc_shap_int.Feature1 = spc_shap_int.Feature1.str.replace(".", "-", regex=False) # replace . with - in the column names
sp_shap_int.Feature1 = sp_shap_int.Feature1.str.replace(".", "-", regex=False)
sc_shap_int.Feature1 = sc_shap_int.Feature1.str.replace(".", "-", regex=False)
pc_shap_int.Feature1 = pc_shap_int.Feature1.str.replace(".", "-", regex=False)
spc_shap_int.Feature2 = spc_shap_int.Feature2.str.replace(".", "-", regex=False)
sp_shap_int.Feature2 = sp_shap_int.Feature2.str.replace(".", "-", regex=False)
sc_shap_int.Feature2 = sc_shap_int.Feature2.str.replace(".", "-", regex=False)
pc_shap_int.Feature2 = pc_shap_int.Feature2.str.replace(".", "-", regex=False)

# map the SNPs, PAVs, and CNVs to genes
def feature2gene(feature):
	try:
		return map_orfs_dict["gene"][feature]
	except KeyError:
		try:
			return map_snps_dict["gene"][feature]
		except KeyError:
			if feature.endswith("NumOfGenes_3"):
				actual_feat = feature.replace("NumOfGenes_3", "NumOfGenes_2")
				return map_orfs_dict["gene"][actual_feat]


snp_shap_int.insert(0, "Gene2", snp_shap_int.apply(lambda x: map_snps_dict["gene"][x["Feature2"]], axis=1))
snp_shap_int.insert(0, "Gene1", snp_shap_int.apply(lambda x: map_snps_dict["gene"][x["Feature1"]], axis=1))
pav_shap_int.insert(0, "Gene2", pav_shap_int.apply(lambda x: map_orfs_dict["gene"][x["Feature2"]], axis=1))
pav_shap_int.insert(0, "Gene1", pav_shap_int.apply(lambda x: map_orfs_dict["gene"][x["Feature1"]], axis=1))
cnv_shap_int.insert(0, "Gene2", cnv_shap_int.apply(lambda x: map_orfs_dict["gene"][x["Feature2"]], axis=1))
cnv_shap_int.insert(0, "Gene1", cnv_shap_int.apply(lambda x: map_orfs_dict["gene"][x["Feature1"]], axis=1))
spc_shap_int.insert(0, "Gene2", spc_shap_int.apply(lambda x: feature2gene(x["Feature2"]), axis=1))
spc_shap_int.insert(0, "Gene1", spc_shap_int.apply(lambda x: feature2gene(x["Feature1"]), axis=1))
sp_shap_int.insert(0, "Gene2", sp_shap_int.apply(lambda x: feature2gene(x["Feature2"]), axis=1))
sp_shap_int.insert(0, "Gene1", sp_shap_int.apply(lambda x: feature2gene(x["Feature1"]), axis=1))
sc_shap_int.insert(0, "Gene2", sc_shap_int.apply(lambda x: feature2gene(x["Feature2"]), axis=1))
sc_shap_int.insert(0, "Gene1", sc_shap_int.apply(lambda x: feature2gene(x["Feature1"]), axis=1))
pc_shap_int.insert(0, "Gene2", pc_shap_int.apply(lambda x: feature2gene(x["Feature2"]), axis=1))
pc_shap_int.insert(0, "Gene1", pc_shap_int.apply(lambda x: feature2gene(x["Feature1"]), axis=1))

# check for missing gene identifiers
for df in [snp_shap_int, pav_shap_int, cnv_shap_int, spc_shap_int, sp_shap_int, sc_shap_int, pc_shap_int]:
    print(df.loc[spc_shap_int.isna().any(axis=1),:].shape) # nothing was missing, no problems with the gene names.

# pc_shap_int.loc[pc_shap_int.Gene2.isna(), "Feature2_Data"].unique() # it's a PAV problem for all the below
# pc_miss_feat2 = pc_shap_int.loc[pc_shap_int.Gene2.isna(), "Feature2"]
# pc_miss_feat1 = pc_shap_int.loc[pc_shap_int.Gene1.isna(), "Feature1"]
# spc_miss_feat2 = spc_shap_int.loc[spc_shap_int.Gene2.isna(), "Feature2"]
# spc_miss_feat1 = spc_shap_int.loc[spc_shap_int.Gene1.isna(), "Feature1"]
# Note: all the features end in _3, so perhaps RF internally changed it from _2 to _3
# pd.DataFrame(pc_miss_feat2).apply(lambda x: map_orfs.loc[map_orfs.orf.str.contains(x.split("-")[0]), "orf"].values[0], axis=1)
# pd.DataFrame(pc_miss_feat2).apply(lambda x: map_orfs.loc[map_orfs.orf.str.contains(x.split("-")[0]), "orf"].values[0], axis=1)

# get the unique gene pairs
snp_shap_int_gp = get_unique_gp(snp_shap_int, "Gene1", "Gene2")
pav_shap_int_gp = get_unique_gp(pav_shap_int, "Gene1", "Gene2")
cnv_shap_int_gp = get_unique_gp(cnv_shap_int, "Gene1", "Gene2")
len(snp_shap_int_gp) # 36756
len(pav_shap_int_gp) # 458
len(cnv_shap_int_gp) # 14629
len(snp_shap_int_gp.union(pav_shap_int_gp).union(cnv_shap_int_gp)) # 43,862 unique interactions

spc_shap_int_gp = get_unique_gp(spc_shap_int, "Gene1", "Gene2")
sp_shap_int_gp = get_unique_gp(sp_shap_int, "Gene1", "Gene2")
sc_shap_int_gp = get_unique_gp(sc_shap_int, "Gene1", "Gene2")
pc_shap_int_gp = get_unique_gp(pc_shap_int, "Gene1", "Gene2")
len(spc_shap_int_gp) # 58862
len(sp_shap_int_gp) # 64537
len(sc_shap_int_gp) # 36266
len(pc_shap_int_gp) # 14460
len(spc_shap_int_gp.union(sp_shap_int_gp).union(sc_shap_int_gp).union(pc_shap_int_gp)) # 68773
len(snp_shap_int_gp.union(pav_shap_int_gp).union(cnv_shap_int_gp).\
	union(spc_shap_int_gp).union(sp_shap_int_gp).union(sc_shap_int_gp).union(pc_shap_int_gp)) # 69523 unique gene-gene interactions

tableS13_data = pd.concat([snp_shap_int, pav_shap_int, cnv_shap_int,
    spc_shap_int, sp_shap_int, sc_shap_int, pc_shap_int], axis=0,
    ignore_index=True)

len(get_unique_gp(tableS13_data, "Gene1", "Gene2")) # sanity check: 69523

# Fix the CNV feature names whose suffix was changed (due to concatenating PAV & CNV columns with the same names)
# tableS13_data["Feature1"] = tableS13_data.apply(lambda x: x["Feature1"] if x["Feature1"] in map_orfs.orf.values else (x["Feature1"].replace("NumOfGenes_3", "NumOfGenes_2") if x["Feature1"].endswith("NumOfGenes_3") else x["Feature1"]), axis=1)
# tableS13_data["Feature2"] = tableS13_data.apply(lambda x: x["Feature2"] if x["Feature2"] in map_orfs.orf.values else (x["Feature2"].replace("NumOfGenes_3", "NumOfGenes_2") if x["Feature2"].endswith("NumOfGenes_3") else x["Feature2"]), axis=1)

# Save all interactions to a file
tableS13_data.loc[tableS13_data.Model=="SNP", ["Feature1_Data", "Feature2_Data"]] = "SNP"
tableS13_data.loc[tableS13_data.Model=="PAV", ["Feature1_Data", "Feature2_Data"]] = "PAV"
tableS13_data.loc[tableS13_data.Model=="CNV", ["Feature1_Data", "Feature2_Data"]] = "CNV"
tableS13_data.isna().sum() # no missing values!
# tableS13_data.to_csv(
# 	"Scripts/Data_Vis/Section_6/shap_interaction/Table_S13_benomyl_benchmark_RF_models_SHAP_interactions.txt",
# 	sep="\t", index=False)

# the PAV model did not perform well, so we will not consider it
ts13_no_pav = tableS13_data[tableS13_data.Model!="PAV"].copy(deep=True)

# check for duplicates in ts13_no_pav
ts13_no_pav["Pair"] = ts13_no_pav.apply(lambda x: tuple(sorted([x.Feature1, x.Feature2])), axis=1)
ts13_no_pav.sort_values("Interaction", ascending=False, inplace=True)
ts13_no_pav.groupby("Pair").count().value_counts() # some pairs appear 2 to 6 times
ts13_no_pav = ts13_no_pav.drop_duplicates(subset=["Pair"], keep="first") # keep the first one (has the highest interaction score)

'''another sanity check (regarding ORFs on line 133):
ts13_no_pav.loc[(ts13_no_pav.Feature1=='572-snap_masked-12142-CPK_4'), 'Feature1_Data'].unique() # PAV, great!
ts13_no_pav.loc[(ts13_no_pav.Feature2=='572-snap_masked-12142-CPK_4'), 'Feature2_Data'].unique() # PAV :)
ts13_no_pav.loc[(ts13_no_pav.Feature1=='793-augustus_masked-CPI_4-15676'), 'Feature1_Data'].unique() # PAV :)
ts13_no_pav.loc[(ts13_no_pav.Feature2=='793-augustus_masked-CPI_4-15676'), 'Feature2_Data'].unique() # PAV :)
ts13_no_pav.loc[(ts13_no_pav.Feature1=='2670-YDL225W'), 'Feature1_Data'].unique() # CNV & PAV, not good!
ts13_no_pav.loc[(ts13_no_pav.Feature2=='2670-YDL225W'), 'Feature2_Data'].unique() # CNV & PAV :(
ts13_no_pav.loc[(ts13_no_pav.Feature1=='4200-YHL007C'), 'Feature1_Data'].unique() # CNV & PAV :(
ts13_no_pav.loc[(ts13_no_pav.Feature2=='4200-YHL007C'), 'Feature2_Data'].unique() # CNV & PAV :(
'''
ts13_no_pav.loc[(ts13_no_pav.Feature1=='2670-YDL225W'), 'Feature1_Data'] = "CNV"
ts13_no_pav.loc[(ts13_no_pav.Feature2=='2670-YDL225W'), 'Feature2_Data'] = "CNV"
ts13_no_pav.loc[(ts13_no_pav.Feature1=='4200-YHL007C'), 'Feature1_Data'] = "CNV"
ts13_no_pav.loc[(ts13_no_pav.Feature2=='4200-YHL007C'), 'Feature2_Data'] = "CNV"

ts13_no_pav.groupby(["Feature1_Data", "Feature2_Data"]).count() # variant-variant interactions; will be greater than the gene-gene interactions
len(get_unique_gp(ts13_no_pav, "Gene1", "Gene2")) # excluding the PAV model, 69523 unique gene-gene interactions

# Need to collapse the variant-variant interaction types to 6 types:
ts13_no_pav["GI_Type"] = ts13_no_pav.Feature1_Data + "-" + ts13_no_pav.Feature2_Data
ts13_no_pav.GI_Type.replace('CNV-PAV', 'PAV-CNV', inplace=True)
ts13_no_pav.GI_Type.replace('PAV-SNP', 'SNP-PAV', inplace=True)
ts13_no_pav.GI_Type.replace('CNV-SNP', 'SNP-CNV', inplace=True)
ts13_no_pav.GI_Type.nunique() # 6 :)

snp_snp = ts13_no_pav.loc[(ts13_no_pav.Feature1_Data=="SNP") & (ts13_no_pav.Feature2_Data=="SNP")]
len(get_unique_gp(snp_snp, "Gene1", "Gene2")) # 67403 unique gene-gene represented by SNP-SNP
pav_pav = ts13_no_pav.loc[(ts13_no_pav.Feature1_Data=="PAV") & (ts13_no_pav.Feature2_Data=="PAV")]
len(get_unique_gp(pav_pav, "Gene1", "Gene2")) # 18905 PAV-PAV
cnv_cnv = ts13_no_pav.loc[(ts13_no_pav.Feature1_Data=="CNV") & (ts13_no_pav.Feature2_Data=="CNV")]
len(get_unique_gp(cnv_cnv, "Gene1", "Gene2")) # 12537 CNV-CNV
snp_pav = ts13_no_pav.loc[((ts13_no_pav.Feature1_Data=="SNP") & (ts13_no_pav.Feature2_Data=="PAV")) |
			((ts13_no_pav.Feature1_Data=="PAV") & (ts13_no_pav.Feature2_Data=="SNP"))]
len(get_unique_gp(snp_pav, "Gene1", "Gene2")) # 45459 SNP-PAV
snp_cnv = ts13_no_pav.loc[((ts13_no_pav.Feature1_Data=="SNP") & (ts13_no_pav.Feature2_Data=="CNV")) |
			((ts13_no_pav.Feature1_Data=="CNV") & (ts13_no_pav.Feature2_Data=="SNP"))]
len(get_unique_gp(snp_cnv, "Gene1", "Gene2")) # 7688 SNP-CNV
pav_cnv = ts13_no_pav.loc[((ts13_no_pav.Feature1_Data=="PAV") & (ts13_no_pav.Feature2_Data=="CNV")) |
			((ts13_no_pav.Feature1_Data=="CNV") & (ts13_no_pav.Feature2_Data=="PAV"))]
len(get_unique_gp(pav_cnv, "Gene1", "Gene2")) # 4327 PAV-CNV

### Plot B. UpSet diagram of variant-variant shap interactions------------------#
from venn import venn
sets = {"SNP-SNP": get_unique_gp(snp_snp, "Gene1", "Gene2"), 
	"PAV-PAV": get_unique_gp(pav_pav, "Gene1", "Gene2"),
	"CNV-CNV": get_unique_gp(cnv_cnv, "Gene1", "Gene2"),
	"SNP-PAV": get_unique_gp(snp_pav, "Gene1", "Gene2"),
	"SNP-CNV": get_unique_gp(snp_cnv, "Gene1", "Gene2"),
	"PAV-CNV": get_unique_gp(pav_cnv, "Gene1", "Gene2")}
venn(sets)
plt.title("Number of gene-gene interactions identified by SHAP Interactions")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7B_num_unique_GI_in_shap_interactions.pdf")
plt.close()

from upsetplot import UpSet
from upsetplot import from_memberships
from collections import defaultdict
def pair2sets(sets_dict):
	# Invert the dictionary: map each gene pair to all sets it belongs to.
	# A set refers to all the variant-variant type categories.
	pair_to_sets = defaultdict(set)
	for set_name, pairs in sets_dict.items():
		for pair in pairs:
			pair_to_sets[pair].add(set_name)
	# create the memberships list
	memberships = list(pair_to_sets.values())
	# create a count for each membership (all 1s if unique)
	data = from_memberships(memberships, data=[1]*len(memberships))
	return data


UpSet(pair2sets(sets), subset_size="count", show_percentages="{:.1%}").plot()
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7B_num_unique_GI_in_shap_interactions_upset_plot.pdf")
plt.close()
del sets

'''I wanted to add the interaction scores to the UpSet plot, but it didn't work
pair_to_sets = defaultdict(set)
pair_to_scores = {}
for _, row in ts13_no_pav.iterrows():
    pair = frozenset([row["Gene1"], row["Gene2"]])
    pair_to_sets[pair].add(row["GI_Type"])
    # If a pair appears multiple times, you might want to average, max, or sum
    if pair in pair_to_scores:
        pair_to_scores[pair].append(row["Interaction"])
    else:
        pair_to_scores[pair] = [row["Interaction"]]


# Finalize the memberships and interaction scores
memberships = []
scores = []
for pair in pair_to_sets:
    memberships.append(pair_to_sets[pair])
    # For simplicity, use max score per pair; adjust as needed (mean, sum, etc.)
    scores.append(max(pair_to_scores[pair]))


# Build UpSet plot data
upset_data = from_memberships(memberships, data=scores)
# Plot UpSet diagram with interaction score
tmp = pair2sets(sets)
#upset_data = pd.concat([tmp, upset_data], axis=1, ignore_index=False)
upset = UpSet(upset_data, subset_size="sum", show_counts=True, show_percentages="{:.1%}").plot()
# upset.add_catplot(kind="strip", color="black")
# upset.plot()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7B_num_unique_GI_in_shap_interactions_upset_plo    t.pdf")
plt.close()'''
### ---------------------------------------------------------------------------#

# Which variant-variant types tend to have higher interaction scores?
from scipy.stats import mannwhitneyu
variant_pairs = {"snp_snp": snp_snp,"pav_pav": pav_pav, "cnv_cnv": cnv_cnv,
				 "snp_pav": snp_pav, "snp_cnv": snp_cnv, "pav_cnv": pav_cnv}
res = {"Comparison":[], "Statistic":[], "p-value":[]}
for name1, data1 in variant_pairs.items():
	for name2, data2 in variant_pairs.items():
		if name1 != name2:
			stat, p = mannwhitneyu(data1["Interaction"], data2["Interaction"], alternative="greater")
			print(f"{name1} vs {name2}: U={stat:.2f}, p={p:.4e}")
			res["Comparison"].append(f"{name1} > {name2}?")
			res["Statistic"].append(f"{stat:.2f}")
			res["p-value"].append(f"{p:.4e}")
			del stat, p

pd.DataFrame(res).to_csv(
	"Scripts/Data_Vis/Section_6/shap_interaction/Table_S14_benomyl_benchmark_RF_models_SHAP_interactions_mwu_results.txt",
	sep="\t", index=False)

### Plot ?. SHAP interaction score violin plots -------------------------------#
fig, ax = plt.subplots(nrows=2, ncols=1)
sns.violinplot(data=ts13_no_pav, inner="box", fill=False, x="GI_Type",
	log_scale=10, y="Interaction", hue="Model", ax=ax[0])
sns.violinplot(data=ts13_no_pav, inner="box", fill=False, x="GI_Type",
	log_scale=10, y="Interaction", hue="GI_Type", ax=ax[1])
plt.axhline(ts13_no_pav["Interaction"].quantile(0.95), color="red", linestyle="--")
plt.axhline(ts13_no_pav["Interaction"].quantile(0.99), color="blue", linestyle="--")
plt.ylabel("SHAP interaction score")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7B_violin_log10_summed_shap_interaction_variant-variant.pdf")
plt.close()

rows_above_95th = ts13_no_pav[ts13_no_pav["Interaction"] > ts13_no_pav["Interaction"].quantile(0.95)]
rows_above_99th = ts13_no_pav[ts13_no_pav["Interaction"] > ts13_no_pav["Interaction"].quantile(0.99)]
rows_above_95th["Model"].value_counts()
# Model
# SNP + CNV          3117
# SNP                2635
# SNP + PAV + CNV    1482
# SNP + PAV           652
# CNV                 388
# PAV + CNV           243
rows_above_99th["Model"].value_counts()
# Model
# SNP + CNV          787
# SNP                534
# SNP + PAV + CNV    166
# CNV                108
# PAV + CNV           69
# SNP + PAV           40
### ---------------------------------------------------------------------------#

###### Enrichment of known genetic interactions in the SHAP interactions #######
# Remove PAV model interactions since the model was not predictive, these are probably random interactions
ts13_no_pav.shape # (170331, 8); 458 rows from PAV model and additional duplicate feature pairs removed

### BioGRID database genetic interactions
biogrid = dt.fread("~/Shiu_Lab/Project/Data/BioGRID/yeast_gi_biogrid.txt").to_pandas()
biogrid = biogrid.iloc[:,[5,6,7,8,11,13,14,17,36]]
biogrid.columns = ["Systematic Name Interactor A", "Systematic Name Interactor B",
                   "Standard Name Interactor A", "Standard name Interactor B",
                   "Evidence", "Author", "PMID", "Throughput", "Organism"]
biogrid = biogrid.loc[biogrid.Organism=="Saccharomyces cerevisiae (S288c)",:]
biogrid = biogrid.loc[biogrid.Evidence.str.strip().isin(["Synthetic Growth Defect",\
            "Synthetic Lethality", "Synthetic Rescue", "Negative Genetic",\
            "Positive Genetic"]),:] # remove overexpression gene pairs

biogrid_gp = get_unique_gp(biogrid, "Systematic Name Interactor A", "Systematic Name Interactor B")
len(biogrid_gp) # 438546

### Costanzo et al. 2021 benomyl genetic interactions
costanzo = pd.read_excel(
	"~/Shiu_Lab/Project/Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx",\
	engine="openpyxl", sheet_name="Genome-scale_Benomyl")

# filter the genetic interaction network using strict or intermediate criteria
costanzo_ben_strict = costanzo.loc[(costanzo.mean_condition_epsilon.abs() > 0.12) & \
	(costanzo.condition_p_value < 0.05),:] # benomyl interactions with high confidence
costanzo_ben_strict.insert(0, "array_gene", costanzo_ben_strict.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_ben_strict.insert(0, "query_gene", costanzo_ben_strict.apply(lambda x: x.query_orf.split("_")[0], axis=1))

costanzo_ben_intermediate = costanzo.loc[(costanzo.mean_condition_epsilon.abs() > 0.08) & \
	(costanzo.condition_p_value < 0.05),:] # benomyl interactions with intermediate confidence
costanzo_ben_intermediate.insert(0, "array_gene", costanzo_ben_intermediate.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_ben_intermediate.insert(0, "query_gene", costanzo_ben_intermediate.apply(lambda x: x.query_orf.split("_")[0], axis=1))

costanzo_ctrl_strict = costanzo.loc[(costanzo.mean_reference_epsilon.abs() > 0.12) & \
	(costanzo.reference_p_value < 0.05),:] # control condition interactions with high confidence
costanzo_ctrl_strict.insert(0, "array_gene", costanzo_ctrl_strict.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_ctrl_strict.insert(0, "query_gene", costanzo_ctrl_strict.apply(lambda x: x.query_orf.split("_")[0], axis=1))

costanzo_ctrl_intermediate = costanzo.loc[(costanzo.mean_reference_epsilon.abs() > 0.08) & \
	(costanzo.reference_p_value < 0.05),:] # control condition interactions with intermediate confidence
costanzo_ctrl_intermediate.insert(0, "array_gene", costanzo_ctrl_intermediate.apply(lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_ctrl_intermediate.insert(0, "query_gene", costanzo_ctrl_intermediate.apply(lambda x: x.query_orf.split("_")[0], axis=1))

# costanzo_diff_strict = costanzo.loc[(costanzo.mean_differential_epsilon.abs() > 0.12) & \
# 	(costanzo.differential_p_value < 0.05),:] # differential interaction network with high confidence
# costanzo_diff_strict.insert(0, "array_gene", costanzo_diff_strict.apply(lambda x: x.array_orf.split("_")[0], axis=1))
# costanzo_diff_strict.insert(0, "query_gene", costanzo_diff_strict.apply(lambda x: x.query_orf.split("_")[0], axis=1))

# costanzo_diff_intermediate = costanzo.loc[(costanzo.mean_differential_epsilon.abs() > 0.08) & \
# 	(costanzo.differential_p_value < 0.05),:] # differential interaction network with intermediate confidence
# costanzo_diff_intermediate.insert(0, "array_gene", costanzo_diff_intermediate.apply(lambda x: x.array_orf.split("_")[0], axis=1))
# costanzo_diff_intermediate.insert(0, "query_gene", costanzo_diff_intermediate.apply(lambda x: x.query_orf.split("_")[0], axis=1))

# get the unique gene pairs
costanzo_ben_strict_gp = get_unique_gp(costanzo_ben_strict, "query_gene", "array_gene")
costanzo_ben_intermediate_gp = get_unique_gp(costanzo_ben_intermediate, "query_gene", "array_gene")
costanzo_ctrl_strict_gp = get_unique_gp(costanzo_ctrl_strict, "query_gene", "array_gene")
costanzo_ctrl_intermediate_gp = get_unique_gp(costanzo_ctrl_intermediate, "query_gene", "array_gene")
# costanzo_diff_strict_gp = get_unique_gp(costanzo_diff_strict, "query_gene", "array_gene")
# costanzo_diff_intermediate_gp = get_unique_gp(costanzo_diff_intermediate, "query_gene", "array_gene")

len(costanzo_ben_strict_gp) # 3472 unique interactions
len(costanzo_ctrl_strict_gp) # 3417
# len(costanzo_diff_strict_gp) # 1499
len(costanzo_ben_intermediate_gp) # 6715
len(costanzo_ctrl_intermediate_gp) # 6170
#len(costanzo_diff_intermediate_gp) # 2682

len(biogrid_gp.union(costanzo_ctrl_strict_gp).union(costanzo_ben_strict_gp)) # 441,520 unique GIs
len(biogrid_gp.union(costanzo_ctrl_intermediate_gp).union(costanzo_ben_intermediate_gp)) # 445,317

### Plot S10 Experimentally verified GI ----------------------------------------#
from venn import venn
sets = {"Costanzo Control": costanzo_ctrl_strict_gp,
        "Costanzo Benomyl": costanzo_ben_strict_gp,
        "BioGRID": biogrid_gp}
venn(sets) # total: 441520 unique GIs
plt.title("Experimentally verified genetic interactions")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Experimentally_verified_GIs_strict_venn_diagram.pdf")
plt.close()
del sets

sets = {"Costanzo Control": costanzo_ctrl_intermediate_gp,
        "Costanzo Benomyl": costanzo_ben_intermediate_gp,
        "BioGRID": biogrid_gp}
venn(sets) # total: 441520 unique GIs
plt.title("Experimentally verified genetic interactions")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Experimentally_verified_GIs_intermediate_venn_diagram.pdf")
plt.close()
del sets
### ---------------------------------------------------------------------------#

#### How many Costanzo & BioGRID interactions were identified by SHAP?
len(get_unique_gp(ts13_no_pav, "Gene1", "Gene2").intersection(biogrid_gp)) # 6391
len(get_unique_gp(ts13_no_pav, "Gene1", "Gene2").intersection(costanzo_ben_strict_gp)) # 1
len(get_unique_gp(ts13_no_pav, "Gene1", "Gene2").intersection(costanzo_ctrl_strict_gp)) # 59

verified_GIs_strict = {"ben":{}, "ctrl":{}, "biogrid":{}}
for name, data in variant_pairs.items(): # from line 633
	# strict threshold
	data_gp = get_unique_gp(data, "Gene1", "Gene2")
	print(f"Number of {name} interactions in benomyl GI network:\
		{len(data_gp.intersection(costanzo_ben_strict_gp))} and the control GI network:\
		{len(data_gp.intersection(costanzo_ctrl_strict_gp))} using the strict threshold.")
	verified_GIs_strict["ben"][name] = data_gp.intersection(costanzo_ben_strict_gp) 
	verified_GIs_strict["ctrl"][name] = data_gp.intersection(costanzo_ctrl_strict_gp)
	print(f"Number of {name} interactions in BioGRID:\
		{len(data_gp.intersection(biogrid_gp))}")
	verified_GIs_strict["biogrid"][name] = data_gp.intersection(biogrid_gp)
	#
	# intermediate threshold
	print(f"Using the intermediate threshold, benomyl GI network:\
		{len(data_gp.intersection(costanzo_ben_intermediate_gp))} and the control GI network:\
		{len(data_gp.intersection(costanzo_ctrl_intermediate_gp))}")
	del data_gp

'''
Number of snp_snp interactions in benomyl GI network:        1 and the control GI network:        59 using the strict threshold.
Number of snp_snp interactions in BioGRID:        6186 (unique gene-gene interactions)
Using the intermediate threshold, benomyl GI network:        3 and the control GI network:        74
Number of pav_pav interactions in benomyl GI network:        1 and the control GI network:        35 using the strict threshold.
Number of pav_pav interactions in BioGRID:        2053
Using the intermediate threshold, benomyl GI network:        1 and the control GI network:        41
Number of cnv_cnv interactions in benomyl GI network:        0 and the control GI network:        7 using the strict threshold.
Number of cnv_cnv interactions in BioGRID:        1273
Using the intermediate threshold, benomyl GI network:        0 and the control GI network:        8
Number of snp_pav interactions in benomyl GI network:        1 and the control GI network:        43 using the strict threshold.
Number of snp_pav interactions in BioGRID:        4454
Using the intermediate threshold, benomyl GI network:        3 and the control GI network:        57
Number of snp_cnv interactions in benomyl GI network:        1 and the control GI network:        1 using the strict threshold.
Number of snp_cnv interactions in BioGRID:        604        
Using the intermediate threshold, benomyl GI network:        1 and the control GI network:        1
Number of pav_cnv interactions in benomyl GI network:        1 and the control GI network:        1 using the strict threshold.
Number of pav_cnv interactions in BioGRID:        356
Using the intermediate threshold, benomyl GI network:        1 and the control GI network:        3
'''

# Add columns to Table S13
ts13_no_pav.insert(9, "In_BioGRID", ts13_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in biogrid_gp else 0, axis=1))
ts13_no_pav.insert(10, "In_Costanzo_Benomyl_Strict", ts13_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in costanzo_ben_strict_gp else 0, axis=1))
ts13_no_pav.insert(11, "In_Costanzo_Ctrl_Strict", ts13_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in costanzo_ctrl_strict_gp else 0, axis=1))
ts13_no_pav.insert(12, "In_Costanzo_Benomyl_Intermediate", ts13_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in costanzo_ben_intermediate_gp else 0, axis=1))
ts13_no_pav.insert(13, "In_Costanzo_Ctrl_Intermediate", ts13_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in costanzo_ctrl_intermediate_gp else 0, axis=1))
'''sanity check:
tmp = ts13_no_pav.apply(lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in biogrid_gp else 0, axis=1)
tmp = pd.concat([ts13_no_pav[['Gene1', 'Gene2', 'GI_Type']], tmp], axis=1, ignore_index=True)
len(get_unique_gp(tmp[tmp[3]==1], 0, 1)) # 6391
len(get_unique_gp(tmp.loc[(tmp[3]==1) & (tmp[2]=="SNP-SNP")], 0, 1)) # 6186
len(get_unique_gp(tmp.loc[(tmp[3]==1) & (tmp[2]=="PAV-PAV")], 0, 1)) # 2053
''' # looking good!

ts13_no_pav.to_csv(
	"Scripts/Data_Vis/Section_6/shap_interaction/Table_S13_benomyl_benchmark_RF_models_SHAP_interactions_expanded_pav_model_removed.txt",
	sep="\t", index=False)

### Plot C. Num of experimentally verified GIs identified by SHAP interactions-#
## Bar plot
summary_df = ts13_no_pav.groupby("GI_Type").sum()[["In_BioGRID", "In_Costanzo_Benomyl_Strict", "In_Costanzo_Ctrl_Strict"]]
# The below gives the number of variant-variant interactions, not unique gene-gene interactions
#          In_BioGRID  In_Costanzo_Benomyl_Strict  In_Costanzo_Ctrl_Strict
# GI_Type                                                                 
# CNV-CNV        1273                           0                        7
# PAV-CNV         359                           0                        1
# PAV-PAV        2086                           1                       35
# SNP-CNV         618                           0                        1
# SNP-PAV        5874                           1                       60
# SNP-SNP        6186                           1                       59

# Get the number of unique gene-gene interactions represented by the variant-variant interactions
# for source, inner_dict in verified_GIs_strict.items():
#     print(f"=== {source.upper()} ===")
#     for variant_pair_type, interaction_set in inner_dict.items():
#         print(f"{variant_pair_type}: {len(interaction_set)} interactions")
unique_counts = [
    {"Source": source, "Variant_Pair": pair, "Num_Interactions": len(gps)}
    for source, d in verified_GIs_strict.items()
    for pair, gps in d.items()
]
unique_counts_df = pd.DataFrame(unique_counts)
unique_counts_df = unique_counts_df.pivot(index="Variant_Pair", columns="Source", values="Num_Interactions")

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(9, 4))
summary_df["In_BioGRID"].plot(kind="bar", ax=axes[0], color="skyblue")
unique_counts_df["biogrid"].plot(kind="bar", ax=axes[0], color="darkblue", alpha=0.5)
axes[0].set_title("BioGRID")
axes[0].set_ylabel("Number of validated interactions")
axes[0].set_xlabel("Variant pair type")
summary_df["In_Costanzo_Benomyl_Strict"].plot(kind="bar", ax=axes[1], color="gold")
unique_counts_df["ben"].plot(kind="bar", ax=axes[1], color="darkorange", alpha=0.5)
axes[1].set_title("Costanzo Benomyl")
summary_df["In_Costanzo_Ctrl_Strict"].plot(kind="bar", ax=axes[2], color="lightgreen")
unique_counts_df["ctrl"].plot(kind="bar", ax=axes[2], color="darkgreen", alpha=0.5)
axes[2].set_title("Costanzo Control")
plt.suptitle("Number of experimentally verified genetic interactions identified by SHAP interactions")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7C_num_verified_GI_in_shap_interactions_bar_plot.pdf")
plt.close()

## Violin plot
fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True, figsize=(3, 7.5))
sns.violinplot(data=ts13_no_pav, x="In_BioGRID", y="Interaction",
	log_scale=10, fill=False, ax=ax[0])
sns.violinplot(data=ts13_no_pav, x="In_Costanzo_Benomyl_Strict",
	y="Interaction", log_scale=10, fill=False, ax=ax[1])
sns.violinplot(data=ts13_no_pav, x="In_Costanzo_Ctrl_Strict",
	y="Interaction", log_scale=10, fill=False, ax=ax[2])

# annotate stats
for i,col in enumerate(["In_BioGRID", "In_Costanzo_Benomyl_Strict", "In_Costanzo_Ctrl_Strict"]):
	annotator = Annotator(ax[i], [(0, 1)], data=ts13_no_pav,
		x=col, y="Interaction", order=[0, 1])
	annotator.configure(test="Mann-Whitney", text_format="star", loc="outside", comparisons_correction="bonferroni")
	annotator.apply_and_annotate()

for a in ax:
	a.axhline(ts13_no_pav["Interaction"].quantile(0.95), color="red", linestyle="--")
	a.axhline(ts13_no_pav["Interaction"].quantile(0.99), color="blue", linestyle="--")
'''
p-value annotation legend:
      ns: 5.00e-02 < p <= 1.00e+00
       *: 1.00e-02 < p <= 5.00e-02
      **: 1.00e-03 < p <= 1.00e-02
     ***: 1.00e-04 < p <= 1.00e-03
    ****: p <= 1.00e-04

0 vs. 1: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.459e-10 U_stat=1.300e+09
(<Axes: xlabel='In_BioGRID', ylabel='Interaction'>, [<statannotations.Annotation.Annotation object at 0x14b23f33f650>])

0 vs. 1: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:5.639e-01 U_stat=3.046e+05
(<Axes: xlabel='In_Costanzo_Benomyl_Strict', ylabel='Interaction'>, [<statannotations.Annotation.Annotation object at 0x14b231ea0290>])

0 vs. 1: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:2.186e-02 U_stat=1.531e+07
(<Axes: xlabel='In_Costanzo_Ctrl_Strict', ylabel='Interaction'>, [<statannotations.Annotation.Annotation object at 0x14b18942f350>])
'''

plt.ylabel("SHAP interaction score")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_7B_violin_log10_summed_shap_interaction_num_validated.pdf")
plt.close()

## Test for enrichment of validated interactions
# a. # of validated shap interactions
# b. # of non-validated shap interactions
# c. # of validated not in shap interactions
# d. # of non-validated and not in shap interactions

ts13_no_pav["Validated"] = ts13_no_pav[["In_BioGRID", "In_Costanzo_Benomyl_Strict", "In_Costanzo_Ctrl_Strict"]].any(axis=1).astype(int)
all_validated = biogrid_gp.union(costanzo_ctrl_strict_gp).union(costanzo_ben_strict_gp)

dpath = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"
benomyl_features = pd.read_csv(f"{dpath}/Features_one_variant_per_gene_benomyl_500ugml_snp_pav_cnv.txt", sep=",", header=None)
benomyl_features = pd.DataFrame(benomyl_features.iloc[0,:])
benomyl_features = pd.DataFrame(np.unique(benomyl_features))
benomyl_features[0] = benomyl_features[0].str.replace("^X", "", regex=True)
benomyl_features[0] = benomyl_features[0].str.replace("\.", "-", regex=True)
benomyl_features = benomyl_features.iloc[1:,:] # exclude header
all_pairs = pd.DataFrame(list(combinations(benomyl_features[0].values, 2)), columns=["Feature1", "Feature2"])
all_pairs.insert(2, "Gene1", all_pairs.apply(lambda x: feature2gene(x["Feature1"]), axis=1))
all_pairs.insert(3, "Gene2", all_pairs.apply(lambda x: feature2gene(x["Feature2"]), axis=1))
# all_pairs.isna().sum() # all good!

# Fisher's Exact Test
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
def enrichment_direction(k, n, C, G):
    """determine direction of enrichment
    if >= 1: + (overrepresented)
    if < 1: - (underrepresented)
    k: number of pairs in shap and are experimentally validated GIs
    n: total number of genes in shap
    C: total number of genes (in shap + background) and are known GIs
    G: total number of genes (in shap + background)"""
    return((k/C)/(n/G))

a = len(get_unique_gp(ts13_no_pav[ts13_no_pav.Validated==1], "Gene1", "Gene2")) # 6,395 gene-gene
# b = len(get_unique_gp(ts13_no_pav, "Gene1", "Gene2").difference(all_validated)) # 63,128 (this works too)
b = len(get_unique_gp(ts13_no_pav[ts13_no_pav.Validated==0], "Gene1", "Gene2")) # 63,128
c = len(all_validated.difference(get_unique_gp(ts13_no_pav, "Gene1", "Gene2"))) # 435,125 = 441,520 - 6,395
not_in_shap = get_unique_gp(all_pairs, "Gene1", "Gene2").difference(get_unique_gp(ts13_no_pav, "Gene1", "Gene2"))
not_in_validated = get_unique_gp(all_pairs, "Gene1", "Gene2").difference(all_validated) # 64,417
d = len(not_in_validated.intersection(not_in_shap)) # 1289 = 64,417 - 63,128
odds, pval = fisher_exact([[a, b],[c, d]]) # odds ratio = 0.00030009406022613993 ; pvalue = 0.0
enrichment_direction(k=a, n=a+b, C=a+c, G=a+b+c+d) # direction < 1, thus, underrepresented [known GIs not enriched]
### ---------------------------------------------------------------------------#

### Draw a genetic interaction network-----------------------------------------#
def draw_networks(best_rf_fs_gi, key, d, model_dir, env, data_type):
	# Draw the networks of the known GIs identified by SHAP-based interactions
	
	ts13_no_pav.loc[ts13_no_pav.Interaction > ts13_no_pav.Interaction.quantile(0.95)]
	# Get the identified GIs
	# gi = best_rf_fs_gi[key][data_type][env]
	
	# Load the SHAP-based interaction scores
	shap = pd.read_csv(
		glob.glob(f"{d}/{data_type}/{model_dir}/shap_interaction_scores_{data_type.lower()}_{env}_*_summed.txt")[0],
		sep="\t")
	
	# map the feature names to the gene names
	if data_type == "SNP":
		shap["Gene1"] = shap.Feature1.map(map_snps.set_index("snp").gene)
		shap["Gene2"] = shap.Feature2.map(map_snps.set_index("snp").gene)
		# drop interactions with at least one intergenic snp, we are only looking for gene-gene interactions
		shap = shap.loc[(shap.Gene1 != "intergenic") & (shap.Gene2 != "intergenic"),:]
		# drop snps with multiple gene matches
		shap = shap.loc[~(shap.Gene1.str.contains(",")) & ~(shap.Gene2.str.contains(",")),:]
	else:
		shap.Feature1 = shap.Feature1.apply(lambda x: re.sub("^X", "", x))
		shap.Feature1 = shap.Feature1.apply(lambda x: re.sub("\.", "-", x))
		shap.Feature2 = shap.Feature2.apply(lambda x: re.sub("^X", "", x))
		shap.Feature2 = shap.Feature2.apply(lambda x: re.sub("\.", "-", x))
		shap["Gene1"] = shap.Feature1.map(map_orfs.set_index("orf").gene)
		shap["Gene2"] = shap.Feature2.map(map_orfs.set_index("orf").gene)
		shap["Gene1"] = shap["Gene1"].fillna(shap["Feature1"])
		shap["Gene2"] = shap["Gene2"].fillna(shap["Feature2"])
	
	# get the known GIs identified by SHAP-based interactions
	idx = []
	for pair in gi:
		pair = list(pair)
		idx.append(shap.loc[(shap.Gene1 == pair[0]) & (shap.Gene2 == pair[1]) |\
			(shap.Gene1 == pair[1]) & (shap.Gene2 == pair[0])].index[0])
	
	identified_gi = shap.loc[idx,:]
	
	# Build the network
	G = nx.from_pandas_edgelist(identified_gi, "Gene1", "Gene2", "Interaction",
		create_using=nx.Graph())
	
	# Extract edge colors
	edge_weights = identified_gi["Interaction"].values
	norm = mcolors.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))  # Normalize color scale
	cmap = cm.Blues  # Choose colormap
	sm = cm.ScalarMappable(norm=norm, cmap=cmap)  # Create a ScalarMappable for the colorbar
	sm.set_array([])  # Required for colorbar to work
	
	# Draw the network
	fig = plt.figure(figsize=(10,10))
	gs = gridspec.GridSpec(1, 2, width_ratios=[10, 0.5]) # allocate space for colorbar
	ax = plt.subplot(gs[0]) # main plot area
	cax = plt.subplot(gs[1]) # colorbar axis
	pos = nx.spring_layout(G)  # Compute node positions
	nx.draw_networkx_edges(G, pos, edge_color=edge_weights, edge_cmap=cmap, edge_vmin=min(edge_weights), edge_vmax=max(edge_weights), ax=ax)
	nx.draw_networkx_nodes(G, pos, node_size=12, ax=ax)
	nx.draw_networkx_labels(G, pos, font_size=12, ax=ax)
	
	# color the nodes that correspond to benchmark genes in red
	if data_type == "SNP":
		benchmark_genes = map_snps.loc[map_snps[["Benomyl", "Caffeine", "CuSO4", "Sodium_meta-arsenite"]].sum(axis=1) > 0,"gene"].unique()
	if data_type in ["PAV", "CNV"]:
		benchmark_genes = map_ords.loc[map_orfs[["Benomyl", "Caffeine", "CuSO4", "Sodium_meta-arsenite"]].sum(axis=1) > 0,"gene"].unique()
	benchmark_genes = set(benchmark_genes).intersection(G.nodes()) # search for benchmark_genes values in G.nodes()
	node_colors = ["red" if node in benchmark_genes else "blue" for node in G.nodes()]
	nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=12, ax=ax)
	
	cbar = plt.colorbar(sm, cax=cax) # Add a colorbar legend
	cbar.set_label("SHAP Interaction Score")
	ax.set_title(f"{data_type} {env} {model_dir}")
	plt.tight_layout()
	
	tosave = "Scripts/Data_Vis/Section_6/shap_interaction/"
	plt.savefig(f"{tosave}/shap_interaction_network_{data_type.lower()}_{env}_{model_dir}_{key}.pdf")
	plt.close()
	
	return G
### ---------------------------------------------------------------------------#

### Which verified GIs were identified by multiple variant-variant interactions?
from collections import defaultdict
overlap_counts = defaultdict(lambda: defaultdict(set))

for sample, categories in verified_GIs_strict.items():
	for category, pair_set in categories.items():
		for pair in pair_set:
			overlap_counts[pair][sample].add(category)

rows = []
for pair, sample_dict in overlap_counts.items():
	gene1, gene2 = sorted(pair)
	for sample, categories in sample_dict.items():
		if len(categories) > 1:
			rows.append({
				"gene1": gene1, "gene2": gene2,
				"sample": sample, "n_categories": len(categories),
				"categories": "; ".join(sorted(categories))})
			# print(f"{pair} appears in {len(categories)} categories for {sample}: {sorted(categories)}")

df = pd.DataFrame(rows)
df = df.sort_values(by=['sample', 'n_categories'], ascending=[True, False])
df.groupby("sample").head(20)
'''these are not necessarily the top [or the best] because many gene pairs were
identified by 4 variant pair types, for example. I would have to sort by max
interaction score or an average or something like that.'''
df.shape # (5523, 5)

df.loc[(df.gene1=='YGR078C') & (df.gene2=="YNL153C")] # was found by all three datasets (biogrid, ben, and ctrl)

# Which gene pairs have the highest number of feature interactions?
ts13_no_pav.groupby(["Gene1", "Gene2"]).count().sort_values("Model", ascending=False)

# Which gene pairs have the highest interaction scores?
ts13_no_pav.sort_values("Interaction", ascending=False).iloc[:20,:] # highest shap interaction scores



