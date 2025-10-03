#!/usr/bin/env python3
# 9/16/2025: This script seems to be old in general. I didn't use these feature sets
# in the final SHAP interaction models.

import pickle
import re
import swifter
import os
import pandas as pd
import datatable as dt
import numpy as np
from tqdm import tqdm
from sklearn.metrics import r2_score

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
# TABLE S11
# This is also SHAP interaction
# Actually, the code below might be old, but need to verify-------Sep. 10, 2025
################################################################################
"""After running features selection on the RF models used for SHAP interaction
analysis, concatenate the literature genes to the top features and run a new set
of RF models. Also run the models with the label randomized. The purpose of
these is to obtain SHAP gene-gene interactions between top genes and literature
genes and also the feature SHAP values. 
First kind: top FS features plus unimportant benomyl, caffeine, cuso4, or sma genes
Second kind: top FS features plus unimportant non-lit genes
Third kind: only top benchmark genes
Fourth kind: only top non-benchmark genes
"""

# SNP and ORF gene maps
# some models may include intergenic snps and snps that mapped to multiple genes in B_nl
map_snps = pd.read_csv(
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")

# Read in the RF FS performance results
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models_fs"
# XGB FS performance results
res = pd.read_csv(f"{d}/RESULTS_reg.txt", sep="\t")
res = res.loc[res.DateTime.str.contains("03-26"),]

ben_snp = map_snps.loc[map_snps.Benomyl == 1, ["snp", "gene"]]
ben_orf = map_orfs.loc[map_orfs.Benomyl == 1, ["orf", "gene"]]
caf_snp = map_snps.loc[map_snps.Caffeine == 1, ["snp", "gene"]]
caf_orf = map_orfs.loc[map_orfs.Caffeine == 1, ["orf", "gene"]]
cu_snp = map_snps.loc[map_snps.CuSO4 == 1, ["snp", "gene"]]
cu_orf = map_orfs.loc[map_orfs.CuSO4 == 1, ["orf", "gene"]]
sma_snp = map_snps.loc[map_snps["Sodium_meta-arsenite"] == 1, ["snp", "gene"]]
sma_orf = map_orfs.loc[map_orfs["Sodium_meta-arsenite"] == 1, ["orf", "gene"]]

# Now create the feature lists of top features + individual literature gene lists
d2 = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/models"
# models_to_run_a = open(f"{d}/fs_plus_lit_genes_models_to_run.txt", "a") # list of models to run
# models_to_run_b = open(f"{d}/fs_plus_unimportant_non_lit_genes_models_to_run.txt", "a")
models_to_run_c = open(f"{d}/only_lit_genes_models_to_run.txt", "a")
models_to_run_d = open(f"{d}/only_top_non_lit_genes_models_to_run.txt", "a")

best_fs_models = open(f"{d}/best_rf_fs_models.txt", "a")
test_set = pd.read_csv(
    '~/Shiu_Lab/Project/Data/Peter_2018/Test.txt', header=None)

for data_type in ["snp", "pav", "cnv"]:
    for env in target_envs:
        if env in ["YPDCAFEIN40", "YPDCAFEIN50"]:
            lit_env = "caffeine"
        elif env == "YPDBENOMYL500":
            lit_env = "benomyl"
        elif env == "YPDCUSO410MM":
            lit_env = "cuso4"
        elif env == "YPDSODIUMMETAARSENITE":
            lit_env = "sodium_meta-arsenite"
        print(data_type, env, lit_env)

        # Determine which model has the top validation R2 score
        res_env = res.loc[(res.Tag.str.contains(f"{env}") &
                           res.Tag.str.contains(f"{data_type}")), :]
        top = res_env.loc[res_env["r2_val"] ==
                          res_env["r2_val"].max(), "FeatureNum"].values[0]

        # Write the best training repetition model to a file (need for shap interaction calculation)
        tag = res_env.loc[res_env['r2_val'] ==
                          res_env['r2_val'].max(), 'Tag'].values[0]
        # predicted trait values for each training repetition
        scores = pd.read_csv(f"{d}/{tag}_scores.txt", sep="\t", index_col=0)
        scores_val = scores.loc[~scores.index.isin(test_set[0]), :].drop(
            columns=["Mean", "stdev"])  # validation set
        rep = int(scores_val.iloc[:, 1:].apply(lambda x: r2_score(
            scores_val["Y"], x), axis=0).idxmax().split("_")[1])  # best training repetition
        best_fs_models.write(f"{d}/{tag}_models_rep_{rep-1}.pkl\n")

        # Read in the top model feature importances
        imp = pd.read_csv(
            f"{d}/{data_type}_{env}_max_gini_from_RF_baseline_imp_top_{top}_imp", index_col=0, sep="\t")
        imp = imp.loc[imp.mean_imp != 0.0, :]  # drop unimportant features

        # Get literature genes that map to baseline model
        imp_base = pd.read_csv(
            f"{d2}/{data_type}_{env}_max_gini_from_RF_baseline_imp_imp", index_col=0, sep="\t")

        if data_type != "snp":
            imp.index = imp.apply(lambda x: re.sub("^X", "", x.name), axis=1)
            imp.index = imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
            imp_base.index = imp_base.apply(
                lambda x: re.sub("^X", "", x.name), axis=1)
            imp_base.index = imp_base.apply(
                lambda x: re.sub("\.", "-", x.name), axis=1)

        if lit_env == "benomyl":
            if data_type == "snp":
                # Determine literature and non-literature genes in the top features
                T_l = imp.loc[imp.index.isin(
                    ben_snp.iloc[:, 0]), :]  # snp features
                T_l = T_l["mean_imp"].to_frame().merge(
                    ben_snp, how="left", left_index=True, right_on="snp")
                # snp features  not in benomyl genes
                T_nl = imp.loc[~imp.index.isin(ben_snp.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")

                # Identify the unimportant model features found in ben_snp and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(ben_snp.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
                # features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp") # fs + unimportant lit genes
                features_c = pd.concat([T_l, B_l]).set_index(
                    "snp")  # lit genes only
                features_d = T_nl.copy(deep=True).set_index(
                    "snp")  # important non-lit genes only
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(ben_snp.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
                # features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp") # fs + unimportant non-lit genes

            else:
                # Determine literature and non-literature genes in the top features
                T_l = imp.loc[imp.index.isin(
                    ben_orf.iloc[:, 0]), :]  # orf features
                T_l = T_l["mean_imp"].to_frame().merge(
                    ben_orf, how="left", left_index=True, right_on="orf")
                # orf features not in benomyl genes
                T_nl = imp.loc[~imp.index.isin(ben_orf.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")

                # Identify the unimportant model features found in ben_orf and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(ben_orf.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_a = pd.concat([T_l, T_nl, B_l])
                features_c = pd.concat([T_l, B_l])
                features_d = T_nl.copy(deep=True)
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(ben_orf.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_b = pd.concat([T_l, T_nl, B_nl])

        elif lit_env == "caffeine":
            if data_type == "snp":
                # Determine literature and non-literature genes in the top features
                T_l = imp.loc[imp.index.isin(
                    caf_snp.iloc[:, 0]), :]  # snp features
                T_l = T_l["mean_imp"].to_frame().merge(
                    caf_snp, how="left", left_index=True, right_on="snp")
                # snp features not in caf genes
                T_nl = imp.loc[~imp.index.isin(caf_snp.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")

                # Identify the unimportant model features found in caf_snp and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(caf_snp.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
                # features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp")
                features_c = pd.concat([T_l, B_l]).set_index("snp")
                features_d = T_nl.copy(deep=True).set_index("snp")
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(caf_snp.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
                # features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp")

            else:
                # Determine literature and non-literature genes in the top features
                T_l = imp.loc[imp.index.isin(
                    caf_orf.iloc[:, 0]), :]  # orf features
                T_l = T_l["mean_imp"].to_frame().merge(
                    caf_orf, how="left", left_index=True, right_on="orf")
                # orf features  not in caf genes
                T_nl = imp.loc[~imp.index.isin(caf_orf.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")

                # Identify the unimportant model features found in caf_orf and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(caf_orf.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_a = pd.concat([T_l, T_nl, B_l])
                features_c = pd.concat([T_l, B_l])
                features_d = T_nl.copy(deep=True)
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(caf_orf.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_b = pd.concat([T_l, T_nl, B_nl])

        elif lit_env == "cuso4":
            if data_type == "snp":
                # Determine literature and non-literature genes in the top features
                T_l = imp.loc[imp.index.isin(
                    cu_snp.iloc[:, 0]), :]  # snp features
                T_l = T_l["mean_imp"].to_frame().merge(
                    cu_snp, how="left", left_index=True, right_on="snp")
                # snp features not in cuso4 genes
                T_nl = imp.loc[~imp.index.isin(cu_snp.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")

                # Identify the unimportant model features found in cu_snp and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(cu_snp.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
                # features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp")
                features_c = pd.concat([T_l, B_l]).set_index("snp")
                features_d = T_nl.copy(deep=True).set_index("snp")
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(cu_snp.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")
                # features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp")

            else:
                # Determine literature and non-literature genes in the top features
                T_l = imp.loc[imp.index.isin(
                    cu_orf.iloc[:, 0]), :]  # orf features
                T_l = T_l["mean_imp"].to_frame().merge(
                    cu_orf, how="left", left_index=True, right_on="orf")
                # orf features not in cuso4 genes
                T_nl = imp.loc[~imp.index.isin(cu_orf.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")

                # Identify the unimportant model features found in cu_orf and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(cu_orf.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_a = pd.concat([T_l, T_nl, B_l])
                features_c = pd.concat([T_l, B_l])
                features_d = T_nl.copy(deep=True)
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(cu_orf.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_b = pd.concat([T_l, T_nl, B_nl])

        elif lit_env == "sodium_meta-arsenite":
            if data_type == "snp":
                # Determine literature and non-literature genes in the top features
                # snp features in sma genes
                T_l = imp.loc[imp.index.isin(sma_snp.iloc[:, 0]), :]
                T_l = T_l["mean_imp"].to_frame().merge(
                    sma_snp, how="left", left_index=True, right_on="snp")
                # snp features not in sma genes
                T_nl = imp.loc[~imp.index.isin(sma_snp.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")

                # Identify the unimportant model features found in sma_snp and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(sma_snp.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")  # snp features not top and in sma genes
                # features_a = pd.concat([T_l, T_nl, B_l]).set_index("snp")
                features_c = pd.concat([T_l, B_l]).set_index("snp")
                features_d = T_nl.copy(deep=True).set_index("snp")
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(sma_snp.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_snps[["snp", "gene"]], how="left", left_index=True, right_on="snp")  # snp features not top and not in sma genes
                # features_b = pd.concat([T_l, T_nl, B_nl]).set_index("snp")

            else:
                # Determine literature and non-literature genes in the top features
                T_l = imp.loc[imp.index.isin(
                    sma_orf.iloc[:, 0]), :]  # orf features
                T_l = T_l["mean_imp"].to_frame().merge(
                    sma_orf, how="left", left_index=True, right_on="orf")
                # orf features not in sma genes
                T_nl = imp.loc[~imp.index.isin(sma_orf.iloc[:, 0]), :]
                T_nl = T_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")

                # Identify the unimportant model features found in sma_orf and combine with imp
                B_l = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                   (imp_base.index.isin(sma_orf.iloc[:, 0])), :]
                B_l = B_l["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_a = pd.concat([T_l, T_nl, B_l])
                features_c = pd.concat([T_l, B_l])
                features_d = T_nl.copy(deep=True)
                B_nl = imp_base.loc[(~imp_base.index.isin(imp.index)) &
                                    (~imp_base.index.isin(sma_orf.iloc[:, 0])), :]
                B_nl = B_nl["mean_imp"].to_frame().merge(
                    map_orfs[["orf", "gene"]], how="left", left_index=True, right_on="orf")
                # features_b = pd.concat([T_l, T_nl, B_nl])

        # Ensure that the genes are unique
        if data_type != "snp":
            # Replace missing gene values with the orf ID
            # features_a.gene = features_a.gene.fillna(features_a.orf)
            # features_b.gene = features_b.gene.fillna(features_b.orf)
            features_c.gene = features_c.gene.fillna(features_c.orf)
            features_d.gene = features_d.gene.fillna(features_d.orf)
            # features_a.orf = "X" + features_a.orf.str.replace("-", ".")
            # features_b.orf = "X" + features_b.orf.str.replace("-", ".")
            features_c.orf = "X" + features_c.orf.str.replace("-", ".")
            features_d.orf = "X" + features_d.orf.str.replace("-", ".")

            features_c = features_c.set_index("orf")
            features_d = features_d.set_index("orf")

        # features_a = features_a.loc[features_a.index.isin(features_a.groupby("gene")["mean_imp"].idxmax()),:] # one feature per gene
        # features_b = features_b.loc[features_b.index.isin(features_b.groupby("gene")["mean_imp"].idxmax()),:]
        features_c = features_c.loc[features_c.index.isin(
            features_c.groupby("gene")["mean_imp"].idxmax()), :]
        features_d = features_d.loc[features_d.index.isin(
            features_d.groupby("gene")["mean_imp"].idxmax()), :]
        # EXCLUDE SNPs THAT MAPPED TO MORE THAN ONE GENE, SINCE ONE COULD BE A BENCHMARK GENE
        features_d = features_d.loc[~features_d.gene.str.contains(","), :]
        # if (features_a.gene.nunique() == len(features_a)) & \
        # (features_b.gene.nunique() == len(features_b)) & \
        if (features_c.gene.nunique() == len(features_c)) & \
                (features_d.gene.nunique() == len(features_d)):
            # np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_{lit_env}_lit_genes_info", features_a, fmt="%s")
            # np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_{lit_env}_lit_genes", features_a.iloc[:,1], fmt="%s")
            # models_to_run_a.write(f"{d}/{data_type}_{env}_top_{top}_plus_{lit_env}_lit_genes\n")
            # np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_{lit_env}_genes_info", features_b, fmt="%s")
            # np.savetxt(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_{lit_env}_genes", features_b.iloc[:,1], fmt="%s")
            # models_to_run_b.write(f"{d}/{data_type}_{env}_top_{top}_plus_unimportant_non_{lit_env}_genes\n")

            # Before saving, ensure feature_c and features_d have the same number of features
            np.savetxt(f"{d}/{data_type}_{env}_{lit_env}_lit_genes_info",
                       features_c.reset_index(), fmt="%s")
            np.savetxt(f"{d}/{data_type}_{env}_top_non_{lit_env}_genes_info",
                       features_d.reset_index(), fmt="%s")

            # ensure the same number of features in both models
            if len(features_c) > len(features_d):
                features_c = features_c.sort_values(
                    "mean_imp", ascending=False).iloc[:len(features_d), :]
            elif len(features_c) < len(features_d):
                features_d = features_d.sort_values(
                    "mean_imp", ascending=False).iloc[:len(features_c), :]
            print(len(features_c), len(features_d))

            np.savetxt(f"{d}/{data_type}_{env}_{lit_env}_lit_genes",
                       features_c.index, fmt="%s")
            # models_to_run_c.write(f"{d}/{data_type}_{env}_{lit_env}_lit_genes\n")
            np.savetxt(
                f"{d}/{data_type}_{env}_top_non_{lit_env}_genes", features_d.index, fmt="%s")
            # models_to_run_d.write(f"{d}/{data_type}_{env}_top_non_{lit_env}_genes\n")

        del res_env, top, imp, imp_base, T_l, T_nl, B_l, B_nl, features_c, features_d

# models_to_run_a.close()
# models_to_run_b.close()
models_to_run_c.close()
models_to_run_d.close()
best_fs_models.close()

# Randomize the label and add it to the top features + combined literature gene list datasets
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
pheno_train = pheno.loc[~pheno.index.isin(test[0]), :]
dr = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized"
for i in tqdm(range(100)):
    randomized = pd.DataFrame(index=pheno_train.index)
    for env in target_envs:
        randomized[env] = pheno_train[env].sample(
            frac=1, random_state=i).reset_index(drop=True).values
    randomized = pd.concat([randomized, pheno.loc[pheno.index.isin(
        test[0]), target_envs]])  # add the test set back
    randomized.to_csv(f"{dr}/pheno_randomized_{i}.csv")
    del randomized
