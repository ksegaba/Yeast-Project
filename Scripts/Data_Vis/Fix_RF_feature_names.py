"""Feature names in geno.csv are in the order as they appear in the VCF file,
but, PLINK changed the order when converting the VCF to fastPHASE format with
minor and major allele encodings.

I fixed the code fastPHASE_to_csv.py to get the correct feature order. Now, in
this script, I do the following:

1. I determine what the "wrong" feature names in geno.csv are supposed
to be.
    a. Correct the geno.csv feature names --> geno_corrected.csv
2. Write new feature importance files for complete and optimized RF models.
"""

import re
import datatable as dt
import pandas as pd
from tqdm import tqdm

# --- 1. Determine what the actual feature names are
# Read in geno.csv and the corrected geno.csv
path = "/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018"
wrong_geno = dt.fread(f"{path}/geno.csv").to_pandas()
right_geno = dt.fread(
    f"{path}/0_raw_data/1011Matrix_diploid750_biallelic_maf0.05_09102025_genotypes.csv").to_pandas()

assert wrong_geno.shape == right_geno.shape  # great!
assert wrong_geno.ID.equals(right_geno["index"])  # great!

mapping = pd.concat([pd.Series(wrong_geno.columns[1:]),
                    pd.Series(right_geno.columns[1:])], axis=1)
mapping.columns = ["wrong_geno.csv_feature", "actual_feature"]
mapping.to_csv(f"{path}/0_raw_data/mapping_of_geno.csv_feature_names_to_actual_feature_names_09102025.csv",
               index=False)
mapping_dict = mapping.set_index("wrong_geno.csv_feature")[
    "actual_feature"].to_dict()

# a. Correct the geno.csv feature names --> geno_corrected.csv
wrong_geno.rename(columns=mapping_dict, inplace=True)
wrong_geno.to_csv(f"{path}/geno_corrected.csv", index=False)

# --- 2. Write new feature importance files for complete and optimized RF models
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results"
complete_model_list = pd.read_csv(
    f"{path}/baseline/snp_best_models.txt", header=None)
optimized_model_list = pd.read_csv(
    f"{path}/fs/snp_best_fs_models.txt", header=None)

path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/SNP"
snp_shap_optimized_files = [f"{path}/fs/{file}" for file in os.listdir(
    path + "/fs") if file.startswith("SHAP_values_sorted_average_Y")]
snp_shap_complete_files = [f"{path}/baseline/{file}" for file in os.listdir(
    path + "/baseline") if file.startswith("SHAP_values_sorted_average_Y")]


def correct_feature_names_in_gini(model):
    """Add a column with the actual feature names to the Gini importance files"""
    imp_filename = re.sub("_models_rep_[0-9]+.pkl", "_imp", model)
    path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results"
    # read the gini importance file
    if "baseline" in model:  # RF models built with all features
        imp = pd.read_csv(
            f"{path}/baseline/{imp_filename}", sep="\t", header=0)
        out_filename = f"{path}/baseline/{imp_filename}_with_actual_feature_names_09102025"
    elif "rf_top" in model:  # optimized RF models
        imp = pd.read_csv(
            f"{path}/fs/{imp_filename}", sep="\t", header=0)
        out_filename = f"{path}/fs/{imp_filename}_with_actual_feature_names_09102025"
    # insert the actual feature names and save the new file
    imp.insert(0, "actual_feature", imp["Unnamed: 0"].map(mapping_dict))
    imp.rename(columns={"Unnamed: 0": "wrong_geno.csv_feature"}, inplace=True)
    imp.to_csv(out_filename, sep="\t", index=False)
    del imp_filename, imp, out_filename


complete_model_list[0].apply(correct_feature_names_in_gini)
optimized_model_list[0].apply(correct_feature_names_in_gini)


def correct_feature_names_in_shap(shap_filename):
    """Add a column with the actual feature names to the SHAP values files"""
    # read the shap values file
    shap = pd.read_csv(shap_filename, sep="\t", header=0)
    if "baseline" in shap_filename:  # RF models built with all features
        out_filename = f"{shap_filename.replace('.txt','')}_with_actual_feature_names_09102025"
    elif "rf_fs" in shap_filename:  # optimized RF models
        out_filename = f"{shap_filename.replace('.txt','')}_with_actual_feature_names_09102025"
    # insert the actual feature names and save the new file
    shap.insert(0, "actual_feature", shap["Unnamed: 0"].map(mapping_dict))
    shap.rename(columns={"Unnamed: 0": "wrong_geno.csv_feature"}, inplace=True)
    shap.to_csv(out_filename, sep="\t", index=False)
    del shap_filename, shap, out_filename


for i in tqdm(range(len(snp_shap_optimized_files))):
    correct_feature_names_in_shap(snp_shap_optimized_files[i])
    correct_feature_names_in_shap(snp_shap_complete_files[i])
