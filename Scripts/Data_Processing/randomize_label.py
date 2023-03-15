"""
Script to randomize the phenotype matrix
"""
__author__ = "Kenia Segura Abá"

import os
import pandas as pd

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018")

# read in fitnes data
pheno = pd.read_csv("pheno.csv")

# randomized dataframe
pheno2 = pd.DataFrame(pheno.ID)

# randomize columns individually
for col in pheno.columns[1:]:
    pheno2[col] = pheno[col].sample(frac=1, random_state=123).to_list()

pheno2.head()
pheno2.to_csv("pheno_randomized.csv", index=None)