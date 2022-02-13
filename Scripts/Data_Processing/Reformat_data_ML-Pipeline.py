# Description: Reformat data to suit Shiu Lab ML-Pipeline for regression
# Author: Kenia Segura Aba

import os
import pandas as pd
import datatable as dt
# Load genotype and phenotype matrices
os.chdir("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/")
geno_data="SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv"
pheno_data="Phenotype_value_383_common_accessions_2017_Grimm.csv"
geno = dt.fread(geno_data) # read in genotype data
geno = geno.to_pandas() # convert to pandas dataframe
geno = geno.sort_values(by=geno.columns[0], axis=0) # sort values by sample ID
geno = geno.set_index(geno.columns[0], drop=True)
pheno = pd.read_csv(pheno_data, index_col=0)
#pheno=pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/pheno.csv", index_col=0)
#geno=pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/geno.csv", index_col=0)

# traits
traits = pheno.columns
pheno.head, geno.head, traits

import copy
for i in range(len(traits)):
    trait = traits[i]
    g = copy.deepcopy(geno) # create a copy of geno.csv
    g.insert(0, "Y", pheno[trait]) # insert trait column as "Class" into geno matrix
    g.to_csv(geno_data+"_"+trait+"_rf.csv", index=True)
    #g.to_csv("/mnt/scratch/seguraab/yeast_project/"+geno_data+"_"+trait+"_rf.csv", index=True)
