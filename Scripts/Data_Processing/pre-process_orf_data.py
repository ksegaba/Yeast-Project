#!/usr/bin/env python3

'''Pre-processing of raw open reading frame data from Peter et al., 2018
'''

import os
import pandas as pd
import matplotlib.pyplot as plt

# First I extracted diploid isolates from genesMatrix_CopyNumber.tab
os.system("python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Processing/filter_diploid.py \
    -file genesMatrix_CopyNumber.tab -labels Isolate_ploidy.tsv")

# Keep only isolates with phenotype information
cnv = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/0_raw_data/genesMatrix_CopyNumber.tab_diploid.txt', sep='\t')
pheno = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv')
out = cnv[cnv.STD_name.isin(pheno['ID'])]
out.rename(columns={'STD_name': 'ID'}, inplace=True)
out.to_csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs.csv', index=False)

# Remove CNVs with missing values (7,796 dropped to 7,708 ORFs)
cnv = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs.csv', index_col=0)
cnv.isna().sum().apply(lambda x: x/750).plot(kind='hist', bins=10)
cnv.isna().sum().apply(lambda x: x/750).sum() # 88 cnvs missing in all 750 isolates
cnvs_to_drop = cnv.columns[cnv.isna().sum() == 750]
plt.xlabel('Percentage of isolates with missing data')
plt.savefig('/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Processing/ORFs_with_NAs.png')
plt.close()

# S6 File CNV data
cnv.dropna(axis=1).to_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv")

# Remove PAVs with missing values
pav = pd.read_csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/0_raw_data/genesMatrix_PresenceAbsence.tab', sep='\t')
pav = pav[pav.STD_name.isin(pheno['ID'])] # diploid isolates only
pav.rename(columns={'STD_name': 'ID'}, inplace=True)
pav.isna().sum().apply(lambda x: x/750).sum() # 88 PAVs missing in all 750 isolates
pav.columns[pav.isna().sum() == 750].isin(cnvs_to_drop).sum() # the same CNVs are missing in PAVs

# S5 File PAV data
pav.dropna(axis=1).to_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv", index=False)
