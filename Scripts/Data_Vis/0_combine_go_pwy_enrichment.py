#!/usr/bin/env python3
"""
Table S6
Combine GO/pathway enrichment results into an excel file
Arguments:
    path: path to csv files
    pattern: csv file name pattern to match
    dtype: data type (SNP, ORF, CNV)
    atype: annotation type (go, pwy)
    save: path to save excel file to
Usage:
    python combine_go_pwy_enrichment.py <path> <pattern> <dtype> <atype>
Examples:
    python Scripts/Data_Vis/0_SuppDatFile_14_15_combine_go_pwy_enrichment.py \
        Scripts/Genomic_Prediction_RF/PWY_Enrichment/SNPs_fs/ \
        ORA_PWY_Genes_ SNP pwy Scripts/Data_Vis/Section_4/
    python Scripts/Data_Vis/0_SuppDatFile_14_15_combine_go_pwy_enrichment.py \
        Scripts/Genomic_Prediction_RF/PWY_Enrichment/PAVs_fs/ \
        ORA_PWY_Genes_[A-Z0-9]+_pav PAV pwy Scripts/Data_Vis/Section_4/
    python Scripts/Data_Vis/0_SuppDatFile_14_15_combine_go_pwy_enrichment.py \
        Scripts/Genomic_Prediction_RF/PWY_Enrichment/CNVs_fs/ \
        ORA_PWY_Genes_[A-Z0-9]+_cnv CNV pwy Scripts/Data_Vis/Section_4/

    python Scripts/Data_Vis/0_SuppDatFile_14_15_combine_go_pwy_enrichment.py \
        Scripts/Genomic_Prediction_RF/GO_Enrichment/SNPs_fs/ \
        ORA_Genes_ SNP go Scripts/Data_Vis/Section_4/
    python Scripts/Data_Vis/0_SuppDatFile_14_15_combine_go_pwy_enrichment.py \
        Scripts/Genomic_Prediction_RF/GO_Enrichment/PAVs_fs/ \
        ORA_Genes_[A-Z0-9]+_pav PAV go Scripts/Data_Vis/Section_4/
    python Scripts/Data_Vis/0_SuppDatFile_14_15_combine_go_pwy_enrichment.py \
        Scripts/Genomic_Prediction_RF/GO_Enrichment/CNVs_fs/ \
        ORA_Genes_[A-Z0-9]+_cnv CNV go Scripts/Data_Vis/Section_4/
Returns:
    excel file with all GO/pathway enrichment results
    e.g. SNP_pathway_enrichment.xlsx or SNP_GO_enrichment.xlsx
"""

import os
import sys
import re
import pandas as pd
from pathlib import Path

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

# Arguments
path = sys.argv[1] # path to csv files
pattern = sys.argv[2] # csv file name pattern to match
dtype = sys.argv[3] # data type (SNP, ORF, CNV)
atype = sys.argv[4] # annotation type (go, pwy)
save = sys.argv[5] # path to save excel file to

# Pathway names
pwy_info = pd.read_csv("../Co-function/Data/MetaCyc/All_pathways_S288c_names.txt", sep="\t")
print(pwy_info.columns)

# Create and save excel file
if (atype=="pwy" or atype=="pathway"):
    excel_file = Path(save) / f"{dtype}_pathway_enrichment.xlsx"
    writer = pd.ExcelWriter(excel_file)
    pd.DataFrame({}).to_excel(writer, sheet_name="All") # create empty sheet
    writer.close()
    with pd.ExcelWriter(excel_file, mode="a", if_sheet_exists="overlay") as writer:
        for csvfilename in os.listdir(path):
            if re.search(pattern, str(csvfilename)):
                print(csvfilename)
                df = pd.read_csv(Path(path) / csvfilename, sep="\t")
                print(df.shape)
                df = df.loc[df["p.val"]<=0.05,:] # remove non-significant pathways
                print(df.shape)
                df['PWY'] = df['PWY'].str.strip()
                df = df.merge(pwy_info, left_on="PWY", right_on="Object ID")
                sheet_name = os.path.splitext(csvfilename)[0].split("_")[3]
                df.to_excel(writer,sheet_name=sheet_name)
                df.insert(0, "Environment", sheet_name)
                df.to_excel(writer, sheet_name="All", index=False, header=False, startrow=writer.sheets['All'].max_row)

if (atype=="go" or atype=="GO"):
    excel_file = Path(save) / f"{dtype}_GO_enrichment.xlsx"
    writer = pd.ExcelWriter(excel_file)
    pd.DataFrame({}).to_excel(writer, sheet_name="All") # create empty sheet
    writer.close()
    with pd.ExcelWriter(excel_file, mode="a", if_sheet_exists="overlay") as writer:
        for csvfilename in os.listdir(path):
            if re.search(pattern, str(csvfilename)):
                print(csvfilename)
                df = pd.read_csv(Path(path) / csvfilename, sep="\t")
                print(df.shape)
                df = df.loc[df["p.val"]<=0.05,:] # remove non-significant pathways
                print(df.shape)
                sheet_name = os.path.splitext(csvfilename)[0].split("_")[2]
                df.to_excel(writer,sheet_name=sheet_name)
                df.insert(0, "Environment", sheet_name)
                df.to_excel(writer, sheet_name="All", index=False, header=False, startrow=writer.sheets['All'].max_row)


# Call the following lines on the command line to create the GO and pathway enrichment excel files
# cd ~/Shiu_Lab/Project/Scripts/Data_Vis/
# python 0_SuppDatFile_14_15_combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/GO_Enrichment/SNPs_fs ORA_Genes_[A-Z0-9]+_ SNP go Scripts/Data_Vis/Section_4
# python 0_SuppDatFile_14_15_combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/GO_Enrichment/PAVs_fs ORA_Genes_[A-Z0-9]+_ PAV go Scripts/Data_Vis/Section_4
# python 0_SuppDatFile_14_15_combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/GO_Enrichment/CNVs_fs ORA_Genes_[A-Z0-9]+_ CNV go Scripts/Data_Vis/Section_4
# python 0_SuppDatFile_14_15_combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/PWY_Enrichment/SNPs_fs ORA_PWY_Genes_[A-Z0-9]+_ SNP pwy Scripts/Data_Vis/Section_4
# python 0_SuppDatFile_14_15_combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/PWY_Enrichment/PAVs_fs ORA_PWY_Genes_[A-Z0-9]+_ PAV pwy Scripts/Data_Vis/Section_4
# python 0_SuppDatFile_14_15_combine_go_pwy_enrichment.py Scripts/Genomic_Prediction_RF/PWY_Enrichment/CNVs_fs ORA_PWY_Genes_[A-Z0-9]+_ CNV pwy Scripts/Data_Vis/Section_4
