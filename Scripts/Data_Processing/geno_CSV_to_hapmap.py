#!/usr/bin/env python3
"""
Convert imputed biallelic genotypes CSV to Hapmap format

Note: The homozygous major allele genotype is encoded as -1, heterozygous as 0,
    and homozygous minor allele as 1.

Example:
    path=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
    file=geno.csv
    vcf=1011Matrix_diploid750_biallelic_maf0.05_11032022.recode.vcf
    save=geno.csv.hmp.txt
    python geno_CSV_to_hapmap.py -path $path -file $file -vcf $vcf -save $save

Arguments:
    path (str): path to working directory.
    file (str): name of the imputed biallelic genotypes CSV in working directory.
    vcf (str): name VCFv4.1 file in working directory.
    save (str): name to save Hapmap file as in working directory.

Returns:
    Saves a Hapmap file in the working directory.
"""

__author__ = "Kenia Segura Abá"

import sys, os
import argparse
import warnings
import datatable as dt
import pandas as pd
import numpy as np
from tqdm import tqdm


def warn(*args, **kwargs):
    pass
warnings.warn = warn


def geno_CSV_to_hmp(df, hmp):
    """ Convert fastPHASE imputed biallelic genotypes file to Hapmap """
    for row in tqdm(range(df.shape[0])):
        tmp = np.where(df.iloc[row,:]==-1, vcf["REF"]+vcf["REF"], \
            np.where(df.iloc[row,:]==0, vcf["REF"]+vcf["ALT"], \
            np.where(df.iloc[row, :]==1, vcf["ALT"]+vcf["ALT"], \
            "NA")))
        # check for no NAs
        if sum(tmp==np.nan)==0:
            hmp[df.index[row]] = tmp # add column to hmp

    return (0)


if __name__=="__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Test script to convert vcf to matrix format")
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-path", type=str, help="Path to working directory.", required=True)
    req_group.add_argument("-file", type=str, help="Name of fastPHASE imputed _genotype.out file in the working directory.", required=True)
    req_group.add_argument("-vcf", type=str, help="Name of VCFv4.1 file in the working directory.", required=True)
    req_group.add_argument("-save", type=str, help="Name to save Hapmap file as.", required=True)
    if len(sys.argv)==1:
        parser.print_help()
        sys.exist(0)
    args = parser.parse_args()

    # Set working directory
    os.chdir(args.path)
    
    # Read in fastPHASE imputed genotypes file & PLINK mapping file
    df = dt.fread(args.file) # CSV file
    df = df.to_pandas()
    df.set_index(df.columns[0], inplace=True)
    vcf = dt.fread(args.vcf, skip_to_line=59) # VCF v4.1
    vcf = vcf[:,0:5] # subset vcf columns
    vcf = vcf.to_pandas() # convert to pandas
    
    # Create Hapmap dataframe
    cols = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
            "protLSID", "assayLSID", "panel", "QCcode"]
    hmp = pd.DataFrame(columns=cols, index=range(0,vcf.shape[0]))
    hmp.set_index(vcf["#CHROM"].str.cat(vcf["POS"].astype(str), sep="_"), inplace=True) # set index
    hmp["chrom"] = vcf["#CHROM"].str.extract(r"([0-9]+)$", expand=True).values # chromosome number
    hmp["pos"] = vcf["POS"].values # SNP position
    hmp["strand"] = "+" # designate as the 5' end of the DNA strand
    hmp["alleles"] = vcf["REF"].str.cat(vcf["ALT"], sep="/").values # join alleles into one column
    
    # Add genotypes to Hapmap dataframe
    geno_CSV_to_hmp(df=df, hmp=hmp) # convert to CSV

    # Save hmp to file
    hmp.to_csv(args.save, sep="\t", index=False, chunksize=1000)
    # hmp = dt.Frame(hmp)
    # hmp.to_csv(args.save, sep="\t", quoting="none") # sep not working
