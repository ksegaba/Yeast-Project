#!/usr/bin/env python3

"""
Convert fastPHASE v1.4.8 imputed biallelic genotypes file to CSV matrix format

Note: The minor allele of the fastPHASE file is encoded as 1 and 
    the major allele is encoded as 2.

Example:
    path=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/
    file=1011Matrix_diploid750_biallelic_maf0.05_11032022_genotypes.out
    map=1011Matrix_diploid750_biallelic_maf0.05_11032022.plink.map
    save=1011Matrix_diploid750_biallelic_maf0.05_11032022_genotypes.csv
    python fastPHASE_to_csv.py -path $path -file $file -map $map -save $save
    
Args:
    path (str): Path to working directory.
    file (str): Name of fastPHASE imputed _genotype.out file.
    map (str): Name of PLINK map file containing chromosome and bp position of genotypes.
    save (str): Name to save matrix file as.

Returns:
    Saves a comma-separate values file in the working directory. Genotypes are
    encoded as -1 (homozygous major), 0 (heterozygous), 1 (homozygous minor)[1].

Reference:
    [1] Kim, S., & Lee, C. Y. (2020). A consistent approach to the genotype 
        encoding problem in a genome-wide association study of continuous 
        phenotypes. PloS one, 15(7), e0236139.
"""

__author__ = "Kenia Segura Abá"

import sys, os
import argparse
import warnings
import pandas as pd
import datatable as dt
from tqdm import tqdm
import pyspark


def warn(*args, **kwargs):
    pass
warnings.warn = warn


def fastPHASE_to_csv(inf, map, save):
    """ Convert fastPHASE imputed biallelic genotypes file to CSV """
    
    # read line by line
    inl = inf.readline() # read the first line of the file
    
    geno = {} # dict of IDs and recoded genotypes
    
    with tqdm(total=751) as pbar: # progress bar
        while inl.strip() != "END GENOTYPES":
            if inl.strip() == "BEGIN GENOTYPES":
                pass

            if inl.strip().startswith("#"):
                # fetch individual IDs
                ID = inl.strip().split(" ")[2] # individual ID
                if ID not in geno.keys(): # check if ID exists
                    geno[ID] = {}

                # fetch genotypes for all SNPs
                inl = inf.readline().strip().split(" ") # allele 1
                inl2 = inf.readline().strip().split(" ") # allele 2

                # recode SNP genotype
                for i in range(len(inl)):
                    snp = map.iloc[i,1] # SNP name
                    snp = snp.replace(":", "_") # change to chr_pos format
                    if snp not in geno[ID].keys():
                        if inl[i]=="2" and inl2[i]=="2": # homozygous major
                            geno[ID][snp] = -1
                        elif inl[i]!=inl2[i]: # heterozygous
                            geno[ID][snp] = 0
                        elif inl[i]=="1" and inl2[i]=="1": # homozygous minor
                            geno[ID][snp] = 1
            inl = inf.readline()
            pbar.update(1)
    pbar.close()

    # save recoded genotype matrix to file
    geno = pd.DataFrame(geno)
    geno = geno.T
    geno.rename(columns={'index':'ID'}, inplace=True)
    geno.reset_index(inplace=True)
    geno = dt.Frame(geno)
    geno.to_csv(save, quoting="none")
    return(0)


if __name__=="__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Test script to convert vcf to matrix format")
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-path", type=str, help="Path to working directory.", required=True)
    req_group.add_argument("-file", type=str, help="Name of fastPHASE imputed _genotype.out file in the working directory.", required=True)
    req_group.add_argument("-map", type=str, help="Name of PLINK map file containing chromosome and bp position of genotypes.")
    req_group.add_argument("-save", type=str, help="Name to save matrix file as.", required=True)
    if len(sys.argv)==1:
        parser.print_help()
        sys.exist(0)
    args = parser.parse_args()

    # Set working directory
    os.chdir(args.path)
    
    # Read in fastPHASE imputed genotypes file & PLINK mapping file
    inf = open(args.file, "r")
    map = pd.read_csv(args.map, sep="\t", header=None)
    fastPHASE_to_csv(inf=inf, map=map, save=args.save) # convert to CSV
    inf.close()