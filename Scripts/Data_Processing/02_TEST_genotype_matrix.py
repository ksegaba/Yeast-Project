#!/usr/bin/env python3

"""
Test script to convert a filtered VCF file containing only biallelic sites 
to matrix format and save it to a file. This code corresponds to 
02_genotype_matrix.py written by Peipei Wang.

Example:
    $ path=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/
    $ vcf=1011Matrix_diploid750_biallelic_maf0.05_10172022.recode.vcf
    $ save=1011Matrix_diploid750_biallelic_maf0.05_10172022.tsv
    $ python 02_TEST_genotype_matrix.py -path $path -vcf $vcf -save $save

Args:
    path (str): Path to working directory.
    vcf (str): Name of VCFv4.1 file in the working directory.
    save (str): Name to save matrix file as.

Returns:
    Saves a tab-separated values file in the working directory.
"""
__author__ = "Kenia Segura Abá"

import sys, os
import argparse
import warnings
import pandas as pd
import numpy as np


def warn(*args, **kwargs):
    pass
warnings.warn = warn


def test_01_genotype_matrix(vcf, save):
    """ Test for 02_genotype_matrix.py by Peipei Wang """
    # outfile to save matrix
    out = open(save, "w")
    
    # check: first line matches to "##fileformat=VCFv4.1"
    inl = vcf.readline() # read the first line of the file
    if inl.strip()=="##fileformat=VCFv4.1": 
        print("Passed check: first line matches to '##fileformat=VCFv4.1'")
    
    # check: read line by line, does it begin with ## or #CHROM
    while inl:
        if inl.startswith("##"):
            pass
        elif inl.startswith("#CHROM"):
            print(inl)
            out.write(inl[1:])
            out.flush()
            print("Passed check: read line by line, line begins with #CHROM")
        else:
            break
        inl = vcf.readline()

    # check biallelic snps
    while inl:
        tmp = inl.split("\t") # columns are tab-delimeted
        ref = tmp[3] # reference allele
        alt = tmp[4] # alternate allele
        # check: ALT (index 4) is not multiallelic (comma-delimeted)
        if len(alt.split(",")) >= 2: # 2 or more alternate alleles
            print("Failed: ALT is multiallelic")
            inl = vcf.readline()
            continue
        # check: FORMAT (index 8) is GT (vcf v4.1) and grab corresponding genotypes
        elif tmp[8].split(":")[0]=="GT": # data types in FORMAT are colon-delimeted
            print("Passed: ALT is only one allele")
            print("Passed: 1st element in FORMAT is GT")
            geno = tmp[9:] # genotypes of individuals start from index 9 onward
            geno = [g.split(":")[0] for g in geno] # genotypes only
            
            # check encoding
            # reference allele encoded as 0 in genotypes
            # alternate allele encoded as 1 in genotypes
            geno = [g.replace('0',ref) for g in geno]
            geno = [g.replace('1',alt) for g in geno]

            # check no alleles were lost
            if len(geno)==len(tmp[9:]):
                outl = tmp[:5] + geno
                print("Passed: ref and alt alleles are encoded correctly")

            # write outl to file
            out.write("\t".join(outl)+"\n")
            out.flush()

        inl = vcf.readline()
    out.close()

    return 0


if __name__=="__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Test script to convert vcf to matrix format")
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-path", type=str, help="Path to working directory.", required=True)
    req_group.add_argument("-vcf", type=str, help="Name of VCFv4.1 file in the working directory.", required=True)
    req_group.add_argument("-save", type=str, help="Name to save matrix file as.", required=True)
    if len(sys.argv)==1:
        parser.print_help()
        sys.exist(0)
    args = parser.parse_args()

    # Set working directory
    os.chdir(args.path)
    
    # Read in VCF file
    vcf = open(args.vcf, "r")
    test_01_genotype_matrix(vcf=vcf, save=args.save) # run check
    vcf.close()