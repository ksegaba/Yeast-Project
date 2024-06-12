"""
Convert CSV genotype file to Variant Call Format (VCF)

The marker IDs (columns) in the CSV are in the <chromosome>_<position> format
and genotypes are encoded as -1 (homozygous reference allele), 0 (heterozygous), 
1 (homozygous alternate allele). The CSV file contains a subset of the markers
found in the original VCF file. The original VCF file is used to extract the header
for the new VCF file to be made.

Args:
    gvcf (str): path to original VCF file from which the csv was generated
    csv (str): path to CSV file of genotypes
    out (str): path to save VCF file as
"""
__author__ = "Kenia E. Segura Abá"

import sys
import argparse
import warnings
import datatable as dt
from tqdm import tqdm


def warn(*args, **kwargs):
    pass
warnings.warn = warn


def csv2vcf(csv, gvcf, save_path):
    """Convert CSV file to VCF format"""
    with open(save_path, "w") as vcf:
        print(f"Writing header to {save_path} and getting SNP information from original VCF...")
        snp_info = {}
        for i,line in enumerate(tqdm(gvcf)):
            if line.startswith("##"): # write header lines
                vcf.write(f"{line.strip()}\n")
            else:
                # save ID, REF, ALT, QUAL, FILTER, INFO, FORMAT for snps in csv
                if "_".join(line.strip().split("\t")[0:2]) in csv.names[1:]:
                    snp_info["_".join(line.strip().split("\t")[0:2])] = "\t".join(line.strip().split("\t")[2:9])
                elif len(snp_info)==len(csv.names[1:]):
                    break
                else:
                    pass
        cols="\t".join(map(str, csv[:,0].to_list()[0])) # individual IDs
        vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{cols}\n")
        print(f"Writing marker information to {save_path}...")
        for i in tqdm(range(len(csv.names[1:]))):
            snp = csv.names[1:][i]
            chr = snp.split("_")[0][10:] # to exclude 'chromosome' from name
            pos = snp.split("_")[1]
            # convert marker genotypes to strings
            col = [f"{csv[i, snp]}" for i in range(csv.shape[0])] # snp column
            col = ["0/0" if x=="-1" else x for x in col]
            col = ["0/1" if x=="0" else x for x in col]
            col = ["1/1" if x=="1" else x for x in col]
            genotypes="\t".join(map(str, col))
            # write the snp to the new VCF
            vcf.write(f"{chr}\t{pos}\t{snp_info[snp]}\t{genotypes}\n")
    print("Done!")

if __name__=="__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Convert VCF to CSV format")
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-gvcf", type=str, help="path to original VCF file from which the csv was generated")
    req_group.add_argument("-csv", type=str, help="path to CSV file of genotypes")
    req_group.add_argument("-out", type=str, help="path to save VCF file as")
    if len(sys.argv)==1:
        parser.print_help()
        sys.exist(0)
    args = parser.parse_args()

    # Read in CSV file and gVCF file
    gvcf = open(args.gvcf, "r")
    csv = dt.fread(args.csv)

    # Convert to VCF
    csv2vcf(csv, gvcf, args.out)
    gvcf.close()
