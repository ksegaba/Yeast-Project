"""
Convert CSV genotype file to Variant Call Format (VCF)

The marker IDs (columns) in the CSV are in the <chromosome>_<position> format
and genotypes are encoded as -1 (homozygous reference allele), 0 (heterozygous), 
1 (homozygous alternate allele).

Args:
    csv (str): path to CSV file
    out (str): path to save VCF file as
"""
__author__ = "Kenia E. Segura Abá"

import sys
import argparse
import warnings
import datatable as dt

def warn(*args, **kwargs):
    pass
warnings.warn = warn

def csv2vcf(csv, save_path):
    """Convert CSV file to VCF format"""
    with open(save_path, "w") as vcf:
        # write header
        vcf.write("##fileformat=VCFv4.3\n")
        cols="\t".join(map(str, csv[:,0].to_list()[0])) # individual IDs
        vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{cols}\n")
        # write marker information
        for snp in csv.names[1:]:
            chr = snp.split("_")[0]
            pos = snp.split("_")[1]
            # convert marker genotypes to strings
            col = csv[:, dt.as_type(csv[snp], str)] # snp column
            col.replace({"-1":"0/0", "0":"0/1", "1":"1/1"})
            genotypes="\t".join(map(str, col.to_list()[0]))
            vcf.write(f"{chr}\t{pos}\t.\t.\t.\t.\t.\t.\tGT\t{genotypes}\n")

if __name__=="__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Convert VCF to CSV format")
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-csv", type=str, help="path to CSV file")
    req_group.add_argument("-out", type=str, help="path to save VCF file as")
    if len(sys.argv)==1:
        parser.print_help()
        sys.exist(0)
    args = parser.parse_args()

    # Read in CSV file
    csv = dt.fread(args.csv)

    # Convert to VCF
    csv2vcf(file, args.out)
