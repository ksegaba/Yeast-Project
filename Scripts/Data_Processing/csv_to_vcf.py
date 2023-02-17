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
        # write header (from unfiltered vcf file used to generate csv file)
        vcf.write("##fileformat=VCFv4.1\n")
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write("##contig=<ID=1,length=230218>\n")
        vcf.write("##contig=<ID=2,length=813184>\n")
        vcf.write("##contig=<ID=3,length=316620>\n")
        vcf.write("##contig=<ID=4,length=1531933>\n")
        vcf.write("##contig=<ID=5,length=576874>\n")
        vcf.write("##contig=<ID=6,length=270161>\n")
        vcf.write("##contig=<ID=7,length=1090940>\n")
        vcf.write("##contig=<ID=8,length=562643>\n")
        vcf.write("##contig=<ID=9,length=439888>\n")
        vcf.write("##contig=<ID=10,length=745751>\n")
        vcf.write("##contig=<ID=11,length=666816>\n")
        vcf.write("##contig=<ID=12,length=1078177>\n")
        vcf.write("##contig=<ID=13,length=924431>\n")
        vcf.write("##contig=<ID=14,length=784333>\n")
        vcf.write("##contig=<ID=15,length=1091291>\n")
        vcf.write("##contig=<ID=16,length=948066>\n")
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
            vcf.write(f"{chr.lstrip('chromosome')}\t{pos}\t{snp}\t.\t.\t.\t.\t.\tGT\t{genotypes}\n")

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
    csv2vcf(csv, args.out)
