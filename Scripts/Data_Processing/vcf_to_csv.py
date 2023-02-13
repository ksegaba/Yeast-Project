"""
Convert VCF file to CSV format

Args:
    vcf (str): path to VCF file
    out (str): path to CSV file to save as
"""
__author__ = "Kenia E. Segura Abá"

import sys
import argparse
import warnings
import re
import pandas_plink

def warn(*args, **kwargs):
    pass
warnings.warn = warn

if __name__=="__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Convert VCF to CSV format")
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-vcf", type=str, help="path to VCF file")
    req_group.add_argument("-out", type=str, help="path to CSV file to save as")
    if len(sys.argv)==1:
        parser.print_help()
        sys.exist(0)
    args = parser.parse_args()
    
    # Convert to CSV