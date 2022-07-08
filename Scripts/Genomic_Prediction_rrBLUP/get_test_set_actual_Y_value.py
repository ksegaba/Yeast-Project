#!/opt/software/Python/3.6.4-foss-2018a/bin/python
"""
Description: Parse phenotype matrix to obtain actual trait values of the test set
Author: Kenia Segura Abá
"""
import sys,os,argparse
import pandas as pd
import warnings

def warn(*args, **kwargs):
    pass
warnings.warn = warn

def main():
    parser = argparse.ArgumentParser(description="This code is for parsing phenotype matrix to extract test set actual trait values")
    # Required
    req_group = parser.add_argument_group(title="Required Input")
    req_group.add_argument("-p", help="phenotype matrix", required=True)
    req_group.add_argument("-t", help="list of individuals in test set", required=True)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # Load in data
    p_file = args.p
    t_file = args.t
    pheno = pd.read_csv(p_file, header=0, index_col=None,sep=",")
    test = pd.read_csv(t_file, header=None, index_col=None)

    # Merge data based on test
    final = pd.merge(test, pheno, left_on=test.columns[0], right_on=pheno.columns[0])
    final.to_csv("pheno_test.csv", header=True, index=True)

if __name__ == '__main__':
	main()