'''
filter genotype matrix/ORF copy number matrix to only include diploid isolates
Input: genotype matrix/ORF copy number matrix and list/matrix of labels
Author: Kenia Segura Abá
'''
import sys,os,argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='This code is for converting the genotype gvcf to matrix format')
    # Required
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument('-file', help='genotype matrix file', required=True)
    req_group.add_argument('-labels', help='labels file containing individual/isolate ploidy information', required=True)
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    inf = pd.read_csv(args.file, sep='\t') 
    print("input", inf.shape)
    labels = pd.read_csv(args.labels, sep='\t') # ex. Isolate_ploidy.tsv
    print("labels", labels.shape)

    # select on diploid isolates
    diploid = labels.loc[labels['Ploidy'] == '2']
    print("diploid", diploid.shape)

    # filter inf to extract only diploid
    if 'Chr' in inf.columns: # genotype matrix
        temp = inf[['Chr','Pos','Ref','Alt']]
        inf_diploid = inf[inf.columns.intersection(diploid.Isolate)]
        print("temp", temp.shape, "\ninf_diploid", inf_diploid.shape)
        # Merge dataframes by index
        out = temp.merge(inf_diploid, left_index=True, right_index=True)
        print("out", out.shape)
    else: # ORF copy number matrix
        out = inf[inf.STD_name.isin(diploid.Isolate)]
        print("out", out.shape)
    
    # save subset to a file
    out.to_csv(args.file+'_diploid.txt', header=True, sep='\t', index=False)

if __name__ == '__main__':
	main()
