"""
Cross-environment transfer learning with fitted random forest models.

Comparisons:
- YPDBENOMYL500 to all the others
- YPDCAFEIN40 to all the others
- YPDCUSO410MM to all the others

Command:
dir1 = /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
dir2 = /mnt/gs21/scratch/seguraab/yeast_project/yeast_rf_results
Cross_env_transfer.py -dir1 $dir1 -X geno.csv -dir2 $dir2

Author: Kenia Segura Abá
"""

import sys, os, argparse
import warnings
import joblib
import datatable as dt
import pandas as pd
import numpy as np
from itertools import permutations

def warn(*args, **kwargs):
    pass

warnings.warn = warn

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Cross-environment transfer')
    # Required input
    req_group = parser.add_argument_group(title='Required Input')
    req_group.add_argument('-dir1', help='path to data directory', required=True)
    req_group.add_argument('-X', help='input X data', required=True)
    #req_group.add_argument('-y', help='y data', required=True)
    #req_group.add_argument('-test', help='test instances', required=True)
    req_group.add_argument('dir2', help='path to RF output', required=True)

    #req_group.add_argument('-model', help='', required=True)
    
   


    # Read in data
    dir = "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018"
    X = dt.fread("%s/geno.csv" % dir) # Genotype data
    X = X.to_pandas()
    X.set_index(X.columns[0], inplace=True)
    y = pd.read_csv("%s/pheno.csv" % dir, index_col=0) # Phenotype data
    test = pd.read_csv("%s/Test.txt" % dir, header=None) # Test instances

    # All possible environment pairs
    # permut = permutations(Y.columns, 2)
    # p = [i for i in permut]
    
    # RF results
    dir = "/mnt/gs21/scratch/seguraab/yeast_project/yeast_rf_results"
    results = pd.read_csv("%s/RESULTS_reg.txt" % dir, sep="\t")

    # Optimal # of features for each single-environment snp model
    snp_feats = {'YPACETATE':16, 'YPDCAFEIN40':512, 'YPDHU':128, 
                'YPETHANOL':1024,'YPD14':256, 'YPDCAFEIN50':1024, 
                'YPDKCL2M':256, 'YPGALACTOSE':128, 'YPD40':1024, 
                'YPDCHX05':1024, 'YPDLICL250MM':128, 'YPGLYCEROL':32, 
                'YPD42':512, 'YPDCHX1':512, 'YPDMV':256, 'YPRIBOSE':128, 
                'YPD6AU':1024, 'YPDCUSO410MM':1024, 'YPDNACL15M':512, 
                'YPSORBITOL':256, 'YPDANISO10':128, 'YPDDMSO':256, 
                'YPDNACL1M':1024, 'YPXYLOSE':256, 'YPDANISO20':1024, 
                'YPDETOH':32, 'YPDNYSTATIN': 16, 'YPDANISO50': 256, 
                'YPDFLUCONAZOLE':256, 'YPDSDS':1024, 'YPDBENOMYL200':512, 
                'YPDFORMAMIDE4':256, 'YPDSODIUMMETAARSENITE':256, 
                'YPDBENOMYL500':1024, 'YPDFORMAMIDE5':128}
    orf_feats =
    cno_feats = 

    # Environment pairs
    envs_A = ['YPDBENOMYL500', 'YPDCAFEIN40', 'YPDCUSO410MM']
    envs_B = Y.columns.to_list()

    A = 'YPDCAFEIN40'
    B = 'YPD40'
    
    # data frame to collect transfer results
    out = pd.DataFrame(columns=['Env A', 'Env B', 'R2', 'PCC']) 
    for A in envs_A:
        # Load model fitted to environment A
        mod_A = joblib.load('%s/%s_exp_rf_%i_models.pkl' % (dir, A, snp_feats[A]))
        
        # Features of environment A
        fs_A = pd.read_csv('%s/feat_exp_rf_%s_top_%i' % (dir, A, snp_feats[A]), header=None)
        
        for B in envs_B:
            if A != B:
                # Features of environment B
                fs_B = pd.read_csv('%s/feat_exp_rf_%s_top_%i' % (dir, B, snp_feats[B]), header=None)
                overlap = set(fs_A) & set(fs_B) # check if top features match

                if len(overlap) == 0: # feature sets differ

                    # Testing set of environment B
                    X_test = X.loc[X.index.isin(test[0])]
                    X_test = X_test.loc[:,fs_A[0]]

                    # Apply transfer to environment B
                    y_pred = mod_A.predict(X_test)

                    # Compare results to model fitted to environment B
                    y_test = y.loc[y.index.isin(test[0])]
                    cor = np.corrcoef(np.array(y_test[B]), y_pred)

if __name__ == "__main__":
    main()