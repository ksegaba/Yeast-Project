#!/usr/bin/env python3
import pandas as pd

# Determine which models are missing from the RF results summary file
def get_missing(dirpath):
    res = pd.read_csv(f'{dirpath}/RESULTS_reg.txt', sep='\t') # RF results summary file
    
    # get the list of feature selection files that did not run
    paths = {}
    lit = ['caffeine', 'caffeine', 'sodium_meta-arsenite', 'benomyl', 'cuso4']
    envs = ['YPDCAFEIN40', 'YPDCAFEIN50', 'YPDSODIUMMETAARSENITE', 'YPDBENOMYL500', 'YPDCUSO410MM']
    
    for dat in ['snp', 'pav', 'cnv']:
        for i, env in enumerate(envs):
            for j in range(100):
                if 'top_non_lit' in dirpath:
                    prefix = f'{dat}_{env}_top_non_{lit[i]}_genes_randomized_{j}'
                else:
                    prefix = f'{dat}_{env}_{lit[i]}_lit_genes_randomized_{j}'
                paths[prefix] = (prefix in res.ID.values)
    
    len(paths.keys()) # sanity check should be 1500, passed
    
    with open(f'{dirpath}_missing_runs.txt', 'w') as f:
        for key in paths.keys():
            if not paths[key]:
                f.write(f'{key}\n')


get_missing('/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/only_lit_genes_models')
get_missing('/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/randomized/only_top_non_lit_genes_models')

