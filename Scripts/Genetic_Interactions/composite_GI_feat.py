
"""
Using biallelic SNPs as features, make composite SNP-SNP features that represent
genetic interactions (GI) based on the yeast global genetic interaction network
(Costanzo et al. 2016)

Kenia Segura Abá 

Command to run on 4 cores: python encode_GI_feat.py
"""
import os
import time
import datatable as dt
import pandas as pd
import numpy as np
from alive_progress import alive_bar
from multiprocessing import Pool, cpu_count

def make_GI_feats(geno, genes, genes2, All, dir, save):
    with alive_bar(len(geno.columns), bar = "circles", spinner = "dots_waves") as bar: # progress bar
        # Map SNPs to GIs
        names = [] # list containing GI feature names
        dict = {} # dictionary containing GI features to convert to dataframe
        for snp in geno.columns:
            if snp == "ID": # skip isolate ID label in geno
                #df_GI.append(geno["ID"]) # add isolate IDs column
                dict["ID"] = geno["ID"]
                names.append("ID")
                pass
            else:
                # first find corresponding gene
                if snp in genes:
                    gene = genes[snp][3] # gene in which snp is in
                    if gene == "intergenic": # exclude snps labelled as intergenic
                        pass
                    else:
                        # then find GI(s)
                        if All["Query Strain ID"].str.contains(gene).any(): # check if gene has mutant query strain
                            GIs = All[All["Query Strain ID"].str.contains(gene)]
                            # make name for composite GI features
                            for ID in GIs["Array Strain ID"].str.split("_"):
                                gene2 = ID[0]
                                snp2 = []
                                # find corresponding snp(s) to array gene
                                if gene2 in genes2:
                                    vals = genes2[gene2]
                                    for s in vals:
                                        snp2.append(s[0])
                                else: # no snp info
                                    pass
                                # combine snp and snp2 columns into a composite feature
                                for s2 in snp2:
                                    name = "%s_%s-%s_%s"%(gene,snp,gene2,s2) # name of feature
                                    feat = geno[snp]+"_"+geno[s2] # pandas series
                                    #print("geno[snp] length: ", geno[snp].shape) # 750
                                    #print("geno[s2] length: ", geno[s2].shape) # 750
                                    #print("Feature length: ", feat.size) # 750
                                    dict[name] = feat
                                    names.append(name) 

                else:
                    pass
            bar()

    # one-hot encode df_GI?
    return names, dict

def list_to_file(df_GI): # this is not outputting like I want it too
    dir = "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018"
    with open("%s/geno_GI_All.csv"%dir, "w") as out:
        for col in df_GI:
            for elem in col:
                out.write(elem+",")
            out.write("\n")
    out.close()

def main():
    
    dir = "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018"
    dir2 = "/mnt/home/seguraab/Shiu_Lab/Project/Data/Costanzo_S1"
    
    print("Reading data...")
    # Read in data
    geno = dt.fread("%s/geno.csv"%dir) # genotyp matrix (SNP features)
    geno = geno.to_pandas() # convert to pandas dataframe
    geno = geno.applymap(str) # convert int32 to str types

    genes = pd.read_csv("%s/biallelic_snps_diploid_and_S288C_genes.txt"%dir, header=None, index_col=0) # SNPs mapped to S288C genes
    genes2 = genes.copy() # copy into a new dataframe
    genes = genes.to_dict("index") # convert to dictionary
    genes2 = genes2.reset_index() # remove snp IDs as index
    col_names = [3,0,1,2] # new column order
    genes2 = genes2.reindex(columns = col_names) # reorder columns
    genes2 = genes2.set_index(3) # set index to gene IDs
    genes2 = genes2.groupby(level=0).apply(lambda x: x.to_dict("records")).to_dict() # convert to dictionary

    #ExE = pd.read_csv("%s/SGA_ExE_1per.txt"%dir2, header=0) # Essential vs Essential GIs (1%)
    #ExN = pd.read_csv("%s/SGA_ExN_1per.txt"%dir2, header=0) # Essential vs Nonessential GIs (1%)
    #NxN = pd.read_csv("%s/SGA_NxN_1per.txt"%dir2, header=0) # Nonessential vs Nonessential GIs (1%)
    All = pd.read_csv("%s/SGA_combined_1per.txt"%dir2, header=0) # The All GIs (1%) *Not the same as the above combined!

    # Execute
    print("Composing features...")
    names, dict = make_GI_feats(geno, genes, genes2, All, dir, "geno_GI_All.csv")
    print("dict", len(dict))
    #print(dict)    

    del geno # to prevent bus error
    print("Creating dataframe from dict...")
    # pass dict or numpy array to pd.DataFrame
    start = time.perf_counter()
    d1 = pd.DataFrame(data=dict, index=None) # I need to check if this is what's losing several isolate rows
    print("Final dataframe shape:", d1.shape)
    end = time.perf_counter()
    print(f"Elapsed Time: {end-start} seconds")

    print("Saving features to file...")
    # Save new input dataframe to use in models
    start = time.perf_counter()
    d1.to_csv("%s/geno_GI_All.csv"%dir) # or else this took way too long and timed out and that's why I'm missing rows.
    end = time.perf_counter()
    print(f"Elapsed Time: {end-start} seconds")

    # Save list to file in parallel
    #p = Pool(cpu_count())
    #p.map(list_to_file, (df_GI))
    #p.close()
    #p.join()

if __name__ == "__main__":
    main()