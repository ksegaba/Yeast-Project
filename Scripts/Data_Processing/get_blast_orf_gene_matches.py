#!/usr/bin/env python3
"""
Obtain an ORF to gene mapping from BLAST output files.
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project/")

def fill_missing(row, fill_col='organism', lookup_col='gene', lookup_dict={}):
    if pd.isnull(row[fill_col]) and row[lookup_col] in lookup_dict:
        return lookup_dict[row[lookup_col]]
    return row[fill_col]


def get_final(row, lookup_col1='pident_x', lookup_col2='pident_y', fill_col1='gene_x', fill_col2='gene_y', reverse=False):
    if row[lookup_col1] > row[lookup_col2]:
        if reverse:
            return row[fill_col2]
        else:
            return row[fill_col1] 
    elif row[lookup_col1]==row[lookup_col2]:
        if row[fill_col1]==row[fill_col2]:
            return row[fill_col1]
        else:
            return " // ".join([row[fill_col1], row[fill_col2]])
    else:
        return row[fill_col2]


if __name__ == "__main__":
    ########################## INITIAL BLASTX RESULTS ##########################
    #### Filter ORF features mapped to S288C genes file
    s288c_blastx = pd.read_csv("Data/S288C_reference_genome_R64-3-1_20210421/S288C_orf_peter_blastx.txt", skiprows=1, sep="\t")
    strict = s288c_blastx.groupby("qacc").apply(lambda group:
        group[(group.evalue == group.evalue.min()) & (group["pident"] >= 95) &
            (group.pident == group.pident.max())]).reset_index(drop=True)
    strict["organism"] = "Saccharomyces cerevisiae S288C"

    #### Filter ORF features mapped to Arabidopsis thaliana genes
    tair10_blastx = pd.read_csv("Data/Arabidopsis_Genome_TAIR10.1/peter_orf_TAIR10_blastx.txt", sep="\t")
    tair10_strict = tair10_blastx.groupby("qacc").apply(lambda group:
        group[(group.evalue == group.evalue.min()) & (group["pident"] >= 95) &
            (group.pident == group.pident.max())]).reset_index(drop=True)
    tair10_strict["organism"] = "Arabidopsis thaliana TAIR10.1"

    #### Filter ORF features mapped to Drosophila melanogaster genes
    iso1mt_blastx = pd.read_csv("Data/Drosophila_Genome_R6_ISO1MT/peter_orf_ISO1MT_blastx.txt", sep="\t")
    iso1mt_strict = iso1mt_blastx.groupby("qacc").apply(lambda group:
        group[(group.evalue == group.evalue.min()) & (group["pident"] >= 95) &
            (group.pident == group.pident.max())]).reset_index(drop=True)
    iso1mt_strict["organism"] = "Drosophila melanogaster Release 6 plus ISO1 MT"

    #### Filter ORF features mapped to Human genes
    human_blastx = pd.read_csv("Data/Human_Genome_GRCh38.p14/peter_orf_human_blastx.txt", sep="\t")
    human_strict = human_blastx.groupby("qacc").apply(lambda group:
        group[(group.evalue == group.evalue.min()) & (group["pident"] >= 95) &
            (group.pident == group.pident.max())]).reset_index(drop=True) # no matches returned

    #### Filter ORF features mapped to Neurospora crassa genes
    nc12_blastx = pd.read_csv("Data/Neurospora_OR74A_Genome_NC12/peter_orf_NC12_blastx.txt", sep="\t")
    nc12_strict = nc12_blastx.groupby("qacc").apply(lambda group:
        group[(group.evalue == group.evalue.min()) & (group["pident"] >= 95) &
            (group.pident == group.pident.max())]).reset_index(drop=True)
    nc12_strict["organism"] = "Neurospora crassa NC12"

    #### Filter ORF features mapped to non-redundant database
    nr_blastx = pd.read_csv("Data/BLAST_nr_db/peter_orf_nr_blastx.txt", sep="\t")
    nr_strict = nr_blastx.groupby("qacc").apply(lambda group:
        group[(group.evalue == group.evalue.min()) & (group["pident"] >= 95) &
            (group.pident == group.pident.max())]).reset_index(drop=True)
    # pd.Series(nr_strict["sacc"].unique()).to_csv("Data/BLAST_nr_db/nr_genes_for_sgd_08222023.txt", index=False, header=False)

    #### Map nr database identifiers with gene systematic names
    # Note 03/16/2025: I uploaded nr_genes_for_sgd_08222023.txt to NCBI datasets
    # and SGD gene lists tool (called AllianceMine) to obtain gene systematic names.
    ncbi_matches = pd.read_csv("Data/BLAST_nr_db/nr_genes_for_sgd_08222023_ncbi_datasets_matches.tsv", sep="\t", index_col=False)
    sgd_matches = pd.read_csv("Data/BLAST_nr_db/nr_genes_for_sgd_08222023_sgd_matches.tsv", sep="\t", index_col=False)
    
    # unify ncbi and sgd S. cerevisiae matches
    sgd_matches["Nomenclature ID"] = sgd_matches.apply(lambda x: x["Gene.primaryIdentifier"].split(":")[1], axis=1)
    ncbi_matches = ncbi_matches.merge(sgd_matches, how="left", on="Nomenclature ID")
    ncbi_matches["Gene.systematicName"] = ncbi_matches.apply(lambda x: x["Gene.secondaryIdentifier"] if pd.notna(x["Gene.secondaryIdentifier"]) else x["Symbol"], axis=1)
    '''sanity check
    >>> ncbi_matches["Taxonomic Name"].unique()
    array(['Saccharomyces cerevisiae S288C', 'Saccharomyces paradoxus',
       'Saccharomyces eubayanus', 'Zygosaccharomyces rouxii',
       'Nakaseomyces glabratus', 'Vanderwaltozyma polyspora DSM 70294',
       'Naumovozyma dairenensis CBS 421',
       'Tetrapisispora phaffii CBS 4417', 'Torulaspora delbrueckii',
       'Drosophila guanche', 'Naumovozyma castellii'], dtype=object)
    >>> ncbi_matches["Taxonomic Name"].isna().sum()
    0'''
    ncbi_sgd_map = dict(zip(ncbi_matches["Input"], ncbi_matches["Gene.systematicName"])) # nr identifier to gene name map
    gene_org_map = dict(zip(ncbi_matches["Gene.systematicName"], ncbi_matches["Taxonomic Name"])) # gene name to organism map
    
    
    # combine all the blastx results into one file
    matches = pd.concat([strict, tair10_strict, iso1mt_strict, nc12_strict, nr_strict], ignore_index=True)
    matches["Gene.systematicName"] = matches.apply(lambda x: ncbi_sgd_map[x["sacc"]] if x["sacc"] in ncbi_sgd_map.keys() else x["sacc"], axis=1)
    matches["Organism"] = matches.apply(lambda x: gene_org_map[x["Gene.systematicName"]] if x["Gene.systematicName"] in gene_org_map.keys() else x["organism"], axis=1)
    matches["Organism"].unique() # looks good, no repeated org. names
    matches.Organism.isna().sum() # 2588 identifiers with no organism match, they weren't found by ncbi datasets
    matches.loc[matches.Organism.isna(), "Gene.systematicName"].\
        to_csv("Data/BLAST_nr_db/orf_gene_blastx_matches_CORRECTED_missing.tsv", index=False) # I may need to search these manually on ncbi
    matches.to_csv("Data/BLAST_nr_db/orf_gene_blastx_matches_CORRECTED.tsv", sep="\t", index=False)
    
    
    # Get the SGD Identifiers at least, it's more complicated when the organism is different
    '''On the command line:
    mapfile -t query_genes < orf_gene_blastx_matches_CORRECTED_missing.tsv
    echo ${#query_genes[*]} # 2589
    for i in {1..2588}; do
        query="${query_genes[${i}]}"
        /mnt/home/seguraab/ncbi_edirect/efetch -db protein \
            -id ${query_genes[${i}]} -format tabular | grep "/db_xref=\"SGD:" | \
            awk -v q="$query" '{print q, $0}' >> orf_gene_blastx_matches_CORRECTED_missing_output_entrez.tsv
    done
    '''
    matches = pd.read_csv("Data/BLAST_nr_db/orf_gene_blastx_matches_CORRECTED.tsv", sep="\t")
    id_map = pd.read_csv("Data/BLAST_nr_db/orf_gene_blastx_matches_CORRECTED_missing_output.tsv", sep="\t", header=None)
    id_map.drop_duplicates(keep="first", inplace=True)
    id_map.set_index(0, inplace=True)
    to_add = pd.read_csv("Data/BLAST_nr_db/orf_gene_blastx_matches_CORRECTED_missing_output_sgd_results.tsv", sep="\t", index_col=0)
    for i, r in matches.iterrows():
        if r["Gene.systematicName"] in id_map.index:
            print(i)
            matches.iloc[i, 18] = to_add.loc[id_map.loc[r["Gene.systematicName"]].values[0], "Gene.secondaryIdentifier"] # replace "Gene.systematicName"
            matches.iloc[i, 19] = "Saccharomyces cerevisiae S288C" # replace "Organism"
        else:
            continue
    matches.to_csv("Data/BLAST_nr_db/orf_gene_blastx_matches_CORRECTED_final.tsv", sep="\t", index=False)
    
    '''Old code below: 03/16/2025 - cannot delete because I obtained the list of extra organisms to do a reciprocal tblastx on
    ### Combine SGD and NCBI nr_strict gene matches to blastx output
    sgd = pd.read_csv("Data/BLAST_nr_db/nr_genes_sgd_matches_08222023.txt", sep="\t") # nr to sgd results
    sgd_matches = sgd[sgd.reason!="UNRESOLVED"] # genes with no sgd match
    sgd_matches_map = dict(zip(sgd_matches["input"], sgd_matches["secondaryIdentifier"])) # nr blastx sacc to sgd gene IDs
    sgd_gene_org_map = dict(zip(sgd_matches["secondaryIdentifier"], sgd_matches["organism.shortName"])) # sgd gene IDs to organism
    
    ncbi = pd.read_csv("Data/BLAST_nr_db/nr_genes_ncbi_matches_08222023.txt", sep="\t") # nr to ncbi results
    ncbi.rename(columns={"Scientific name":"organism", "Symbol":"gene"}, inplace=True)
    # ncbi["organism"].unique() # organisms from the nr database
    ncbi_matches_map = dict(zip(ncbi["Input"], ncbi["gene"])) # nr blastx sacc to ncbi gene IDs
    ncbi_gene_org_map = dict(zip(ncbi["gene"], ncbi["organism"])) # ncbi gene IDs to organism
    
    matches = pd.concat([strict, tair10_strict, iso1mt_strict, nc12_strict, nr_strict], ignore_index=True) # combine blastx results
    matches["orf"] = matches["qacc"]
    matches["gene"] = matches["sacc"]
    matches.gene = matches.gene.replace(sgd_matches_map) # replace nr blastx sacc with sgd gene IDs
    matches.gene = matches.gene.replace(ncbi_matches_map) # replace nr blastx sacc with ncbi gene IDs
    matches["new_organism"] = matches.apply(fill_missing,args=('organism', 'gene', ncbi_gene_org_map), axis=1) # fill missing organism information
    matches["new_organism"] = matches.apply(fill_missing, args=('new_organism', 'gene', sgd_gene_org_map), axis=1)
    # pd.Series(matches.new_organism.unique()).to_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/nr_blastx_organisms.txt", index=False) # organisms from which all the ORF-gene matches come from
    
    matches[matches.new_organism.isna()]['qacc'].unique() # 1516 ORFs with no organism information
    sum(matches[matches.new_organism.isna()]['sacc'].isin(sgd.input)) # 2447 out of 2468 rows are in the sgd dataframe
    genes2lookup = matches[matches.new_organism.isna()]['sacc']
    # sgd[sgd.input.isin(genes2lookup)]["input"].to_csv("Data/BLAST_nr_db/2375_unresolved_sgd.txt", index=False) # 2375 rows, all are unresolved; none found on neigther ncbi datasets nor sgd databases
    # genes2lookup[~genes2lookup.isin(sgd.input)].to_csv("Data/BLAST_nr_db/19_not_in_sgd.txt", index=False) # 19 genes not in sgd dataframe
    
    sgd_19_genes = pd.read_csv("Data/BLAST_nr_db/19_not_in_sgd.txt", sep="\t") # I manually searched these on ncbi protein database to get the gene and organism information
    sgd_19_genes_map = dict(zip(sgd_19_genes["sacc"], sgd_19_genes["gene"]))
    sgd_19_genes_org_map = dict(zip(sgd_19_genes["gene"], sgd_19_genes["organism"]))
    matches.gene = matches.gene.replace(sgd_19_genes_map) # replace nr blastx sacc with sgd gene IDs
    
    matches.new_organism = matches.apply(fill_missing, args=('new_organism', 'gene', sgd_19_genes_org_map), axis=1)
    len(matches[matches.new_organism.isna()]['qacc'].unique()) # 1508 ORFs with no organism information (some ORFs mapped to resolved genes, hence why 1508+6527!=7708)
    len(matches[~matches.new_organism.isna()]['qacc'].unique()) # 6527 ORFs with organism information
    matches.new_organism.replace(
        {"S. cerevisiae":"Saccharomyces cerevisiae S288C",
         "Saccharomyces paradoxus":"Saccharomyces paradoxus CBS432",
         "Drosophila guanche":"Drosophila guanche DGUA_6",
         "Torulaspora delbrueckii":"Torulaspora delbrueckii CBS 1146 ASM24337v1",
         "Naumovozyma dairenensis CBS 421":"Naumovozyma dairenensis CBS 421 ASM22711v2",
         "Tetrapisispora phaffii CBS 4417":"Tetrapisispora phaffii CBS 4417 ASM23690v1",
         "Zygosaccharomyces rouxii":"Zygosaccharomyces rouxii CBS 732 ASM2636v1",
         "Naumovozyma castellii CBS 4309":"Naumovozyma castellii CBS 4309 ASM23734v1",
         "[Candida] glabrata":"Nakaseomyces glabrata CBS 138 ASM254v2",
         "Vanderwaltozyma polyspora DSM 70294":"Vanderwaltozyma polyspora DSM 70294 ASM15003v1",
         "Saccharomyces eubayanus":"Saccharomyces eubayanus FM1318 SEUB3.0",
         "Apteryx mantelli mantelli":"Apteryx mantelli mantelli AptMant0"}, inplace=True) # replace organism name
    # matches.to_csv("Data/Peter_2018/orf_gene_blastx_matches.txt", sep="\t", index=False) # save to file without removing duplicates
    '''
    
    ######################## RECIPROCAL TBLASTX RESULTS ########################
    ### Read reciprocal tblastx results files
    # blastx = pd.read_csv("Data/Peter_2018/orf_gene_blastx_matches.txt", sep="\t") # initial blastx results; commented out 03/16/2025
    blastx = pd.read_csv("Data/BLAST_nr_db/orf_gene_blastx_matches_CORRECTED_final.tsv", sep="\t")
    s288c_tblastx = pd.read_csv("Data/S288C_reference_genome_R64-3-1_20210421/S288C_to_peter_orf_tblastx.txt", sep="\t")
    s288c_tblastx["organism"] = "Saccharomyces cerevisiae S288C"
    tair10_tblastx = pd.read_csv("Data/Arabidopsis_Genome_TAIR10.1/DNA/TAIR10_to_peter_orf_tblastx.txt", sep="\t")
    tair10_tblastx["organism"] = "Arabidopsis thaliana TAIR10.1"
    iso1mt_tblastx = pd.read_csv("Data/Drosophila_Genome_R6_ISO1MT/DNA/ISO1MT_to_peter_orf_tblastx.txt", sep="\t")
    iso1mt_tblastx["organism"] = "Drosophila melanogaster Release 6 plus ISO1 MT"
    nc12_tblastx = pd.read_csv("Data/Neurospora_OR74A_Genome_NC12/DNA/NC12_to_peter_orf_tblastx.txt", sep="\t")
    nc12_tblastx["organism"] = "Neurospora crassa NC12"
    nt_tblastx = pd.read_csv("Data/BLAST_nr_db/nt_to_peter_orf_tblastx.txt", sep="\t")
    sp_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/S_paradoxus_to_peter_orf_tblastx.txt", sep="\t")
    sp_tblastx["organism"] = "Saccharomyces paradoxus CBS432"
    sp_tblastx["qacc"] = sp_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2])) # replace gene IDs
    dg_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/D_guanche_to_peter_orf_tblastx.txt", sep="\t")
    dg_tblastx["organism"] = "Drosophila guanche DGUA_6"
    dg_tblastx["qacc"] = dg_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    td_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/T_delbrueckii_to_peter_orf_tblastx.txt", sep="\t")
    td_tblastx["organism"] = "Torulaspora delbrueckii CBS 1146 ASM24337v1"
    td_tblastx["qacc"] = td_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    nd_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/N_dairenensis_to_peter_orf_tblastx.txt", sep="\t")
    nd_tblastx["organism"] = "Naumovozyma dairenensis CBS 421 ASM22711v2"
    nd_tblastx["qacc"] = nd_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    tp_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/T_phaffii_to_peter_orf_tblastx.txt", sep="\t")
    tp_tblastx["organism"] = "Tetrapisispora phaffii CBS 4417 ASM23690v1"
    tp_tblastx["qacc"] = tp_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    zr_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/Z_rouxii_to_peter_orf_tblastx.txt", sep="\t")
    zr_tblastx["organism"] = "Zygosaccharomyces rouxii CBS 732 ASM2636v1"
    zr_tblastx["qacc"] = zr_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    nc_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/N_castellii_to_peter_orf_tblastx.txt", sep="\t")
    nc_tblastx["organism"] = "Naumovozyma castellii CBS 4309 ASM23734v1"
    nc_tblastx["qacc"] = nc_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    ng_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/N_glabratus_to_peter_orf_tblastx.txt", sep="\t")
    ng_tblastx["organism"] = "Nakaseomyces glabrata CBS 138 ASM254v2"
    ng_tblastx["qacc"] = ng_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    vp_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/V_polyspora_to_peter_orf_tblastx.txt", sep="\t")
    vp_tblastx["organism"] = "Vanderwaltozyma polyspora DSM 70294 ASM15003v1"
    vp_tblastx["qacc"] = vp_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    se_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/S_eubayanus_to_peter_orf_tblastx.txt", sep="\t")
    se_tblastx["organism"] = "Saccharomyces eubayanus FM1318 SEUB3.0"
    se_tblastx["qacc"] = se_tblastx.qacc.apply(lambda x: "_".join(x.split("_cds_")[1].split("_")[0:2]))
    amm_tblastx =pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/A_mantelli_mantelli_to_peter_orf_tblastx.txt", sep="\t")
    amm_tblastx["organism"] = "Apteryx mantelli mantelli AptMant0"

    ### Gene ID maps
    # pd.Series(pd.concat([nt_tblastx.qacc, sp_tblastx.qacc, dg_tblastx.qacc, td_tblastx.qacc,
    #            nd_tblastx.qacc, tp_tblastx.qacc, zr_tblastx.qacc, nc_tblastx.qacc,
    #            ng_tblastx.qacc, vp_tblastx.qacc, se_tblastx.qacc, amm_tblastx.qacc])\
    # .unique()).to_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_genes_for_ncbi.txt", index=False, header=False) # genes from which all the ORF-gene matches come from
    
    # I saved the ncbi datasets results to tblastx_genes_for_ncbi_results.tsv
    ncbi_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_genes_for_ncbi_results.tsv", sep="\t")
    ncbi_tblastx_map = dict(zip(ncbi_tblastx["Input"], ncbi_tblastx["Symbol"])) # ncbi tblastx qacc to ncbi gene IDs map
    ncbi_tblastx_org_map = dict(zip(ncbi_tblastx["Symbol"], ncbi_tblastx["Scientific name"])) # ncbi tblastx gene IDs to organism map
    
    # Results I got from SGD AllianceMine (03/17/2025)
    # I did this because in the ncbi tblastx file, the S. cerevisiae gene names are the protein names, not the systematic names.
    sgd_tblastx = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_genes_for_ncbi_sgd_results.tsv", sep="\t")
    sgd_tblastx_map = dict(zip(sgd_tblastx["Input"], sgd_tblastx["Gene.systematicName"]))
    sgd_tblastx.organism.replace({"S. cerevisiae":"Saccharomyces cerevisiae",
        "H. sapiens":"Homo sapiens", "D. melanogaster":"Drosophila melanogaster"}, inplace=True) # to match ncbi
    sgd_tblastx_org_map = dict(zip(sgd_tblastx["Gene.systematicName"], sgd_tblastx["organism"]))
    
    ### Subset blastx matches that match the reciprocal tblastx results
    tblastx_combined = pd.concat([s288c_tblastx, tair10_tblastx, iso1mt_tblastx,
                                  nc12_tblastx, nt_tblastx, sp_tblastx, dg_tblastx,
                                  td_tblastx, nd_tblastx, tp_tblastx, zr_tblastx,
                                  nc_tblastx, ng_tblastx, vp_tblastx, se_tblastx,
                                  amm_tblastx], ignore_index=True) # combine results
    tblastx_strict = tblastx_combined.groupby("sacc").apply(lambda group:
        group[(group.evalue == group.evalue.min()) & (group["pident"] >= 95) &
            (group.pident == group.pident.max())]).reset_index(drop=True) # apply strict filter
    # tblastx_strict.qacc.to_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_genes_for_ncbi.tsv", sep="\t", index=False) # added 3/17/2025
    
    # NCBI datasets and SGD results files (03/17/2025)
    ncbi_res = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_genes_for_ncbi_results.txt", sep="\t")
    sgd_res = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_genes_for_ncbi_sgd_results.txt", sep="\t")
    ncbi_res["Nomenclature authority ID"] = ncbi_res.apply(lambda x: x["Symbol"] if not pd.notna(x["Nomenclature authority ID"]) else x["Nomenclature authority ID"], axis=1)
    ncbi_res_map = dict(zip(ncbi_res["Input"], ncbi_res["Nomenclature authority ID"]))
    ncbi_res_org_map = dict(zip(ncbi_res["Nomenclature authority ID"], ncbi_res["Scientific name"]))
    sgd_res["Primary.DBID"] = sgd_res["Primary.DBID"].str.replace("SGD:", "")
    sgd_res_map = dict(zip(sgd_res["Primary.DBID"], sgd_res["Gene.systematicName"]))
    
    # Now, map the tblastx identifiers with the ncbi and sgd identifiers
    tblastx_strict["gene"] = tblastx_strict["qacc"]
    tblastx_strict["orf"] = tblastx_strict["sacc"]
    tblastx_strict.gene = tblastx_strict.gene.replace(ncbi_tblastx_map) # replace nt tblastx qacc with ncbi gene IDs from tblastx results
    tblastx_strict.gene = tblastx_strict.gene.replace(sgd_tblastx_map) # replace nt tblastx qacc with sgd gene IDs from tblastx results
    # tblastx_strict.loc[tblastx_strict.organism.isna(), "gene"].\
    #     to_csv('Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_missing.txt', sep='\t', index=False, header=None) # 3/17/2025
    
    # Get the SGD Identifiers at least, it's more complicated when the organism is different (3/17/2025)
    '''On the command line:
    mapfile -t query_genes < tblastx_strict_missing.txt
    echo ${#query_genes[*]} # 3156
    for i in {0..3155}; do
        query="${query_genes[${i}]}"
        /mnt/home/seguraab/ncbi_edirect/efetch -db nuccore \
            -id ${query_genes[${i}]} -format tabular | \
            awk -v q="$query" '
            /\/gene=/ {gene=$0} 
            /\/organism=/ {organism=$0} 
            END {print q, gene, organism}' >> tblastx_strict_missing_output_entrez.tsv
    done
    '''
    id_map = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_missing_output_entrez.tsv",
        sep="\t", header=None, names=["Input", "Gene", "Organism"])
    id_map.drop_duplicates(keep="first", inplace=True)
    id_map.set_index("Input", inplace=True)
    
    # Entrez retrived the standard gene names, so I had to use SGD to get the systematic names. (3/17/2025)
    id_map_genes = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_missing_output_entrez_sgd_results.tsv",
        sep="\t")
    
    # Map the genes with missing organism info to ncbi and sgd results
    idx_to_check = tblastx_strict.loc[tblastx_strict.organism.isna()].index
    for i, r in tblastx_strict.loc[tblastx_strict.organism.isna()].iterrows():
        try: # get the gene standard name from the tblastx identifier
            if r["gene"] in id_map.index:
                gene = id_map.loc[r["gene"], "Gene"]
            else:
                continue # ncbi entrez did not find a gene match
        except KeyError:
            gene = r["gene"] # the tblastx identifier is already the gene standard name
        
        if gene in id_map_genes.Input.values: # get the gene systematic name
            if pd.isna(id_map_genes.loc[id_map_genes["Input"]==gene, "Gene.systematicName"].values[0]):
                tblastx_strict.iloc[r.name, 13] = id_map_genes.loc[id_map_genes["Input"]==gene, "Primary.DBID"].values[0] # gene column
            else:
                tblastx_strict.iloc[r.name, 13] = id_map_genes.loc[id_map_genes["Input"]==gene, "Gene.systematicName"].values[0]
        else:
            print(f"{gene} not found in ncbi or sgd")
        
        try: # get the organism information
            tblastx_strict.iloc[r.name, 12] = id_map_genes.loc[id_map_genes["Input"]==gene, "Organism"] # organism column
        except ValueError:
            tblastx_strict.iloc[r.name, 12] = id_map.loc[r["gene"], "Organism"]
    
    # There are still 77 genes with no systematic name or organism info (need to get from SGD)
    tblastx_strict.iloc[idx_to_check,:]
    tblastx_strict.loc[tblastx_strict.organism.isna()].qacc.\
        to_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_missing_output_entrez_part2.csv", index=False)
    
    part2_sgd = pd.read_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict_missing_output_entrez_part2_sgd_results.csv", index_col=0)
    
    for i, r in tblastx_strict.loc[tblastx_strict.organism.isna()].iterrows():
        if r["qacc"] in part2_sgd.index:
            tblastx_strict.iloc[i, 13] = part2_sgd.loc[part2_sgd.index==r["qacc"], "Gene.systematicName"].values[0]
            tblastx_strict.iloc[i, 12] = part2_sgd.loc[part2_sgd.index==r["qacc"], "Organism"].values[0]
        else:
            print(r["qacc"])
    
    # Manually add the genes and organism for the 3 remaining identifiers (from ncbi datasets gene)
    tblastx_strict.loc[tblastx_strict.qacc=="NR_132242.1","gene"] = "YNCM0016W"
    tblastx_strict.loc[tblastx_strict.qacc=="NR_132213.1","gene"] = "YNCL0016C"
    tblastx_strict.loc[tblastx_strict.qacc=="NR_132211.1","gene"] = "YNCL0014C"
    tblastx_strict.loc[tblastx_strict.qacc.isin(["NR_132242.1", "NR_132213.1", "NR_132211.1"]), "organism"] = "Saccharomyces cerevisiae S288C"
    
    # Lastly, fix the organism names
    tblastx_strict["organism"] = tblastx_strict["organism"].apply(lambda x: str(x) if isinstance(x, (str, np.ndarray)) else x)
    tblastx_strict["organism"].unique()
    tblastx_strict["organism"].replace("S. cerevisiae", "Saccharomyces cerevisiae", inplace=True)
    tblastx_strict["organism"].unique()
    
    tblastx_strict.to_csv("Data/BLAST_nr_db/nr_blastx_match_organisms/orf_gene_tblastx_matches_CORRECTED_final.tsv", sep="\t", index=False)
    
    # Merge blastx results with tblastx results (match equivalent orf, gene pairs)
    final_matches = blastx[["qacc", "sacc", "evalue", "pident", "Gene.systematicName", "Organism"]].\
        merge(tblastx_strict[["qacc", "sacc", "evalue", "pident", "gene", "organism", "orf"]],
              left_on="qacc", right_on="orf", how="outer", suffixes=("_blastx", "_tblastx"))
    final_matches = final_matches.loc[final_matches["Gene.systematicName"]==final_matches["gene"]] # where the blastx gene is the same as the tblastx gene
    
    final_matches.loc[final_matches.orf.duplicated()] # duplicated orf to gene matches exist
    
    # For each orf, keep the row with the smallest evalue and largest pident
    final_matches["evalue_min"] = final_matches[["evalue_blastx", "evalue_tblastx"]].min(axis=1)
    final_matches["pident_max"] = final_matches[["pident_blastx", "pident_tblastx"]].mean(axis=1)
    final_matches_sub = final_matches.groupby("orf").apply(lambda group:\
        group[(group.evalue_min == group.evalue_min.min()) &\
              (group.pident_max == group.pident_max.max())]).reset_index(drop=True)
    
    # Now drop the duplicate rows
    final_matches_sub = final_matches_sub[["orf", "gene", "organism", "evalue_min", "pident_max"]].\
        drop_duplicates()

    org_dupes = set() # some dupes are because the organism names are slightly different
    for i, r in final_matches_sub.loc[final_matches_sub.orf.duplicated(keep=False)].iterrows():
        org_dupes.add(tuple(final_matches_sub.loc[final_matches_sub.orf==r.orf, "organism"]))
    """>>> org_dupes
    {('Saccharomyces cerevisiae S288C', 'Vanderwaltozyma polyspora DSM 70294 ASM15003v1'),
    ('Saccharomyces cerevisiae S288C', 'Saccharomyces paradoxus CBS432'),
    ('Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae S288C'),
    ('Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae'),
    ('Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae', 'Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae')}"""
    
    for i, r in final_matches_sub.loc[final_matches_sub.orf.duplicated(keep=False)].iterrows():
        if tuple(final_matches_sub.loc[final_matches_sub.orf==r.orf, "organism"]) in \
            [('Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae S288C'),
             ('Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae'),
             ('Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae',
              'Saccharomyces cerevisiae S288C', 'Saccharomyces cerevisiae')]:
            final_matches_sub.loc[final_matches_sub.orf==r.orf, "organism"] = "Saccharomyces cerevisiae S288C"
    
    final_matches_sub = final_matches_sub.drop_duplicates(keep="first")
    
    # Save the final orf to gene map
    # 5439-YLL039C is the only orf that mapped to both S288C and V. polyspora, the evalue and pident are the same, so I'll keep the S288C gene.
    final_matches_sub.loc[final_matches_sub.organism != "Vanderwaltozyma polyspora DSM 70294 ASM15003v1"].\
        to_csv("Data/Peter_2018/final_map_orf_to_gene_with_blast_CORRECTED.tsv", index=False, sep="\t")
    
    final_matches_sub.loc[final_matches_sub.organism != "Vanderwaltozyma polyspora DSM 70294 ASM15003v1"].\
        drop(columns=["evalue_min", "pident_max"]).to_csv(
            "Data/Peter_2018/final_map_orf_to_gene_CORRECTED.tsv", index=False, sep="\t")
    
    final_matches_sub.organism.value_counts()
    """>>> final_matches_sub.organism.value_counts()
    organism
    Saccharomyces cerevisiae S288C                    5555
    Saccharomyces paradoxus CBS432                     326
    Saccharomyces cerevisiae                            44
    Torulaspora delbrueckii CBS 1146 ASM24337v1          5
    Drosophila guanche DGUA_6                            4
    Vanderwaltozyma polyspora DSM 70294 ASM15003v1       1
    Name: count, dtype: int64"""
        
    '''Old code below:
    tblastx_strict["gene"] = tblastx_strict["qacc"]
    tblastx_strict["orf"] = tblastx_strict["sacc"]
    tblastx_strict.gene = tblastx_strict.gene.replace(sgd_matches_map) # replace nt tblastx sacc with sgd gene IDs from blastx results
    tblastx_strict.gene = tblastx_strict.gene.replace(ncbi_matches_map) # replace nt tblastx sacc with ncbi gene IDs from blastx results
    tblastx_strict.gene = tblastx_strict.gene.replace(sgd_19_genes_map) # replace nt tblastx sacc with sgd gene IDs from blastx results
    for i in tqdm(range(tblastx_strict.shape[0])):
        for key in ncbi_tblastx_map.keys():
            if tblastx_strict.gene[i] in key:
                tblastx_strict.gene[i] = ncbi_tblastx_map[key] # replace nt tblastx sacc with ncbi gene IDs from tblastx results
    sum(tblastx_strict.gene==tblastx_strict.qacc) # check how many genes were replaced
    # tblastx_strict.to_csv('Data/BLAST_nr_db/nr_blastx_match_organisms/tblastx_strict.txt', sep='\t', index=False) # save to file without removing duplicates
    tblastx_strict["new_organism"] = tblastx_strict.apply(fill_missing,args=('organism', 'gene', ncbi_tblastx_org_map), axis=1) # fill missing organism info
    final_matches = blastx.merge(tblastx_strict, on=["orf", "gene"], how="inner")
    # final_matches.to_csv("Data/Peter_2018/orf_gene_blastx_tblastx_matches.txt", sep="\t", index=False) # save to file without removing duplicates
    
    ### Determine the final gene matches based on lowest e-value and highest percent identity
    # final_matches = blastx.merge(tblastx_strict, on=["orf"], how="inner")
    # final_matches["Final_gene_eval"] = final_matches.apply(get_final, args=("evalue_y", "evalue_x", "gene_x", "gene_y", True), axis=1)
    # final_matches["Final_gene_pident"] = final_matches.apply(get_final, args=("pident_x", "pident_y", "gene_x", "gene_y"), axis=1)
    # final_matches["new_organism"].fillna("missing", inplace=True)
    # final_matches["organism_y"].fillna("missing", inplace=True)
    # final_matches["Final_organism"] = final_matches.apply(get_final, args=("pident_x", "pident_y", "new_organism", "organism_y"), axis=1)
    # final_matches[['orf', 'gene_x', 'new_organism', 'pident_x', 'gene_y', 'organism_y', 'pident_y', 'Final_gene', 'Final_organism']] # compare blastx and tblastx results

    ### Remove duplicate orf to gene matches
    final_matches_sub = final_matches[["orf", "gene", "evalue_x", "evalue_y", "pident_x", "pident_y", "new_organism_x", "new_organism_y"]].drop_duplicates() # subset columns
    final_matches_sub[final_matches_sub.orf.duplicated(keep=False)] # check duplicate orf to gene matches due to difference in pident
    
    ### For each orf, keep the row with the smallest evalue and largest pident
    final_matches_sub["evalue_av"] = final_matches_sub[["evalue_x", "evalue_y"]].mean(axis=1)
    final_matches_sub["pident_av"] = final_matches_sub[["pident_x", "pident_y"]].mean(axis=1)
    final_matches_sub = final_matches_sub.groupby("orf").apply(lambda group:\
        group[(group.evalue_av == group.evalue_av.min()) &\
              (group.pident_av == group.pident_av.max())]).reset_index(drop=True)

    ### Aggregate all gene matches per orf
    sum(final_matches_sub.new_organism_x==final_matches_sub.new_organism_y) # check if the two columns are the same
    final_matches_sub["organism"] = final_matches_sub["new_organism_x"]
    final_matches_sub.drop(columns=["new_organism_x", "new_organism_y"], inplace=True)
    sum(final_matches_sub.organism.isna()) # check if there are any missing organism information
    final_matches_sub.rename(columns={"evalue_x":"evalue_blastx",
                                    "evalue_y":"evalue_tblastx",
                                    "pident_x":"pident_blastx",
                                    "pident_y":"pident_tblastx"}, inplace=True)
    # some orfs were lost in the filtering process, so add them back
    orfs_lost = final_matches[~final_matches.orf.isin(final_matches_sub.orf)]
    orfs_lost.new_organism_x==orfs_lost.new_organism_y
    orfs_lost = orfs_lost[["orf", "gene", "evalue_x", "evalue_y", "pident_x",
                           "pident_y", "new_organism_x"]].drop_duplicates()
    orfs_lost["evalue_av"] = orfs_lost[["evalue_x", "evalue_y"]].mean(axis=1)
    orfs_lost["pident_av"] = orfs_lost[["pident_x", "pident_y"]].mean(axis=1)
    orfs_lost = orfs_lost.groupby("orf").apply(lambda group:\
        group[(group.evalue_av == group.evalue_av.min())]).reset_index(drop=True)
    orfs_lost.rename(columns={"evalue_x":"evalue_blastx",
                              "evalue_y":"evalue_tblastx",
                              "pident_x":"pident_blastx",
                              "pident_y":"pident_tblastx",
                              "new_organism_x": "organism"}, inplace=True)
    final_matches_sub = pd.concat([final_matches_sub, orfs_lost], ignore_index=True)
    # final_matches_sub.to_csv("Data/Peter_2018/final_map_orf_to_gene_with_blast.txt", sep="\t", index=False) # save to file without removing duplicates
    orf2gene = final_matches_sub.groupby("orf").agg({"gene": lambda x: " // ".join(x), "organism": lambda x: " // ".join(x)}).reset_index() # some orfs map to multiple genes
    gene2orf = final_matches_sub.groupby("gene").agg({"orf": lambda x: " // ".join(x)}).reset_index() # some genes map to multiple orfs
    gene2orf = gene2orf.merge(final_matches_sub[["gene", "organism"]].drop_duplicates(), on="gene", how="left") # add organism information
    # orf2gene.to_csv("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t", index=False) # final orf to gene mapping
    # gene2orf.to_csv("Data/Peter_2018/final_map_gene_to_orf.txt", sep="\t", index=False)
'''