#!/bin/python
"""
Match bi-allelic SNPs from diploid yeast isolates to S288C reference genes
Kenia Segura Abá

Input:
    genotype file
    gff file
"""

import sys,os,argparse
import json
import datatable as dt
import pandas as pd
from tqdm import tqdm

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Match bi-allelic SNPs from diploid yeast isolates to S288C reference genes")
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument("-geno", help="genotype csv file", required=True)
    req_group.add_argument("-gff", help="reference genome gff file", required=True)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # Read input files
    geno = dt.fread(args.geno)
    gff = open(args.gff, "r").readlines()

    # List of bi-allelic snps from diploid yeast isolates
    snps = geno.names[1:] # exclude "ID" header

    # Extract all the genes and their information from the S288C reference genome
    G = {} # directory to hold gene names, start, stop, and strand
    for inl in gff: # loop through each line in gff file
        if inl.startswith("#"): # skip these lines
            pass
        elif inl.startswith(">"): # exit loop when genome sequence reached
            break
        else:
            tem = inl.split('\t') # tab delimeted elements are split into a list
            chr = tem[0] # the first element is the chromosome number location of the sequence
            type = tem[2] # the third element is the sequence type
            if chr not in G: # dictionary for each chromosome number
                G[chr] = {}
            if type == "gene": # search genes
                name = tem[8].split('ID=')[1].split(';')[1].split('Name=')[1] # gene name
                start = tem[3] # gene start position
                stop = tem[4] # gene stop position
                strand = tem[6] # forward (+) or reverse (-)
                G[chr][name] = [start, stop, strand]

    out = open("all_genes_S288C.txt", "w") # write to file
    for key, value in G.items(): # key is chromosome
        for k, v in value.items(): # k is gene name
            out.write("%s,%s,%s,%s,%s\n" % (key, k, v[0], v[1], v[2])) # v is start, stop, strand
    out.close()

    # Match SNPs to genes
    map = {}
    for s in tqdm(snps, total=len(snps), desc="Matching SNPs to genes"): # loop through bi-allelic snps
        map[s] = []
        chr = s.split("_")[0] # chromosome number
        pos = int(s.split("_")[1]) # position of snp
        for CHR, value in G.items(): # loop through reference genes
            for gene, v in value.items():
                start = int(v[0])
                stop = int(v[1])
                if (chr == "chromosome1" and CHR == "chrI"): # chromosomes match
                    if (pos >= start and pos <= stop): # snp falls within genic region
                        map[s].append(gene) # map gene to snp
                elif (chr == "chromosome2" and CHR == "chrII"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome3" and CHR == "chrIII"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome4" and CHR == "chrIV"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome5" and CHR == "chrV"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome6" and CHR == "chrVI"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome7" and CHR == "chrVII"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome8" and CHR == "chrVIII"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome9" and CHR == "chrIX"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome10" and CHR == "chrX"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome11" and CHR == "chrXI"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome12" and CHR == "chrXII"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome13" and CHR == "chrXIII"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome14" and CHR == "chrXIV"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome15" and CHR == "chrXV"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                elif (chr == "chromosome16" and CHR == "chrXVI"):
                    if (pos >= start and pos <=stop):
                        map[s].append(gene)
                    
    print(len(map))

    # out = open("biallelic_snps_diploid_and_S288C_genes.txt", "w")
    out = open("biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv", "w")
    for s in tqdm(snps):
        chr, pos = s.split("_")[0], s.split("_")[1]
        # write the snp with all the genes it maps to
        if len(map[s]) == 1:
            out.write("%s\t%s\t%s\t%s\n" % (s, chr, pos, map[s][0]))
        elif len(map[s]) == 0:
            out.write("%s\t%s\t%s\tintergenic\n" % (s, chr, pos))
        else:
            out.write("%s\t%s\t%s\t%s\n" % (s, chr, pos, ",".join(map[s])))
    out.close()


def get_snp_locations():
    """Determine the SNPs to the left, within, and to the right of each gene"""
    
    # S288C reference genome GFF file (to get gene start/end bp positions)
    s288c_gff = dt.fread("Data/S288C_reference_genome_R64-3-1_20210421/saccharomyces_cerevisiae_R64-3-1_20210421.gff",
                        skip_to_line=22, fill=True, max_nrows=28386).to_pandas()
    s288c_gff = s288c_gff[s288c_gff["C2"] == "gene"] # gene coding sequences only
    s288c_gff.columns = ["chr", "source", "type", "start", "end", "score",
                         "strand", "frame", "attribute"]
    s288c_gff.insert(9, "gene", s288c_gff.attribute.\
        apply(lambda x: x.split("ID=")[1].split(";")[0])) # get gene names
    
    print((s288c_gff.end - s288c_gff.start).describe()) # gene length stats
    
    # SNP to gene map
    map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
        sep="\t", names=["chr", "pos", "gene"])
    chr_map = {f"chromosome{i}": f"chr{roman}" for i, roman in enumerate(
        ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", 
         "XI", "XII", "XIII", "XIV", "XV", "XVI"], start=1)}
    map_snps.insert(3, "chr_ref", map_snps.chr.map(chr_map))
    unique_genes = map_snps.gene.str.split(",").explode().unique()
    
    # Initialize SNP location dictionary
    snp_locs = {g:{"upstream":[], "genic":[], "downstream":[]} \
                for g in unique_genes if g != "intergenic"}
    
    # Assign SNPs to categories
    for ref_row in tqdm(s288c_gff.itertuples(), total=len(s288c_gff), desc="Processing genes"): # loop through reference genes
        if ref_row.gene in snp_locs.keys(): # gene is in SNP location dictionary
            
            for row in map_snps.itertuples(): # loop through bi-allelic snps
                if (row.chr_ref == ref_row.chr): # chromosomes match
                    if (row.pos >= ref_row.start and row.pos <= ref_row.end): # snp falls within genic region
                        snp_locs[ref_row.gene]["genic"].append(row.Index)
                    
                    if (row.pos < ref_row.start): # snp falls upstream of gene
                        snp_locs[ref_row.gene]["upstream"].append(row.Index)
                    
                    if (row.pos > ref_row.end): # snp falls downstream of gene
                        snp_locs[ref_row.gene]["downstream"].append(row.Index)
    
    # # Sanity check
    # to_check = []
    # for g in tqdm(snp_locs.keys()):
    #     try:
    #         mask = map_snps.gene.explode().eq(g).groupby(level=0).any()
    #         assert len(snp_locs[g]["genic"]) == len(map_snps[mask])
    #         # assert len(snp_locs[g]["genic"]) == \
    #         #     len(map_snps[map_snps.gene.str.split(",").apply(lambda genes: g in genes)])
    #     except AssertionError:
    #         to_check.append(g)
    #
    # len(to_check)
    # with open("geno_snp_locations_to_check.txt", "w") as f:
    #     f.write("\n".join(to_check))
    # I manually checked randomly and found no errors
    
    # Save SNP locations to file
    # with open("geno_snp_locations.json", "w") as f:
    #     json.dump(snp_locs, f)
    
    # snp_locs = json.load(open("snp_locations.json"))
    
    return snp_locs


if __name__ == '__main__':
    main()
    # get_snp_locations() # this was added Feb. 26, 2025
