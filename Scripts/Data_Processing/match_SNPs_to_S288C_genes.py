#!/bin/python
"""
Match bi-allelic SNPs from diploid yeast isolates to S288C reference genes
Kenia Segura Abá

Input:
    genotype file
    gff file
"""

import sys,os,argparse
import datatable

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
    geno = datatable.fread(args.geno)
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
    for s in snps: # loop through bi-allelic snps
        chr = s.split("_")[0] # chromosome number
        pos = s.split("_")[1] # position of snp
        for CHR, value in G.items(): # loop through reference genes
            for gene, v in value.items():
                start = v[0]
                stop = v[1]
                if (chr == "chromosome1" and CHR == "chrI"): # chromosomes match
                    if (pos >= start and pos <= stop): # snp falls within genic region
                        map[s] = [gene] # map gene to snp
                elif (chr == "chromosome2" and CHR == "chrII"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome3" and CHR == "chrIII"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome4" and CHR == "chrIV"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome5" and CHR == "chrV"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome6" and CHR == "chrVI"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome7" and CHR == "chrVII"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome8" and CHR == "chrVIII"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome9" and CHR == "chrIX"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome10" and CHR == "chrX"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome11" and CHR == "chrXI"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome12" and CHR == "chrXII"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome13" and CHR == "chrXIII"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome14" and CHR == "chrXIV"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome15" and CHR == "chrXV"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
                elif (chr == "chromosome16" and CHR == "chrXVI"):
                    if (pos >= start and pos <=stop):
                        map[s] = [gene]
    print(len(map))

    out = open("biallelic_snps_diploid_and_S288C_genes.txt", "w")
    for s in snps:
        if s not in map.keys(): # if snp s is not in map keys, then it's intergenic
            chr, pos = s.split("_")[0], s.split("_")[1]
            out.write("%s,%s,%s,intergenic\n" % (s, chr, pos))
        else:
            for snp, gene in map.items(): # snp s is genic
                if s == snp:
                    chr, pos = s.split("_")[0], s.split("_")[1]
                    out.write("%s,%s,%s,%s\n" % (s, chr, pos, gene[0])) # snp ID, chromosome, position, gene
    out.close()

if __name__ == '__main__':
    main()