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
    del geno

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
            if type not in ["chromosome", "region"]:
                if chr not in G: # dictionary for each chromosome number
                    G[chr] = {}
                # if type == "gene": # search gene coding sequences
                #     name = tem[8].split('ID=')[1].split(';')[1].split('Name=')[1] # gene name
                #     start = tem[3] # gene start position
                #     stop = tem[4] # gene stop position
                #     strand = tem[6] # forward (+) or reverse (-)
                #     G[chr][name] = [start, stop, strand]
                if type not in G[chr]:
                    G[chr][type] = {}
                if tem[8].startswith("ID="):
                    name = tem[8].split('ID=')[1].split(';')[1].split('Name=')[1] # gene name
                elif tem[8].startswith("Parent="):
                    name = tem[8].split("Name=")[1].split("_")[0] # gene name
                start = tem[3] # gene start position
                stop = tem[4] # gene stop position
                strand = tem[6] # forward (+) or reverse (-)
                G[chr][type][name] = [start, stop, strand]
    
    # out = open("all_genes_S288C.txt", "w") # write to file
    # for key, value in G.items(): # key is chromosome
    #     for k, v in value.items(): # k is gene name
    #         out.write("%s,%s,%s,%s,%s\n" % (key, k, v[0], v[1], v[2])) # v is start, stop, strand
    # out.close()
    
    out = pd.DataFrame.from_dict({(chr, type, gene): info for chr, values in G.items() for type, genes in values.items() for gene, info in genes.items()},
        orient="index", columns=["start", "stop", "strand"])
    out.index = pd.MultiIndex.from_tuples(out.index, names=["chr", "type", "gene"])
    out["start"] = out["start"].astype(int)
    out["stop"] = out["stop"].astype(int)
    out.to_csv("all_genes_types_S288C.txt", sep=",")
    out.reset_index(inplace=True)
    
    # Match SNPs to genes
    chr_map = {f"chromosome{i}": f"chr{roman}" for i, roman in enumerate(
        ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", 
         "XI", "XII", "XIII", "XIV", "XV", "XVI"], start=1)}
    map = {}
    for s in tqdm(snps, total=len(snps), desc="Matching SNPs to genes"): # loop through bi-allelic snps
        # map[s] = []
        chr = s.split("_")[0] # chromosome number
        pos = int(s.split("_")[1]) # position of snp
        # search reference genes to find where snp falls
        chr_short = chr_map[chr]
        s_out = out.loc[(out.chr == chr_short) & (out.start <= pos) & (out.stop >= pos),]
        map[s] = s_out[["type", "gene"]].to_dict(orient="records") # map gene(s) to snp
        
        # for CHR, value in G.items(): # loop through reference genes
        #     for gene, v in value.items():
        #         start = int(v[0])
        #         stop = int(v[1])
        #         if (chr == "chromosome1" and CHR == "chrI"): # chromosomes match
        #             if (pos >= start and pos <= stop): # snp falls within genic region
        #                 map[s].append(gene) # map gene to snp
        #         elif (chr == "chromosome2" and CHR == "chrII"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome3" and CHR == "chrIII"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome4" and CHR == "chrIV"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome5" and CHR == "chrV"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome6" and CHR == "chrVI"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome7" and CHR == "chrVII"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome8" and CHR == "chrVIII"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome9" and CHR == "chrIX"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome10" and CHR == "chrX"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome11" and CHR == "chrXI"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome12" and CHR == "chrXII"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome13" and CHR == "chrXIII"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome14" and CHR == "chrXIV"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome15" and CHR == "chrXV"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
        #         elif (chr == "chromosome16" and CHR == "chrXVI"):
        #             if (pos >= start and pos <=stop):
        #                 map[s].append(gene)
                    
    print(len(map))
    rows = []
    for snp, types in map.items():
        if not types:  # empty list
            # rows.append({"snp": snp, "type": None, "gene": None})
            rows.append([snp, None, None])
        else:
            for type in types:
                # rows.append({"snp": snp, "type": type.get("type"), "gene": type.get("gene")})
                rows.append([snp, type.get("type"), type.get("gene")])

    # # out = open("biallelic_snps_diploid_and_S288C_genes.txt", "w")
    # out = open("biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv", "w")
    # for s in tqdm(snps):
    #     chr, pos = s.split("_")[0], s.split("_")[1]
    #     # write the snp with all the genes it maps to
    #     if len(map[s]) == 1:
    #         out.write("%s\t%s\t%s\t%s\n" % (s, chr, pos, map[s][0]))
    #     elif len(map[s]) == 0:
    #         out.write("%s\t%s\t%s\tintergenic\n" % (s, chr, pos))
    #     else:
    #         out.write("%s\t%s\t%s\t%s\n" % (s, chr, pos, ",".join(map[s])))
    # out.close()
    out2 = pd.DataFrame(rows, columns=["snp", "type", "gene"])
    # drop duplicate rows
    out2.gene.str.contains("_").sum() # 172941 have suffixes like _telomere or _mRNA and others
    out2["gene"] = out2.gene.str.split("_").str[0]
    out2.gene.str.contains("_").sum() # 0
    out2.drop_duplicates(keep="first", inplace=True)
    out2.to_csv("biallelic_snps_diploid_and_S288C_genes_CORRECTED_types_expanded.tsv", sep="\t", index=False)
    
    '''Ensure that biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv and
    biallelic_snps_diploid_and_S288C_genes_CORRECTED_types_expanded.tsv have the same
    contents for the 'gene' type.'''
    og = dt.fread('/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv').to_pandas()
    og.columns = ['snp', 'chr', 'pos', 'gene']
    
    # subset only the 'genic' type rows
    out2_gene = out2.loc[out2.type == 'gene', ['snp', 'gene']]
    og_gene = og.loc[og.gene != 'intergenic', ['snp', 'gene']]
    
    # join the og_gene snp, gene1 and snp, gene2 into the same row
    out2_gene = out2_gene.groupby('snp')['gene'].apply(lambda x: ','.join(x)).reset_index()
    assert out2_gene.shape[0] == og_gene.shape[0] # pass!
    assert set(out2_gene.snp) == set(og_gene.snp) # pass!
    assert set(out2_gene.gene) == set(og_gene.gene) # pass!
    # All checks passed! The 'gene' type contents of both files are the same.


def num_genes_per_type():
    """Determine the number of genes that belong to each type (gene, 5' UTR, mRNA, etc.)"""
    snp_map = dt.fread(
        "/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Scripts/Data_Vis/biallelic_snps_diploid_and_S288C_genes_CORRECTED_types_expanded.tsv").to_pandas()
    
    # reshape snp_map so that each column corresponds to a gene type and the values are the gene names
    snp_map.pivot_table(index="snp", columns="type", values="gene", aggfunc=lambda x: ",".join(x)).to_csv(
        "/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Scripts/Data_Vis/biallelic_snps_diploid_and_S288C_genes_CORRECTED_types_pivoted.tsv", sep="\t")

    # How many genes are in each type category?
    snp_map[["type", "gene"]].drop_duplicates().groupby("type").count()
    #                                         gene
    # type                                    
    #                                     1
    # ARS                                  314
    # ARS_consensus_sequence                36
    # CDS                                 6258
    # LTR_retrotransposon                   33
    # W_region                               2
    # X_element                             26
    # X_element_combinatorial_repeat        20
    # Y_prime_element                       11
    # Y_region                               1
    # Z1_region                              3
    # Z2_region                              2
    # blocked_reading_frame                  5
    # centromere                            14
    # centromere_DNA_Element_I               2
    # centromere_DNA_Element_II             14
    # centromere_DNA_Element_III             2
    # external_transcribed_spacer_region     4
    # five_prime_UTR_intron                 24
    # gene                                6285
    # intein_encoding_region                 2
    # intron                               219
    # long_terminal_repeat                 301
    # mRNA                                6403
    # mating_type_region                     1
    # matrix_attachment_site                 3
    # ncRNA                                 19
    # ncRNA_gene                            19
    # non_transcribed_region                 2
    # noncoding_exon                        55
    # pseudogene                            10
    # rRNA                                   4
    # rRNA_gene                              4
    # recombination_enhancer                 1
    # silent_mating_type_cassette_array      2
    # snRNA                                  2
    # snRNA_gene                             2
    # snoRNA                                30
    # snoRNA_gene                           30
    # tRNA                                  23
    # tRNA_gene                             23
    # telomerase_RNA                         1
    # telomerase_RNA_gene                    1
    # telomere                              28
    # telomeric_repeat                       1
    # transposable_element                  30
    # transposable_element_gene             30
    # uORF                                   1
    
    # How many of the "gene" genes have "mRNA" entries?
    len(set(snp_map.loc[snp_map.type == "gene", "gene"]).intersection(
        set(snp_map.loc[snp_map.type == "mRNA", "gene"]))) # 6285
    
    # What genes have multiple mRNA entries?
    sub = snp_map[["type", "gene"]].drop_duplicates()
    sub.loc[sub.type == "mRNA", "gene"].nunique() # 6403 unique genes have mRNA entries
    
    # What mRNA genes are not in the "gene" type category?
    len(set(snp_map.loc[snp_map.type == "mRNA", "gene"]).difference(
        set(snp_map.loc[snp_map.type == "gene", "gene"]))) # 118
    set(snp_map.loc[snp_map.type == "mRNA", "gene"]).difference(
        set(snp_map.loc[snp_map.type == "gene", "gene"]))
    # {'YOR369C', 'YGR204C-A', 'YDL160C-A', 'YLR325C', 'YPL098C', 'YDR525W-A',
    #  'YJR048W', 'YHR132W-A', 'YDR154C', 'YCR028C-A', 'YGL037C', 'YPR080W',
    #  'YKL058W', 'YER074W-A', 'YOR094W', 'YOR015W', 'YGR215W', 'YDR119W-A',
    #  'YKL145W-A', 'YHR002W', 'YHR193C-A', 'YHR073C-B', 'YJR047C', 'YML007C-A',
    #  'YDR210W', 'YHR049C-A', 'YDR253C', 'YER088W-B', 'YKL152C', 'YOL052C-A',
    #  'YHR050W-A', 'YPL271W', 'YLR167W', 'YDR079W', 'YGL041C-B', 'YOR020W-A',
    #  'YLR075W', 'YGR038W', 'YDR182W-A', 'YJL047C-A', 'YDL067C', 'YCR024C-A',
    #  'YPL038W-A', 'YBR262C', 'YPL225W', 'YLR327C', 'YCL056C', 'YEL054C', 'YDR100W',
    #  'YBL078C', 'YBR191W-A', 'YPR160C-A', 'YMR013C-A', 'YPR133W-A', 'YML101C',
    #  'YOR300W', 'YHR022C-A', 'YCL058W-A', 'YDL133C-A', 'YDL114W-A', 'YGL029W',
    #  'YMR175W', 'YBR009C', 'YCR075W-A', 'YDL184C', 'YPL037C', 'YHR073W-A',
    #  'YGR020C', 'YDL130W', 'YHR131W-A', 'YKL061W', 'YDR382W', 'YPR149W', 'YBR126W-B',
    #  'YMR256C', 'YJL166W', 'YER044C', 'YNR046W', 'YBL029C-A', 'YOR045W', 'YJL062W-A',
    #  'YOR224C', 'YDR418W', 'YDR045C', 'YPL222C-A', 'YPR052C', 'YLR043C', 'YGL054C',
    #  'YNR001W-A', 'YMR046W-A', 'YMR122W-A', 'YNL030W', 'YLR217W', 'YPL189C-A',
    #  'YIL047C-A', 'YDL208W', 'YML009C', 'YPR160W-A', 'YPR063C', 'YPR102C', 'YPR197C',
    #  'YGR286C', 'YGL007C-A', 'YJR034W', 'YER057C', 'YOR302W', 'YDR183C-A', 'YDR322C-A',
    #  'YPR182W', 'YCR010C', 'YLR355C', 'YIL102C-A', 'YGL221C', 'YEL017C-A', 'YER053C-A',
    #  'YDR524C-A', 'YBR010W', 'YLR347W-A'}


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
