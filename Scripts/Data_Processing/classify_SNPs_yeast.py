 """
Author: Peipei Wang
Modified by: Kenia Segura Abá
Date: October 29, 2021

Classify SNP location as either in gene, coding sequence (CDS), 3' untranslated region (UTR), 5' UTR, exonic, intronic, or intergenic.
Input: 
	gff version 3 file containing whole genome data for Saccharomyces cerevisiae
	gVCF version 4.1 file containing single nucleotide polymorphism (SNP) data for Saccharomyces cerevisiae
"""
import sys,os 
gff = open('/mnt/home/seguraab/Shiu_Lab/Project/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff','r').readlines()

# View SNP types
import datatable
df = datatable.fread('/mnt/home/seguraab/Shiu_Lab/Project/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_edited_kenia.gff', fill=True)
types = datatable.unique(df[2]) # subset types column and remove duplicates
print("SNP Types: ")
for i in range(types.shape[0]): # print out each type in gff
	print(types[i,0])

# Create an empty dictionary for each sequence type
types_of_interest = ['gene', 'mRNA', 'CDS', 'five_prime_UTR_intron']
G = {}
for inl in gff: # loop through each line individually
	if inl.startswith('#') or inl.startswith('>'): # skip commented lines (file description) and genome sequence
		pass
	elif inl.startswith('C') or inl.startswith('G') or inl.startswith('A') or inl.startswith('T'): 
		pass
	else: 
		tem = inl.split('\t') # tab delimeted elements are split into a list
		chr = tem[0] # the first element is the chromosome number location of the sequence
		type = tem[2] # the third element is the sequence type
		if type in types_of_interest: # Exclude all other types
			if type == 'gene':
				name = tem[8].split('ID=')[1].split(';')[1].split('Name=')[1]
				G[name] = {}
			if type == 'mRNA':
				name = tem[8].split('ID=')[1].split(';')[1].split('Name=')[1] # has _mRNA at the end
				gene = tem[8].split('ID=')[1].split(';')[0].split('_')[0] 
				G[gene][name] = {}
			if type == 'CDS':
				name = tem[8].split('Parent=')[1].split(';')[1].split('Name=')[1] # has _CDS at the end
				gene = tem[8].split('Parent=')[1].split(';')[1].split('Name=')[1].split('_')[0] 
				mRNA = tem[8].split('Parent=')[1].split(';')[0] # has _mRNA at the end
				G[gene] = {}
				G[gene][mRNA] = {}
				if 'CDS' not in G[gene][mRNA]:
					G[gene][mRNA]['CDS'] = []
				G[gene][mRNA]['CDS'].append(name)
			if type == 'five_prime_UTR_intron':
				name = tem[8].split('Parent=')[1].split(';')[1].split('Name=')[1] # has_five_prime_UTR_intron at the end
				gene = tem[8].split('Parent=')[1].split(';')[1].split('Name=')[1].split('_')[0]
				mRNA = tem[8].split('Parent=')[1].split(';')[0] # has _mRNA at the end
				G[gene] = {}
				G[gene][mRNA] = {}
				if 'five_prime_UTR_intron' not in G[gene][mRNA]:
					G[gene][mRNA]['five_prime_UTR_intron'] = []
				G[gene][mRNA]['five_prime_UTR_intron'].append(name)
#print(G, " End of G")
out = open('G.txt','w')
out.write('Key\tValue\n')
with open('G.txt') as inf:
	for key, value in G.items():
		out.write('%s\t%s\n'%(key, value))
out.close()

# Obtain the start and end position for each sequence type
GFF = {}
for inl in gff:
	if inl.startswith('#') or inl.startswith('>'): # skip commented lines (file description) and genome sequence
		pass
	elif inl.startswith('C') or inl.startswith('G') or inl.startswith('A') or inl.startswith('T'): 
		pass
	else: 
		tem = inl.split('\t')
		chr = tem[0]
		type = tem[2]
		start = int(tem[3]) # start position of the sequence (numbering starts at 1)
		end = int(tem[4]) # end position of the sequence for this particular 
		strand = tem[6] # strand direction: defined as + (forward) or - (reverse).
		if type in types_of_interest: 
			if tem[8].startswith('ID='):
				name = tem[8].split('ID=')[1].split(';')[1].split('Name=')[1]
				if type not in GFF: # new dictionary for each type
					GFF[type] = {}
				if chr not in GFF[type]: # add chromosome number
					GFF[type][chr] = {}
				if type == 'gene':
					if start not in GFF[type][chr]: # add start, end, strand direction, and gene name
						GFF[type][chr][start] = [strand,end,name]
					else:
						print('%s\t%s\t%s'%(chr,start,name))
						GFF[type][chr][start].append([strand,end,name])
			elif tem[8].startswith('Parent='):
				name = tem[8].split('Parent=')[1].split(';')[1].split('Name=')[1]
				if type not in GFF:
					GFF[type] = {}
				if chr not in GFF[type]:
					GFF[type][chr] = {}
				if type == 'CDS' or type =='five_prime_UTR_intron':
					if name not in GFF[type][chr]:
						GFF[type][chr][name] = [strand,end,start]
					else:
						print('%s\t%s\t%s'%(chr,start,name))
						GFF[type][chr][name].append([strand,end,start])
#print(GFF, " End of GFF")
out = open('GFF.txt','w')
out.write('Key\tValue\n')
with open('GFF.txt') as inf:
	for key, value in GFF.items():
		out.write('%s\t%s\n'%(key, value))
out.close()
'''
# Classify mutation type and location from an input gVCFv4.1 file
def judge(ref,alt): # reference and alternate alleles
	"""Count the number of SNPs and indels in gVCF."""
	T = {}
	if ',' in alt: # multiple alternate alleles
		for snp in alt.split(','): # each alt allele is either a snp or indel
			if len(ref) == len(snp) and '*' not in snp: # * means a missing alt allele info
				T['SNP'] = 1
			else: # if length of reference differs from alternate, then there is an indel
				T['indel'] = 1
	else: # single alternate allele
		if len(ref) == len(alt) and '*' not in alt:
			T['SNP'] = 1
		else:
			T['indel'] = 1
	return '/'.join(sorted(T.keys()))

def location(chr,pos):
	"""Classify location of mutation as either genic or intergenic."""
	genic = 0
	if chr in GFF['gene']: # search for chromosome #
		for start in GFF['gene'][chr]: # search each start position in chromosome
			if len(GFF['gene'][chr][start]) >= 3: ### the start position is unique to a gene
				print("GFF['gene'][chr][start]:", GFF['gene'][chr][start], len(GFF['gene'][chr][start]))
				if start <= pos and GFF['gene'][chr][start][1] >= pos: # start and end position of gene
					genic = 1 # mutation position (pos) falls within a gene
					gene_name = GFF['gene'][chr][start][2]
					print("gen_name:", gene_name)
					res = [] #?
					for mRNA_name in G[gene_name]: ### if in exonic and intronic regions in different transcripts, then location = 'splicing'
						print("mRNA_name:", mRNA_name)
						exonic = 0
						five_UTR = 0
						three_UTR = 0
						for cds_name in G[gene_name][mRNA_name]['CDS']:
							print("cds_name:", cds_name)
							if GFF['CDS'][chr][cds_name][2] <= pos and GFF['CDS'][chr][cds_name][1] >= pos:
								exonic = 1
						if exonic == 0:
							if exonic == 0:
								try:
									for three_name in G[gene_name][mRNA_name]['three_prime_UTR']:
										if GFF['three_prime_UTR'][chr][three_name][2] <= pos and GFF['three_prime_UTR'][chr][three_name][1] >= pos:
											three_UTR = 1
									if three_UTR == 0:
										for five_name in G[gene_name][mRNA_name]['five_prime_UTR_intron']:
											try:
												if GFF['five_prime_UTR_intron'][chr][five_name][2] <= pos and GFF['five_prime_UTR_intron'][chr][five_name][1] >= pos:
													five_UTR = 1
											except:
												print('No five UTR for %s'%gene_name)
								except:
									print('No three UTR for %s'%gene_name)
						res.append([exonic,five_UTR,three_UTR])
					mul = 1
					sum = 0
					sum3 = 0
					sum5 = 0
					for result in res:
						mul = mul*result[0]
						sum = sum + result[0]
						sum3 = sum3 + result[2]
						sum5 = sum5 + result[1]
					if mul == 1:
						location = 'exonic'
					if mul == 0 and sum >= 1:
						location = 'splicing'
					if sum == 0:
						if sum3 >= 1 and sum5 == 0:
							location = 'three_UTR'
						if sum5 >= 1 and sum3 == 0:
							location = 'five_UTR'
						if sum3 >= 1 and sum5 >= 1:
							print('something is wrong with %s_%s'%(chr,pos))
						if sum3 == 0 and sum5 == 0:
							location = 'intronic'
					break
			if len(GFF['gene'][chr][start]) > 3: ### the start is not unique to a gene
				for j in range(3,len(GFF['gene'][chr][start])):
					if start < pos and GFF['gene'][chr][start][j][1] > pos:
						genic = 1
						gene_name = GFF['gene'][chr][start][j][2]
						res = []
						for mRNA_name in G[gene_name]: ### if in exonic and intronic regions in different transcripts, then location = 'splicing'
							exonic = 0
							five_UTR = 0
							three_UTR = 0
							for cds_name in G[gene_name][mRNA_name]['CDS']:
								if GFF['CDS'][chr][cds_name][2] <= pos and GFF['CDS'][chr][cds_name][1] >= pos:
									exonic = 1
							if exonic == 0:
								try:
									for three_name in G[gene_name][mRNA_name]['three_prime_UTR']:
										if GFF['three_prime_UTR'][chr][three_name][2] <= pos and GFF['three_prime_UTR'][chr][three_name][1] >= pos:
											three_UTR = 1
									if three_UTR == 0:
										for five_name in G[gene_name][mRNA_name]['five_prime_UTR']:
											try:
												if GFF['five_prime_UTR'][chr][five_name][2] <= pos and GFF['five_prime_UTR'][chr][five_name][1] >= pos:
													five_UTR = 1
											except:
												print('No five UTR for %s'%gene_name)
								except:
									print('No three UTR for %s'%gene_name)
							res.append([exonic,five_UTR,three_UTR])
						mul = 1
						sum = 0
						sum3 = 0
						sum5 = 0
						for result in res:
							mul = mul*result[0]
							sum = sum + result[0]
							sum3 = sum3 + result[2]
							sum5 = sum5 + result[1]
						if mul == 1:
							location = 'exonic'
						if mul == 0 and sum >= 1:
							location = 'splicing'
						if sum == 0:
							if sum3 >= 1 and sum5 == 0:
								location = 'three_UTR'
							if sum5 >= 1 and sum3 == 0:
								location = 'five_UTR'
							if sum3 >= 1 and sum5 >= 1:
								print('something is wrong with %s_%s'%(chr,pos))
							if sum3 == 0 and sum5 == 0:
								location = 'intronic'
						break
		if genic == 0:
			location = 'intergenic'
	else:
		location = 'intergenic'
	return(location)	

	
out = open(sys.argv[1] + '_classification.txt','w')
out.write('Chr\tPos\tRef\tAlt\tType\tallelic\tLocation\n')
with open(sys.argv[1]) as inf: # input gVCF file
	for inl in inf:
		if not inl.startswith('#'):
			tem = inl.split('\t')
			chr = tem[0] # chromosome number
			pos = int(tem[1]) # basepair position of SNP
			ref = tem[3] # reference allele
			alt = tem[4] # alternate allele
			type = judge(ref,alt) # type of mutation (SNP or indel)
			allelic = len(alt.split(',')) + 1 # number of alleles at this location (including ref and alt); 1=monoallelic, 2=biallelic, etc.
			loc = location(chr,pos) # chromosomal basepair location of mutation
			out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chr,pos,ref,alt,type,allelic,loc))
out.close()
'''



