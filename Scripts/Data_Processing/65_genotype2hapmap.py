"""
Convert genotype matrix ({0,1,-1} encoding) to hapmap format
Author: Peipei Wang
Modified by: Kenia Segura Abá
"""
import sys
inp = open(sys.argv[1],"r").readlines()
out = open(sys.argv[2],"w")


def list2string(s): # 01/27/2022 Kenia: Modified encodings from {0,1} to {0,1,-1}
	newstr = "" 
	for ele in s:
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
		if ele == '0':
			ele = 'AT'
		elif ele == '1':
			ele = 'AA'
		elif ele == '-1':
=======
=======
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
=======
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
		if ele == '0': # heterozygous
			ele = 'TA'
		elif ele == '1': # homozygous minor allele
			ele = 'AA'
		elif ele == '-1': # homozygous major allele
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> origin/main
=======
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
=======
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
			ele = 'TT'
		newstr += ele + "\t"
	return newstr

i = 0
for line in inp:
	li = line.strip().split(',')
	if i == 0:
		rs = "rs#"
		alleles = "alleles"
		chrom = "chrom"
		pos = "pos"
		strand = "strand"
		assembly = "assembly#"
		center = 'center'
		protLSID = 'protLSID'
		assayLSID = "assayLSID"
		panel = "panel"
		QCcode = "QCcode"
		genotype = list2string(li[1:])	
		result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(rs,alleles,chrom,pos,strand,assembly,center,protLSID,assayLSID,panel,QCcode,genotype)
		out.write(result)
		i = i+1
	else:
		rs = "Scerevisiae_"+str(i) # 01/27/2022 Kenia: Changed Ath to Scerevisiae
		alleles = "A/T"
		chrom = li[0].split('_')[0]
		pos = li[0].split('_')[1]
		strand = "+"
		assembly = "NA"
		center = 'NA'
		protLSID = 'NA'
		assayLSID = "NA"
		panel = "NA"
		QCcode = "NA"
		genotype = list2string(li[1:])
		result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(rs,alleles,chrom,pos,strand,assembly,center,protLSID,assayLSID,panel,QCcode,genotype)
		out.write(result)
		i = i+1

out.close()