"""
Convert genotypes from one format to 1,0,-1 format
"""
import sys,os
import pandas as pd

df = pd.read_csv('1011Matrix_genotype_matrix',header=0,index_col=2,sep='\t')
df = df.iloc[:,8:]
<<<<<<< HEAD
<<<<<<< HEAD
df = df.replace('0/0', -1)
df = df.replace('0/1', 0)
df = df.replace('1/1', 1)
=======
df = df.replace('0/0', -1) # homozygous reference allele
df = df.replace('0/1', 0) # heterozygous
df = df.replace('1/1', 1) # homozygous alternate allele
>>>>>>> origin/main
=======
df = df.replace('0/0', -1) # homozygous reference allele
df = df.replace('0/1', 0) # heterozygous
df = df.replace('1/1', 1) # homozygous alternate allele
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
df = df.transpose()
df.to_csv('1011GWAS_matrix.txt',header=True, index=True,sep=',')
