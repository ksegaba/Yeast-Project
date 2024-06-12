#!/bin/bash
# Calculate the percent missing data in the genotype matrix that has
# already been filtered using VCFtools. This script also calculates the
# genotyping rate.

cd Shiu_Lab/Project/Results/missing_data

/mnt/home/seguraab/Shiu_Lab/Project/External_software/plink2 \
    --vcf /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/1011Matrix_diploid750_biallelic_maf0.05_11032022.recode.vcf \
    --missing --allow-extra-chr

# Calculate the genotyping rate
/mnt/home/seguraab/Shiu_Lab/Project/External_software/plink2 \
    --vcf /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/1011Matrix_diploid750_biallelic_maf0.05_11032022.recode.vcf \
    --genotyping-rate --allow-extra-chr