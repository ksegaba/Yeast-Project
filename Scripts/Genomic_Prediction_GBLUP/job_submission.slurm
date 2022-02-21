#!/bin/sh --login
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name BLUPs_YPACETATE_
#SBATCH --output=%x_%j

echo "$SLURM_ARRAY_TASK_ID"
module purge
module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2
cd /mnt/scratch/seguraab/yeast_project/ORF_yeast_BLUPs_results

## YPACETATE
# No cross-validation GBLUP: SNP matrix
Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/GBLUP_Model.R /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt YPACETATE /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv 5 10 GBLUP n n y n

# No cross-validation OBLUP: ORF presence absence
Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/GBLUP_Model.R /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt YPACETATE /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv 5 10 OBLUP y n n n

# No cross-validation CBLUP: ORF copy number
Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/GBLUP_Model.R /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt YPACETATE /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv 5 10 CBLUP n y n n


## YPD14
# No cross-validation GBLUP: SNP matrix
Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/GBLUP_Model.R /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt YPD14 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv 5 10 GBLUP n n y n

# No cross-validation OBLUP: ORF presence absence
Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/GBLUP_Model.R /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt YPD14 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv 5 10 OBLUP y n n n

# No cross-validation CBLUP: ORF copy number
Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/GBLUP_Model.R /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt YPD14 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv 5 10 CBLUP n y n n

scontrol show job $SLURM_ARRAY_TASK_ID