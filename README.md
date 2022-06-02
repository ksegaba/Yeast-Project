# Yeast-Project
Prediction of *Saccharomyces cerevisiae* fitness in different environments and cross-environment prediction of fitness using transfer learning.

## Project Resources
__Data__
- Genetic interactions (Costanzo): [Data Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.4291s)  
- Whole genome RNA-seq data for 1,000 isolates: [SRA](https://www.ebi.ac.uk/ena/browser/view/PRJEB13017)
- Genetic marker data (.gvcf): [here](http://1002genomes.u-strasbg.fr/files/)
- Phenotype data: 35 conditions (YPD standard is control media used to normalize fitness values), 4 replicates, fitness = colony size normalized

__Literature__
- Peter et al. 2018: https://doi.org/10.1038/s41586-018-0030-5 (genetic markers, phenotype, whole genome RNA-seq)
- Costanzo et al. 2016: https://doi.org/10.1126/science.aaf1420 (genetic interactions)

## Description of files and directories

|File/Directory                                     |Description                                                                         |
|---------------------------------------------------|---------------------------------------------------------------------------------   |
|__Data__                                           |                                                                                    |
|Costanzo_S1/                                       |Data File S1. Raw genetic interaction datasets: Pair-wise interaction format        |
|Costanzo_S2/                                       |Data File S2. Raw genetic interaction datasets: Matrix format                       |
|Peter_2018/                                        |Yeast diploid isolates' bi-allelic SNP and fitness data for 35 growth environments  |
|S288C_reference_genome_R64-2-1_20150113/           |Reference yeast genome S288C files                                                  |
|All_genes_and_pathways_in_S._cerevisiae_S288c.txt  |Yeast (S288C) genes and which pathways they belong to                               |
|All_pathways_S._cerevisiae_S288c.txt               |Pathways and which yeast (S288C) genes are in them                                  |
|                                                   |                                                                                    |
|__Scripts__                                        |                                                                                    |
|All_Algorithms/                                    |Figures containing performance comparisons for all models                           |
|Data_Processing/                                   |Scripts for preprocessing genotype and phenotype matrices                           |
|Data_Vis/                                          |Scripts for visualizing genotype and phenotype matrices                             |
|Genetic_Interactions/                              |Scripts for preprocessing genetic interaction matrix                                |
|Genomic_Prediction_BayesC/                         |BayesC scripts                                                                      |
|Genomic_Prediction_BL/                             |Bayesian LASSO scripts                                                              |
|Genomic_Prediction_GBLUP/                          |GBLUP, OBLUP, and CBLUP scripts                                                     |
|Genomic_Prediction_GCN/                            |Initial attempt at building a graph convolutional neural network                    |
|Genomic_Prediction_RF/                             |Random forest and results analysis scripts
|Genomic_Prediction_rrBLUP/                         |Scripts I wrote or edited from Peipei's rrBLUP code in Genomic_prediction_in_Switchgrass|
|06_classify_SNPs_switchgrass.py                    |Peipei Wang's original code for classifying Switchgrass SNPs                        |
|06_classify_SNPs_yeast.ipynb                       |Jupyter notebook for development purposes                                           |
|06_classify_SNPs_yeast.py                          |Adapted from Peipei's code to classify Yeast SNPs                                   |
|yeast_sERRBLUP_07232021.R                          |sERRBLUP code                                                                       |
|                                                   |                                                                                    |
|__External_software__                              |See the following section
|__Job_Submission_Scripts__                         |Contains SLURM job submission scripts for each prediction model                     |
|__yeast_rrBLUP_results__                           |Input and output files and figures for rrBLUP modelling                             |
|__yeast_RF_results__                               |Output files and figures for RF modelling                                           |


## Description of external software that I'm using or exploring
|Software                                           |Description                                                                         |
|---------------------------------------------------|---------------------------------------------------------------------------------   |
|fastPHASE                                          |Executable for imputation of missing genotypes from population data                 |
|Genomic_prediction_in_Switchgrass/                 |Peipei Wang's code for rrBLUP                                                       |
|GWAS_NN                                            |Code for "Gene-Gene Interaction Detection with Deep Learning"                       |
|ML-Pipeline/                                       |Shiu Lab Machine Learning Pipeline (RF code)                                        |
|phase.2.1.1.linux                                  |PHASE source code https://stephenslab.uchicago.edu/software.html                    |
|tasseladmin-tassel-5-standalone-8b0f83692ccb       |TASSEL5 for kinship and linkage disequilibrium analysis                             |

Google Docs with information about all scripts and their development:
The google drive path to the file is `Segura Abá_ShiuLab/Projects/Yeast GI Network/`.
- [2021_Lab_Notebook](https://docs.google.com/document/d/16pWLJoNUdrJx2gudEZUvN1ArLeKEfzN8NybQsQqfP6A/edit#) 
- [2022_Lab_Notebook](https://docs.google.com/document/d/1aUd3k6bq0C7dGq2EUltONa7x2MRFOlhzxdJTPKoCPU0/edit?skip_itp2_check=true#)


[Git Tutorial](https://docs.google.com/presentation/d/1w7n3-A0R0Fjd1lxWuwJmA41r_ctRp2U6N-fBBfPwbJs/edit?usp=sharing)
