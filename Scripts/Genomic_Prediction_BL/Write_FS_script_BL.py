#!/bin/python
"""
This code generates slurm job submission scripts for feature selection (FS)
and Bayesian LASSO (BL) using the FS files for each trait. 

Inputs:
- path1 = path to BL pipeline scripts e.g. /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BL
- path2 = path to save FS files & RF output  e.g /mnt/scratch/seguraab/yeast_project/yeast_BL_results
- path3 = path to BL input files (genomic) e.g. /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
- save = path to save job_submission files e.g. /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/BL
- trait = trait of interest e.g. YPACETATE
- start = minimum number of features e.g. 500
- stop = maximum number of features e.g. 50000
- step = step size (generate a FS file at each step) e.g. 500
- runFS = sbatch slurm script to create FS files (y/n), default is y
- runRF = sbatch slurm script to run FS on RF (y/n), default is n

Outputs:
    - Script to generate feature selection files containing a subset of features
    - Script to run BL on each of the feature selection files.

Command: 
conda activate
python Write_FS_script_BL.py -path1 /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BL -path2 /mnt/scratch/seguraab/yeast_project/yeast_BL_results -path3 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018 -save /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/BL

Kenia Segura Abá
01/19/2022
"""

import sys,os,argparse

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

#~~~~~~~~~~~~~# I will have to determine this later, maybe I can use the same as the RF
# And, maybe I don't need FS files for each trait, maybe once is enough?? 
# Also, are the feats the same as the rrBLUP?? Does it matter?
"""def make_FS_files(path1, path2, path3, save, trait, start, stop, step, name1):
    # Write feature selection slurm job submission script
    out1 = open(name1, 'w')
    out1.write(f'#!/bin/sh --login \
                \n#SBATCH --time=8:00:00 \
                \n#SBATCH --ntasks=1 \
                \n#SBATCH --cpus-per-task=2 \
                \n#SBATCH --mem=36G \
                \n#SBATCH --job-name BL_make_FS_{trait} \
                \n#SBATCH --output=%x_%j \
                \ncd {path2}/ \
                \nmodule purge \
                \nmodule load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2 \
                \ntrait={trait}\n')

    for i in range(start,stop+step,step):
        out1.write('python %s/Feature_Selection.py -df %s/geno_BL_${trait}.csv -alg randomforest -n %i -test /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt -y_name Y -sep , -type r -save feat_BL_${trait}_top -scores T\n'%(path1,path3,i))
    out1.write('\nscontrol show job $SLURM_JOB_ID')
    out1.close()"""

"""def make_BL_FS_files(path1, path2, path3, save, trait, start, stop, step, name2):
    # Write RF slurm job submission script
    out2 = open(name2, 'w')
    out2.write(f'#!/bin/sh --login \
                \n#SBATCH --time=4:00:00 \
                \n#SBATCH --ntasks=2 \
                \n#SBATCH --cpus-per-task=2 \
                \n#SBATCH --mem=36G \
                \n#SBATCH --job-name BL_FS_job_{trait} \
                \n#SBATCH --output=%x_%j \
                \ncd {path2} \
                \nmodule purge \
                \nmodule load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2 \
                \ntrait={trait}\n')
    
    for i in range(start,stop+step,step):
        out2.write('Rscript %s/BL_Model.R %s/geno.csv %s/pheno.csv all %s/Test.txt ${trait} %s/CVFs.csv 5 10 BL_geno'%(path1,path3,path3,path3,path3,i,i))
    out2.write('\nscontrol show job $SLURM_JOB_ID')
    out2.close()"""

def make_BL_files(path1, path2, path3, trait, name):
    # Write RF slurm job submission script
    out = open(name, 'w')
    out.write(f'#!/bin/sh --login \
                \n#SBATCH --time=15:00:00 \
                \n#SBATCH --ntasks=1 \
                \n#SBATCH --cpus-per-task=2 \
                \n#SBATCH --mem=36G \
                \n#SBATCH --job-name BL_job_{trait} \
                \n#SBATCH --output=%x_%j \
                \ncd {path2} \
                \nmodule purge \
                \nmodule load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2\n')
    
    out.write('Rscript %s/BL_Model.R %s/geno.csv %s/pheno.csv all %s/Test.txt %s %s/CVFs.csv 5 10 BL_geno %s'%(path1,path3,path3,path3,trait,path3,path2))
    out.write('\nscontrol show job $SLURM_JOB_ID')
    out.close()


def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='This code generates slurm job submission scripts for feature selection (FS) and random forest (RF).')
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument('-path1', help='path to RF pipeline scripts', required=True)
    req_group.add_argument('-path2', help='path to save feature selection files & RF output', required=True)
    req_group.add_argument('-path3', help='path to RF input files (genomic)' ,required=True)
    req_group.add_argument('-save', help='path to save job_submission files' ,required=True)
    """req_group.add_argument('-trait', help='trait of interest e.g. YPACETATE', required=True)
    req_group.add_argument('-start', help='minimum number of features' ,required=True)
    req_group.add_argument('-stop', help='maxiumum number of features' ,required=True)
    req_group.add_argument('-step', help='step size (generate a FS file at each step)' ,required=True)
    req_group.add_argument('-runFS', help='submit FS slurm script (y/n)', default='y')
    req_group.add_argument('-runRF', help='submit RF slurm script (y/n)', default='n')"""
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # Arguments 
    path1 = args.path1 # path to RF pipeline scripts e.g. /mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
    path2 = args.path2 # path to save feature selection files & RF output  e.g /mnt/scratch/seguraab/yeast_project/yeast_BL_results
    path3 = args.path3 # path to RF input files (genomic) e.g. /mnt/scratch/seguraab/yeast_project/yeast_BL_results
    save = args.save # path to save job_submission files e.g. /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts
    #trait = args.trait # trait of interest e.g. YPACETATE
    #start = int(args.start) # minimum number of features e.g. 500
    #stop = int(args.stop) # maximum number of features e.g. 50000
    #step = int(args.step) # step size (generate a FS file at each step) e.g. 500
    
    # Write slurm scripts
    """name1 = '%s/FS_files_for_RF_%s.slurm'%(save, trait)
    name2 = '%s/RF_FS_job_%s.slurm'%(save,trait)
    make_FS_files(path1, path2, path3, save, trait, start, stop, step, name1)
    make_RF_FS_files(path1, path2, path3, save, trait, start, stop, step, name2)"""

    traits = ['YPACETATE','YPDCAFEIN40','YPDHU','YPETHANOL','YPD14','YPDCAFEIN50','YPDKCL2M','YPGALACTOSE','YPD40','YPDCHX05','YPDLICL250MM','YPGLYCEROL','YPD42','YPDCHX1','YPDMV','YPRIBOSE','YPD6AU','YPDCUSO410MM','YPDNACL15M','YPSORBITOL','YPDANISO10','YPDDMSO','YPDNACL1M','YPXYLOSE','YPDANISO20','YPDETOH','YPDNYSTATIN','YPDANISO50','YPDFLUCONAZOLE','YPDSDS','YPDBENOMYL200','YPDFORMAMIDE4','YPDSODIUMMETAARSENITE','YPDBENOMYL500','YPDFORMAMIDE5']
    for trait in traits:
        print(trait)
        # Write slurm scripts
        name = '%s/BL_job_submission_%s.slurm'%(save,trait)
        make_BL_files(path1, path2, path3, trait, name)
        # Submit the slurm scripts
        os.system(f'sbatch {name}')

    # Submit the slurm scripts
    """if (args.runFS == 'y'):
        os.system(f'sbatch {name1}') # to create feature selection files
    if (args.runRF == 'y'):
        os.system(f'sbatch {name2}') # to run feature selection on RF"""
    

if __name__ == '__main__':
    main()