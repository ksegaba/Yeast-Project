#!/bin/python
"""
This code generates slurm job submission scripts for feature selection (FS)
and random forest (RF) using the FS files for each trait. 

Inputs:
- path1 = path to RF pipeline scripts e.g. /mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
- path2 = path to save FS files & RF output  e.g /mnt/scratch/seguraab/yeast_project/yeast_rf_results
- path3 = path to RF input files (genomic) e.g. /mnt/scratch/seguraab/yeast_project/yeast_rf_results
- save = path to save job_submission files e.g. /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts
- trait = trait of interest e.g. YPACETATE
- start = minimum number of features e.g. 500
- stop = maximum number of features e.g. 50000
- step = step size (generate a FS file at each step) e.g. 500
- runFS = sbatch slurm script to create FS files (y/n), default is y
- runRF = sbatch slurm script to run FS on RF (y/n), default is n

Outputs:
    - Script to generate feature selection files containing a subset of features
    - Script to run RF on each of the feature selection files.

Kenia Segura Abá
12/01/2021
"""

import sys,os,argparse
import datatable

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def make_FS_files(path1, path2, path3, save, trait, start, stop, step, name1):
    # Write feature selection slurm job submission script
    out1 = open(name1, 'w')
    out1.write(f'#!/bin/sh --login \
                \n#SBATCH --time=100:00:00 \
                \n#SBATCH --ntasks=3 \
                \n#SBATCH --cpus-per-task=1 \
                \n#SBATCH --mem=20G \
                \n#SBATCH --job-name RF_make_FS_{trait} \
                \n#SBATCH --output=%x_%j \
                \ncd {path2}/ \
                \nmodule purge \
                \nmodule load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4 \
                \ntrait={trait}\n')

    for i in range(start,stop+step,step):
        out1.write('python %s/Feature_Selection.py -df %s/geno_rf_${trait}.csv -alg randomforest -n %i -test /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt -y_name Y -sep , -type r -save feat_rf_${trait}_top -scores T\n'%(path1,path3,i))
    out1.write('\nscontrol show job $SLURM_JOB_ID')
    out1.close()

def make_RF_FS_files(path1, path2, path3, save, trait, start, stop, step, name2):
    # Write RF slurm job submission script
    out2 = open(name2, 'w')
    out2.write(f'#!/bin/sh --login \
                \n#SBATCH --time=100:00:00 \
                \n#SBATCH --ntasks=3 \
                \n#SBATCH --cpus-per-task=1 \
                \n#SBATCH --mem=20G \
                \n#SBATCH --job-name RF_FS_job_{trait} \
                \n#SBATCH --output=RF_FS_job_{trait}_%j \
                \ncd {path2} \
                \nmodule purge \
                \nmodule load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4 \
                \ntrait={trait}\n')
    
    for i in range(start,stop+step,step):
        out2.write('python %s/ML_regression.py -df %s/geno_rf_${trait}.csv -sep , -feat %s/feat_rf_${trait}_top_%i -test /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt -alg RF -n_jobs 12 -n 1 -cv_num 5 -save ${trait}_rf_%i -plots t\n'%(path1,path3,path2,i,i))
    out2.write('scontrol show job $SLURM_JOB_ID')
    out2.close()

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='This code generates slurm job submission scripts for feature selection (FS) and random forest (RF).')
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument('-path1', help='path to RF pipeline scripts', required=True)
    req_group.add_argument('-path2', help='path to save feature selection files & RF output', required=True)
    req_group.add_argument('-path3', help='path to RF input files (genomic)' ,required=True)
    req_group.add_argument('-save', help='path to save job_submission files' ,required=True)
    req_group.add_argument('-trait', help='trait of interest e.g. YPACETATE', required=True)
    req_group.add_argument('-start', help='minimum number of features' ,required=True)
    req_group.add_argument('-stop', help='maxiumum number of features' ,required=True)
    req_group.add_argument('-step', help='step size (generate a FS file at each step)' ,required=True)
    req_group.add_argument('-runFS', help='submit FS slurm script (y/n)', default='y')
    req_group.add_argument('-runRF', help='submit RF slurm script (y/n)', default='n')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # Arguments 
    path1 = args.path1 # path to RF pipeline scripts e.g. /mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
    path2 = args.path2 # path to save feature selection files & RF output  e.g /mnt/scratch/seguraab/yeast_project/yeast_rf_results
    path3 = args.path3 # path to RF input files (genomic) e.g. /mnt/scratch/seguraab/yeast_project/yeast_rf_results
    save = args.save # path to save job_submission files e.g. /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts
    trait = args.trait # trait of interest e.g. YPACETATE
    start = int(args.start) # minimum number of features e.g. 500
    stop = int(args.stop) # maximum number of features e.g. 50000
    step = int(args.step) # step size (generate a FS file at each step) e.g. 500
    
    # Write slurm scripts
    name1 = '%s/FS_files_for_RF_%s.slurm'%(save, trait)
    name2 = '%s/RF_FS_job_%s.slurm'%(save,trait)
    make_FS_files(path1, path2, path3, save, trait, start, stop, step, name1)
    make_RF_FS_files(path1, path2, path3, save, trait, start, stop, step, name2)
    
    # Submit the slurm scripts
    if (args.runFS == 'y'):
        os.system(f'sbatch {name1}') # to create feature selection files
    if (args.runRF == 'y'):
        os.system(f'sbatch {name2}') # to run feature selection on RF

if __name__ == '__main__':
    main()