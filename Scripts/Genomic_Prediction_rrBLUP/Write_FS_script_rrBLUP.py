#!/bin/python
"""
This code generates and submits slurm job submission scripts to run rrBLUP feature selection for each trait. 

Inputs:
- path1 = path to rrBLUP pipeline scripts e.g. /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass
- path2 = path to save rrBLUP output  e.g /mnt/scratch/seguraab/yeast_project/yeast_rrBLUP_results
- path3 = path to rrBLUP input files (genomic) e.g. /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
- save = path to save job_submission files e.g. /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/Feature_Selection
- start = minimum number of features e.g. 500 (based on Marker_top###.txt files)
- stop = maximum number of features e.g. 50000 (based on Marker_top###.txt files)
- step = step size e.g. 500 (based on Marker_top###.txt files)

Outputs:
    - Script to run rrBLUP feature selection

Command:
python Write_FS_script_rrBLUP.py -path1 /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass -path2 /mnt/scratch/seguraab/yeast_project/yeast_rrBLUP_results -path3 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018 -save /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/Feature_Selection -start 500 -stop 50000 -step 500

Kenia Segura Abá
01/22/2022
"""

import sys,os,argparse
import datatable

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn


def make_rrBLUP_FS_files(path1, path2, path3, trait, start, stop, step, name):
    # Write rrBLUP slurm job submission script
    out = open(name, 'w')
    out.write(f'#!/bin/sh --login \
                \n#SBATCH --time=8:00:00 \
                \n#SBATCH --ntasks=1 \
                \n#SBATCH --cpus-per-task=1 \
                \n#SBATCH --mem=36G \
                \n#SBATCH --job-name rrBLUP_FS_job_{trait} \
                \n#SBATCH --output=rrBLUP_FS_job_{trait}_%j \
                \ncd {path2} \
                \nmodule purge \
                \nmodule load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2 \
                \ntrait={trait}\n')
    
    for i in range(start,stop+step,step):
        out.write('Rscript %s/13_rrBLUP_training_test_split_fread_predict_values.r %s/geno.csv %s/pheno.csv %s/Markers_top%i.txt $trait %s/Test.txt 5 10 %s/CVFs.csv rrBLUP_geno_Markers_top%i\n'%(path1,path3,path3,path2,i,path3,path3,i))
    out.write('scontrol show job $SLURM_JOB_ID')
    out.close()

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='This code generates and submits slurm job submission scripts for rrBLUP feature selection for each trait.')
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument('-path1', help='path to rrBLUP pipeline scripts', required=True)
    req_group.add_argument('-path2', help='path to save feature selection files & rrBLUP output', required=True)
    req_group.add_argument('-path3', help='path to rrBLUP input files (genomic)' ,required=True)
    req_group.add_argument('-save', help='path to save job_submission files' ,required=True)
    req_group.add_argument('-start', help='minimum number of features (based on Marker_top###.txt files)' ,required=True)
    req_group.add_argument('-stop', help='maxiumum number of features (based on Marker_top###.txt files)' ,required=True)
    req_group.add_argument('-step', help='step size (based on Marker_top###.txt files)' ,required=True)
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    # Arguments 
    path1 = args.path1 # path to rrBLUP pipeline scripts e.g. /mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
    path2 = args.path2 # path to save feature selection files & rrBLUP output  e.g /mnt/scratch/seguraab/yeast_project/yeast_rrBLUP_results
    path3 = args.path3 # path to rrBLUP input files (genomic) e.g. /mnt/scratch/seguraab/yeast_project/yeast_rrBLUP_results
    save = args.save # path to save job_submission files e.g. /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts
    start = int(args.start) # minimum number of features e.g. 500
    stop = int(args.stop) # maximum number of features e.g. 50000
    step = int(args.step) # step size (generate a FS file at each step) e.g. 500
    
    # Write slurm scripts
    traits = ['YPACETATE','YPDCAFEIN40','YPDHU','YPETHANOL','YPD14','YPDCAFEIN50','YPDKCL2M','YPGALACTOSE','YPD40','YPDCHX05','YPDLICL250MM','YPGLYCEROL','YPD42','YPDCHX1','YPDMV','YPRIBOSE','YPD6AU','YPDCUSO410MM','YPDNACL15M','YPSORBITOL','YPDANISO10','YPDDMSO','YPDNACL1M','YPXYLOSE','YPDANISO20','YPDETOH','YPDNYSTATIN','YPDANISO50','YPDFLUCONAZOLE','YPDSDS','YPDBENOMYL200','YPDFORMAMIDE4','YPDSODIUMMETAARSENITE','YPDBENOMYL500','YPDFORMAMIDE5']
    for trait in traits:
        print(trait)
        name = '%s/rrBLUP_FS_job_%s.slurm'%(save,trait)
        make_rrBLUP_FS_files(path1, path2, path3, trait, start, stop, step, name)
        # Submit the slurm scripts
        os.system(f'sbatch {name}') # to run feature selection on rrBLUP

if __name__ == '__main__':
    main()