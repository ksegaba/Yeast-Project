#!/usr/bin/env python3
"""
This code generates and submits slurm job submission scripts to run rrBLUP feature selection for each trait. 

Inputs:
- path1 = path to rrBLUP pipeline scripts e.g. /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass
- path2 = path to save rrBLUP output  e.g /mnt/scratch/seguraab/yeast_project/SNP_yeast_rrBLUP_results
- path3 = path to rrBLUP input files (genomic) e.g. /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
- save = path to save job_submission files e.g. /mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/Feature_Selection
- start = minimum number of features e.g. 500 (based on Marker_top###.txt files)
- stop = maximum number of features e.g. 50000 (based on Marker_top###.txt files)
- step = step size e.g. 500 (based on Marker_top###.txt files)

Outputs:
    - Script to run rrBLUP feature selection

Command:
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass
path2=/mnt/scratch/seguraab/yeast_project/SNP_yeast_rrBLUP_results
path3=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
save=/mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/Feature_Selection/rrBLUP
python Write_FS_script_rrBLUP.py -path1 $path1 -path2 $path2 -path3 $path3 -save $save -snp True -submit y

Kenia Segura Aba
01/22/2022
"""

import warnings
import sys,os,argparse
import datatable


def warn(*args, **kwargs):
	pass


warnings.warn = warn


def make_rrBLUP_FS_files(path1, path2, path3, trait, start, stop, step, batch, snp=False, orf=False, cno=False):
    """ Write rrBLUP slurm job submission script """
    if snp:
        name = '%s/rrBLUP_FS_%s.sb' % (save, trait)
    if orf: 
        name = '%s/rrBLUP_FS_orf_%s.sb' % (save, trait)
    if cno: 
        name = '%s/rrBLUP_FS_cno_%s.sb' % (save, trait)
    out = open(name, 'w')
    out.write(f'#!/bin/sh --login \
                \n#SBATCH --array={start}-{stop}:{batch} \
                \n#SBATCH --time=12:00:00 \
                \n#SBATCH --ntasks=1 \
                \n#SBATCH --cpus-per-task=1 \
                \n#SBATCH --mem=20G \
                \n#SBATCH --job-name rrBLUP_FS_{trait} \
                \n#SBATCH --output=%x_%j.out \
                \ncd {path2} \
                \nmodule purge \
                \nmodule load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3 \
                \ntrait={trait} \
                \nPIPE={path1} \
                \nDATA={path3} \
                \nstart_ind=$SLURM_ARRAY_TASK_ID \
                \nend_ind=$(( $start_ind + {batch} - {step}))  \
                \necho "this job start from $start_ind to $end_ind" \
                \nfor i in `seq $start_ind {step} $end_ind`; do')
    if snp:
        out.write('\nRscript ${PIPE}/13_rrBLUP_training_test_split_fread_predict_values.r ${DATA}/geno.csv ${DATA}/pheno.csv Markers_top${i}.txt $trait ${DATA}/Test.txt 5 10 ${DATA}/CVFs.csv rrBLUP_geno_Markers_top${i}')
    if orf:
        out.write('\nRscript ${PIPE}/13_rrBLUP_training_test_split_fread_predict_values.r ${DATA}/ORFs_pres_abs.csv ${DATA}/pheno.csv Markers_top${i}.txt $trait ${DATA}/Test.txt 5 10 ${DATA}/CVFs.csv rrBLUP_orf_Markers_top${i}')
    if cno:
        out.write('\nRscript ${PIPE}/13_rrBLUP_training_test_split_fread_predict_values.r ${DATA}/ORFs_no_NA.csv ${DATA}/pheno.csv Markers_top${i}.txt $trait ${DATA}/Test.txt 5 10 ${DATA}/CVFs.csv rrBLUP_cno_Markers_top${i}')
    out.write('\ndone\n \
               \nscontrol show job $SLURM_JOB_ID')
    out.close()
    return name


if __name__ == '__main__':
    # Argument parser
    parser = argparse.ArgumentParser(description='This code generates and submits slurm job submission scripts for rrBLUP feature selection for each trait.')
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument(
        '-path1', help='path to rrBLUP pipeline scripts', required=True)
    req_group.add_argument(
        '-path2', help='path to save feature selection files & rrBLUP output', required=True)
    req_group.add_argument(
        '-path3', help='path to rrBLUP input files (genomic)' ,required=True)
    req_group.add_argument(
        '-save', help='path to save job_submission files' ,required=True)
    # Optional Input
    req_group.add_argument(
        '-start', help='minimum number of features (based on Marker_top###.txt files)', default=500)
    req_group.add_argument(
        '-stop', help='maxiumum number of features (based on Marker_top###.txt files)', default=50000)
    req_group.add_argument(
        '-step', help='step size (based on Marker_top###.txt files)',  default=500)
    req_group.add_argument(
        '-batch', help='number of features each array task gets e.g 2500', default=2500)
    req_group.add_argument(
        '-snp', help='set to True if features are SNPs', default=False)
    req_group.add_argument(
        '-orf', help='set to True if features are ORF presence/absence', default=False)
    req_group.add_argument(
        '-cno', help='set to True if features are copy number', default=False)
    req_group.add_argument(
        '-submit', help='submit slurm scripts (y/n)', default='n')
    # Help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args() # Read arguments

    # Arguments 
    path1 = args.path1
    path2 = args.path2
    path3 = args.path3
    save = args.save
    start = int(args.start)
    stop = int(args.stop)
    step = int(args.step)
    batch = int(args.batch)

    # Write slurm scripts
    traits = ['YPACETATE','YPDCAFEIN40','YPDHU','YPETHANOL','YPD14','YPDCAFEIN50','YPDKCL2M','YPGALACTOSE','YPD40','YPDCHX05','YPDLICL250MM','YPGLYCEROL','YPD42','YPDCHX1','YPDMV','YPRIBOSE','YPD6AU','YPDCUSO410MM','YPDNACL15M',
              'YPSORBITOL','YPDANISO10','YPDDMSO','YPDNACL1M','YPXYLOSE','YPDANISO20','YPDETOH','YPDNYSTATIN','YPDANISO50','YPDFLUCONAZOLE','YPDSDS','YPDBENOMYL200','YPDFORMAMIDE4','YPDSODIUMMETAARSENITE','YPDBENOMYL500','YPDFORMAMIDE5']
    for trait in traits:
        print(trait)
        name = make_rrBLUP_FS_files(path1, path2, path3, trait, start, stop, 
                                    step, batch, snp=args.snp, orf=args.orf, cno=args.cno)
        # Submit the slurm scripts
        if (args.submit == 'y'):
            os.system(f'sbatch {name}')


