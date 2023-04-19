#!/usr/bin/env python3
"""
This code generates slurm job submission scripts for feature selection (FS)
and random forest (RF) using the FS files for each trait. 

Required Arguments:
    path1 (str): path to RF pipeline scripts
    path2 (str): path to save FS files & RF output
    path3 (str): path to RF input files
    save (str): path to save job_submission files

Optional Arguments:
    suffix (str): suffix to add to file names, default is ''
    start (int): minimum number of features, default is 1000
    stop (int): maximum number of features, default is 40000
    step (int): step size (generate a FS file at each step), default is
    batch (int): number of features each array task gets, default is 5000
    base (int): base (generate a FS file from base**start to base**stop), default is 2
    runFS (str): make FS files slurm scripts (y/n), default is n
    runExpFS (str): make FS slurm scripts based on exponentially picking features, default is n
    runRF (str): make RF slurm scripts (y/n), default is n
    snp (bool): set to True if features are SNPs, default is False
    orf (bool): set to True if features are ORF presence/absence, default is False
    cno (bool): set to True if features are ORF copy number, default is False
    submit (str): submit slurm scripts (y/n), default is n

Returns:
    [1] SLURM scripts to generate feature selection files containing a subset of features
    [2] SLURM scripts to run RF on each of the feature selection files.

Examples:
# For SNP-based RF
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
path2=/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results
path3=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
save=/mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/Feature_Selection/RF
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -suffix _250to5000 -start 250 -stop 5000 -step 250 -batch 250 -runFS y -snp True -submit y
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -start 1000 -stop 40000 -step 1000 -batch 5000 -runFS y -snp True -submit y
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -start 1000 -stop 40000 -step 1000 -batch 5000 -runRF y -snp True -submit y
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -start 1 -stop 11 -base 2 -runExpFS y -snp True -submit y

# For ORF/CNO-based RF
path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
path2=/mnt/scratch/seguraab/yeast_project/ORF_SNP_yeast_RF_results
path3=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
save=/mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/FS/RF
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -start 250 -stop 7708 -step 250 -batch 1000 -runRF y -orf True -submit y
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -start 1 -stop 10 -base 2 -runExpFS y -orf True -submit y
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -start 250 -stop 7708 -step 250 -batch 1000 -runRF y -cno True -submit y
python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/Write_FS_script_RF.py -path1 ${path1} -path2 ${path2} -path3 ${path3} -save ${save} -start 1 -stop 10 -base 2 -runExpFS y -cno True -submit y

# Date
12/01/2021
"""

__author__ = "Kenia Segura Abá"

import warnings
import sys
import os
import argparse
from tkinter import N


def warn(*args, **kwargs):
    pass


warnings.warn = warn


def make_FS_files(path1, path2, path3, save, suffix, trait, start, stop, step, batch, snp=False, orf=False, cno=False):
    """ Write feature selection slurm job submission script """
    if snp:
        name1 = '%s/FS_files_for_RF_%s%s.slurm' % (save, trait, suffix)
    if orf:
        name1 = '%s/FS_files_for_RF_orf_%s%s.slurm' % (save, trait, suffix)
    if cno:
        name1 = '%s/FS_files_for_RF_cno_%s%s.slurm' % (save, trait, suffix)
    out1 = open(name1, 'w')
    out1.write(f'#!/bin/sh --login \
                \n#SBATCH --array={start}-{stop}:{batch} \
                \n#SBATCH --time=100:00:00 \
                \n#SBATCH --ntasks=3 \
                \n#SBATCH --cpus-per-task=1 \
                \n#SBATCH --mem=20G \
                \n#SBATCH --job-name RF_make_FS_{trait} \
                \n#SBATCH --output=%x_%j.out \
                \ncd {path2}/ \
                \nmodule purge \
                \nmodule load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4 \
                \ntrait={trait} \
                \npath1={path1} \
                \npath3={path3} \
                \nstart_ind=$SLURM_ARRAY_TASK_ID \
                \nend_ind=$(( $start_ind + {batch} - {step}))  \
                \necho "this job start from $start_ind to $end_ind" \
                \nfor i in `seq $start_ind {step} $end_ind`; do')

    if snp:
        out1.write('\npython ${path1}/Feature_Selection.py -df ${path3}/geno.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $i -test ${path3}/Test.txt -sep , -type r -save feat_rf_${trait}_top -scores T\n')
    if orf:
        out1.write('\npython ${path1}/Feature_Selection.py -df ${path3}/ORFs_pres_abs.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $i -test ${path3}/Test.txt -sep , -type r -save feat_rf_${trait}_top -scores T\n')
    if cno:
        out1.write('\npython ${path1}/Feature_Selection.py -df ${path3}/ORFs_no_NA.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $i -test ${path3}/Test.txt -sep , -type r -save feat_rf_${trait}_top -scores T\n')
    out1.write('done\nscontrol show job $SLURM_JOB_ID')
    out1.close()
    return name1


def make_exp_FS_files(path1, path2, path3, save, suffix, trait, start, stop, base, snp=False, orf=False, cno=False):
    """ Write feature selection slurm job submission script """
    if snp:
        name1 = '%s/FS_exp_files_for_RF_%s%s.slurm' % (save, trait, suffix)
    if orf:
        name1 = '%s/FS_exp_files_for_RF_orf_%s%s.slurm' % (save, trait, suffix)
    if cno:
        name1 = '%s/FS_exp_files_for_RF_cno_%s%s.slurm' % (save, trait, suffix)

    feats = [base**i for i in range(start, stop)]  # No. of features
    out1 = open(name1, 'w')
    out1.write(f'#!/bin/sh --login \
                \n#SBATCH --array={",".join(map(str,feats))} \
                \n#SBATCH --time=3:00:00 \
                \n#SBATCH --ntasks=3 \
                \n#SBATCH --cpus-per-task=1 \
                \n#SBATCH --mem=20G \
                \n#SBATCH --job-name RF_make_FS_{trait} \
                \n#SBATCH --output=%x_%j.out \
                \ncd {path2}/ \
                \nmodule purge \
                \nmodule load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4 \
                \ntrait={trait} \
                \npath1={path1} \
                \npath2={path2} \
                \npath3={path3}\n')
    if snp:
        out1.write('python ${path1}/Feature_Selection.py -df ${path3}/geno.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -sep , -type r -save feat_exp_rf_${trait}_top -scores T\n')
        out1.write('python ${path1}/ML_regression.py -df ${path3}/geno.csv -df2 ${path3}/pheno.csv -y_name ${trait} -sep , -feat ${path2}/feat_exp_rf_${trait}_top_$SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -alg RF -n_jobs 12 -n 10 -cv_num 5 -save ${trait}_exp_rf_$SLURM_ARRAY_TASK_ID -plots t\n\n')
    if orf:
        out1.write('python ${path1}/Feature_Selection.py -df ${path3}/ORFs_pres_abs.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -sep , -type r -save feat_exp_rf_orf_${trait}_top -scores T\n')
        out1.write('python ${path1}/ML_regression.py -df ${path3}/ORFs_pres_abs.csv -df2 ${path3}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_$SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_orf_$SLURM_ARRAY_TASK_ID -plots t\n\n')
    if cno:
        out1.write('python ${path1}/Feature_Selection.py -df ${path3}/ORFs_no_NA.csv -df2 ${path3}/pheno.csv -y_name ${trait} -alg randomforest -n $SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -sep , -type r -save feat_exp_rf_cno_${trait}_top -scores T\n')
        out1.write('python ${path1}/ML_regression.py -df ${path3}/ORFs_no_NA.csv -df2 ${path3}/pheno.csv -y_name ${trait} -sep , -feat ${path3}/feat_exp_rf_${trait}_top_$SLURM_ARRAY_TASK_ID -test ${path3}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_exp_rf_cno_$SLURM_ARRAY_TASK_ID -plots t\n\n')
    out1.write('\nscontrol show job $SLURM_JOB_ID')
    out1.close()
    return name1


def make_RF_files(path1, path2, path3, save, suffix, trait, start, stop, step, batch, snp=False, orf=False, cno=False):
    """ Write RF slurm job submission script """
    if snp:
        name2 = '%s/RF_FS_job_%s%s.sb' % (save, trait, suffix)
    if orf:
        name2 = '%s/RF_FS_job_orf_%s%s.sb' % (save, trait, suffix)
    if cno:
        name2 = '%s/RF_FS_job_cno_%s%s.sb' % (save, trait, suffix)
    out2 = open(name2, 'w')
    out2.write(f'#!/bin/sh --login \
                \n#SBATCH --array={start}-{stop}:{batch} \
                \n#SBATCH --time=100:00:00 \
                \n#SBATCH --ntasks=3 \
                \n#SBATCH --cpus-per-task=2 \
                \n#SBATCH --mem=20G \
                \n#SBATCH --job-name RF_FS_job_{trait} \
                \n#SBATCH --output=%x_%A_%a.out \
                \ncd {path2} \
                \nmodule purge \
                \nmodule load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4 \
                \ntrait={trait} \
                \npath1={path1} \
                \npath2={path2} \
                \npath3={path3} \
                \nstart_ind=$SLURM_ARRAY_TASK_ID \
                \nend_ind=$(( $start_ind + {batch} - {step}))  \
                \necho "this job start from $start_ind to $end_ind" \
                \nfor i in `seq $start_ind {step} $end_ind`; do')
    if snp:
        out2.write('\n\tpython ${path1}/ML_regression.py -df ${path3}/geno.csv -df2 ${path3}/pheno.csv -y_name ${trait} -sep , -feat ${path2}/feat_rf_${trait}_top_$i -test ${path3}/Test.txt -alg RF -n_jobs 12 -n 10 -cv_num 5 -save ${trait}_rf_$i -plots t')
    if orf:
        out2.write('\n\tpython ${path1}/ML_regression.py -df ${path3}/ORFs_pres_abs.csv -df2 ${path3}/pheno.csv -y_name ${trait} -sep , -feat ${path2}/feat_rf_${trait}_top_$i -test ${path3}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_orf_$i -plots t')
    if cno:
        out2.write('\n\tpython ${path1}/ML_regression.py -df ${path3}/ORFs_no_NA.csv -df2 ${path3}/pheno.csv -y_name ${trait} -sep , -feat ${path2}/feat_rf_${trait}_top_$i -test ${path3}/Test.txt -alg RF -n_jobs 12 -n 100 -cv_num 5 -save ${trait}_cno_$i -plots t')
    out2.write('\ndone\n \
                \nscontrol show job $SLURM_JOB_ID')
    out2.close()
    return name2


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
        description='This code generates slurm job submission scripts for feature selection (FS) and random forest (RF).')
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument(
        '-path1', type=str, help='path to RF pipeline scripts', required=True)
    req_group.add_argument(
        '-path2', type=str, help='path to save feature selection files & RF output', required=True)
    req_group.add_argument(
        '-path3', type=str, help='path to RF input files (genomic)', required=True)
    req_group.add_argument(
        '-save', help='path to save job_submission files', required=True)
    # Optional input
    opt_group = parser.add_argument_group(title='OPTIONAL INPUT')
    opt_group.add_argument(
        '-suffix', type=str, help='suffix to add to file names', default='')
    opt_group.add_argument(
        '-start', type=int, help='minimum number of features', default=1000)
    opt_group.add_argument(
        '-stop', type=int, help='maxiumum number of features', default=40000)
    opt_group.add_argument(
        '-step', type=int, help='step size (generate a FS file at each step)', default=1000)
    opt_group.add_argument(
        '-batch', type=int, help='number of features each array task gets e.g. 5000', default=5000)
    opt_group.add_argument(
        '-base', type=int, help='base (generate a FS file from base ** start to base ** stop) e.g. 2', default=2)
    opt_group.add_argument(
        '-runFS', type=str, help='make FS files slurm scripts (y/n)', default='n')
    opt_group.add_argument(
        '-runExpFS', type=str, help='make FS files slurm scripts based on exponentially picking features, default is n', default='n')
    opt_group.add_argument(
        '-runRF', type=str, help='make RF slurm scripts (y/n)', default='n')
    opt_group.add_argument(
        '-snp', type=bool, help='set to True if features are SNPs', default=False)
    opt_group.add_argument(
        '-orf', type=bool, help='set to True if features are ORF presence/absence', default=False)
    opt_group.add_argument(
        '-cno', type=bool, help='set to True if features are copy number', default=False)
    opt_group.add_argument(
        '-submit', type=str, help='submit slurm scripts (y/n)', default='n')
    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()  # Read arguments

    # Arguments
    path1 = args.path1
    path2 = args.path2
    path3 = args.path3
    save = args.save
    suffix = args.suffix
    start = args.start
    stop = args.stop
    step = args.step
    batch = args.batch
    base = args.base

    # Write slurm scripts
    traits = ['YPACETATE', 'YPDCAFEIN40', 'YPDHU', 'YPETHANOL', 'YPD14', 'YPDCAFEIN50', 'YPDKCL2M', 'YPGALACTOSE', 'YPD40', 'YPDCHX05', 'YPDLICL250MM', 'YPGLYCEROL', 'YPD42', 'YPDCHX1', 'YPDMV', 'YPRIBOSE', 'YPD6AU', 'YPDCUSO410MM', 'YPDNACL15M',
              'YPSORBITOL', 'YPDANISO10', 'YPDDMSO', 'YPDNACL1M', 'YPXYLOSE', 'YPDANISO20', 'YPDETOH', 'YPDNYSTATIN', 'YPDANISO50', 'YPDFLUCONAZOLE', 'YPDSDS', 'YPDBENOMYL200', 'YPDFORMAMIDE4', 'YPDSODIUMMETAARSENITE', 'YPDBENOMYL500', 'YPDFORMAMIDE5']
    for trait in traits:
        print(trait)
        # Submit the slurm scripts
        if (args.runFS == 'y'):
            name = make_FS_files(path1, path2, path3, save, suffix, trait, start,
                                 stop, step, batch, snp=args.snp, orf=args.orf, cno=args.cno)
            if (args.submit == 'y'):
                # to create feature selection files
                os.system(f'sbatch {name}')
        if (args.runExpFS == 'y'):
            name = make_exp_FS_files(path1, path2, path3, save, suffix, trait,
                                     start, stop, base, snp=args.snp, orf=args.orf, cno=args.cno)
            if (args.submit == 'y'):
                # to create feature selection files
                os.system(f'sbatch {name}')
        if (args.runRF == 'y'):
            name = make_RF_files(path1, path2, path3, save, suffix, trait, start,
                                 stop, step, batch, snp=args.snp, orf=args.orf, cno=args.cno)
            if (args.submit == 'y'):
                os.system(f'sbatch {name}')  # to run feature selection on RF


if __name__ == '__main__':
    main()
