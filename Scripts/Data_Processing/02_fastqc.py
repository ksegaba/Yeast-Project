import sys,os
for files in os.listdir('./'):
	if files.endswith('fastq'):
		out = open('Fastqc_%s.sh'%files,'w')
		out.write('#!/bin/sh --login\n#SBATCH --time=1:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=40G\n#SBATCH --job-name Fastqc_%s.sh\n#SBATCH -e Fastqc_%s.sh.e\n#SBATCH -o Fastqc_%s.sh.o\ncd /mnt/scratch/peipeiw/Allelic_specific/00_rawdata/\n'%(files,files,files))
		out.write('module load Python\nmodule load GCC/7.3.0-2.30  OpenMPI/3.1.1\nmodule load Bowtie2/2.3.4.2\n')
		out.write('fastqc -f fastq %s'%files)
		out.close()
