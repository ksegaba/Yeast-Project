import sys,os,argparse

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for calculating the PCC for a gene expression matrix')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-path', help='path to the work directory', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	
	path = args.path
	for files in os.listdir('./'):
		if files.endswith('_1.fastq'):
			name = files.split('_1.fastq')[0]
			out = open('Trim_%s.sh'%(name),'w')
			out.write('#!/bin/sh --login\n#SBATCH --time=4:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=20G\n#SBATCH --job-name Trim_%s.sh\n#SBATCH -e Trim_%s.sh.e\n#SBATCH -o Trim_%s.sh.o\ncd %s\n'%(name,name,name,path))
			out.write('java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE %s_1.fastq %s_2.fastq %s_1trimmed.fastq %s_2trimmed.fastq %s_1Utrimmed.fastq %s_2Utrimmed.fastq ILLUMINACLIP:$/mnt/home/peipeiw/Documents/Yeast_Kenia/PE_adapters.txt:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30'%(name,name,name,name,name,name))
			out.close()
			os.system('sbatch Trim_%s.sh'%name)

if __name__ == '__main__':
	main()