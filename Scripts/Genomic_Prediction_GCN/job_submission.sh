#!/bin/bash --login
#SBATCH --time=40:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=36G
#SBATCH --job-name sequential_model_initial
#SBATCH --output=%x_%j

cd /mnt/scratch/seguraab/yeast_project/ 

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load CUDA/10.0.130 cuDNN/7.5.0.56-CUDA-10.0.130
module load Python/3.6.4
source ~/tf-1.13.1-env/bin/activate

python /mnt/home/seguraab/Shiu_Lab/Project/DL_Genomic_Prediction/sequential_model.py




