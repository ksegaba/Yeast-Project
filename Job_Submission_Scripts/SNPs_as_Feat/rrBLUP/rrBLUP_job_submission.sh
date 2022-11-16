#!/bin/bash --login
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # 4
#SBATCH --cpus-per-task=1 # 4
#SBATCH --mem-per-cpu=36G
#SBATCH --job-name step11_rerun_0192022
#SBATCH --output=%x_%j

cd /mnt/scratch/seguraab/yeast_project/yeast_rrBLUP_results

module purge 

#echo "eleven: create test set"
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
#python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_rrBLUP/10_holdout_test_stratified.py pheno.csv YPACETATE 6

#echo "twelve: build rrBLUP model using training set"
#echo "...split geno and pheno into training and testing sets"
#module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_rrBLUP/11_split_geno_pheno_fread.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt

#echo "...generate the cross-validation scheme file using the pheno training data"
#module purge
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
#python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_rrBLUP/07_make_CVs.py -file /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno_training.csv -cv 5 -number 10

#echo "...build rrBLUP model using training geno and pheno data"
#module purge
#module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2
#~~~~~~~# Note: re-ran this on 01/19/2022
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_training.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno_training.csv all all 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv rrBLUP_geno
Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/09_rrBLUP_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_training.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno_training.csv all ${trait} 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv rrBLUP_geno



#echo "step 13: create Marker files for feature selection"
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
# the marker files do not change regardless of trait ## 08/16/22 this doesn't make sense. the markers are sorted by coef. how could it be that each trait puts the markers in the same order?
#~~~~~~~# Note: re-ran this on 01/19/2022
#python /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/12_select_markers_according_to_abs_coef.py -coef Coef_rrBLUP_geno_YPACETATE.csv -start 500 -stop 50500 -step 500




echo "step 14a: obtain test R2 values for each condition and when using a different number of markers (ranging from 500 to 50500)"
TRAITS=( YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDDMSO YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDNYSTATIN YPDANISO50 YPDFLUCONAZOLE YPDSDS YPDBENOMYL200 YPDFORMAMIDE4 YPDSODIUMMETAARSENITE YPDBENOMYL500 YPDFORMAMIDE5 )
module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2
#for trait in ${TRAITS[*]} # for each trait (condition)
#do
#echo "${trait} baseline rrBLUP"
#~~~~~~~# Note: re-ran this on 01/19/2022
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/13_rrBLUP_training_test_split_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all $trait /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv rrBLUP_geno
#done

FILES=$(ls Markers-*)
for trait in ${TRAITS[*]}
do
for file in $FILES # each marker file corresponding to the trait
do
echo "${file} feature selection rrBLUP"
Rscript /mnt/home/seguraab/Shiu_Lab/Project/External_software/Genomic_prediction_in_Switchgrass/13_rrBLUP_training_test_split_fread_predict_values.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv $file $trait /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt 5 10 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv rrBLUP_geno_${file} 
done
done

#echo "step14a-2: Create feature selection R2 plots"
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
#python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_rrBLUP/Create_R2_vs_Features.py

#echo "step14a-3: Create R2 vs Conditions plots"
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4
#python /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_rrBLUP/Create_R2_vs_Conditions.py

#echo "step 14b: calculate narrow-sense heritability for each condition after feature selection"
#TRAITS=( YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDDMSO YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDNYSTATIN YPDANISO50 YPDFLUCONAZOLE YPDSDS YPDBENOMYL200 YPDFORMAMIDE4 YPDSODIUMMETAARSENITE YPDBENOMYL500 YPDFORMAMIDE5 )
#FILES=( YPACETATE_top1000 YPDCAFEIN40_top10000 YPDHU_top750 YPETHANOL_top9750 YPD14_top8950 YPDCAFEIN50_top64456 YPDKCL2M_top7950 YPGALACTOSE_top750 YPD40_top2750 YPDCHX05_top8650 YPDLICL250MM_top8550 YPGLYCEROL_top4750 YPD42_top9650 YPDCHX1_top10000 YPDMV_top9050 YPRIBOSE_top8650 YPD6AU_top8250 YPDCUSO410MM_top4000 YPDNACL15M_top1250 YPSORBITOL_top1500 YPDANISO10_top9850 YPDDMSO_top8050 YPDNACL1M_top4000 YPXYLOSE_top9650 YPDANISO20_top3250 YPDETOH_top10000 YPDNYSTATIN_top9950 YPDANISO50_top8350 YPDFLUCONAZOLE_top1750 YPDSDS_top9450 YPDBENOMYL200_top7750 YPDFORMAMIDE4_top10000 YPDSODIUMMETAARSENITE_top5950 YPDBENOMYL500_top9750 YPDFORMAMIDE5_top250 )
#module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2
#n=0 # counter to loop through FILES
#for trait in ${TRAITS[*]} # for each condition
#do
#echo "${trait}"
#echo "Markers-${FILES[n]}_.txt"
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/Genomic_Prediction_rrBLUP/calculate_heritability.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv all ${trait}
#Rscript /mnt/home/seguraab/Shiu_Lab/Project/Genomic_Prediction_rrBLUP/calculate_heritability.r /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv Markers-${FILES[n]}_.txt ${trait}
#((n+=1)) # increment counter
#done

# step 15: obtain y_pred for each condition using the "optimum" # of features to make figure of y_pred vs y_actual for the test set 
TRAITS=( YPACETATE YPDCAFEIN40 YPDHU YPETHANOL YPD14 YPDCAFEIN50 YPDKCL2M YPGALACTOSE YPD40 YPDCHX05 YPDLICL250MM YPGLYCEROL YPD42 YPDCHX1 YPDMV YPRIBOSE YPD6AU YPDCUSO410MM YPDNACL15M YPSORBITOL YPDANISO10 YPDDMSO YPDNACL1M YPXYLOSE YPDANISO20 YPDETOH YPDNYSTATIN YPDANISO50 YPDFLUCONAZOLE YPDSDS YPDBENOMYL200 YPDFORMAMIDE4 YPDSODIUMMETAARSENITE YPDBENOMYL500 YPDFORMAMIDE5 )
FILES=(YPACETATE_top1000 YPDCAFEIN40_top10000 YPDHU_top750 YPETHANOL_top9750 YPD14_top8950 YPDCAFEIN50_top64456 YPDKCL2M_top7950 YPGALACTOSE_top750 YPD40_top2750 YPDCHX05_top8650 YPDLICL250MM_top8550 YPGLYCEROL_top4750 YPD42_top9650 YPDCHX1_top10000 YPDMV_top9050 YPRIBOSE_top8650 YPD6AU_top8250 YPDCUSO410MM_top4000 YPDNACL15M_top1250 YPSORBITOL_top1500 YPDANISO10_top9850 YPDDMSO_top8050 YPDNACL1M_top4000 YPXYLOSE_top9650 YPDANISO20_top3250 YPDETOH_top10000 YPDNYSTATIN_top9950 YPDANISO50_top8350 YPDFLUCONAZOLE_top1750 YPDSDS_top9450 YPDBENOMYL200_top7750 YPDFORMAMIDE4_top10000 YPDSODIUMMETAARSENITE_top5950 YPDBENOMYL500_top9750 YPDFORMAMIDE5_top250 )
module load GCC/8.3.0  OpenMPI/3.1.4
module load R/4.0.2
n=0
for trait in ${TRAITS[*]} # for each trait (condition)
do
echo $n
echo $trait
echo ${FILES[n]} # for each marker file corresponding to that trait

((n+=1))
done

# echo "step 15: sEERBLUP on geno.csv for each condition"

scontrol show job $SLURM_JOB_ID







