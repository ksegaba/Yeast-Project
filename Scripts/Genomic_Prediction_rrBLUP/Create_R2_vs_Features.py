#!/opt/software/Python/3.6.4-foss-2018a/bin/python
"""
Description: Automated feature selection for all 35 conditions to 
assess the contribution of additive genetic variance to the phenotypic variance
Input: rrBLUP R2 test results for each set of selected markers for all 35 conditions
Output: Graph of phenotype variation vs each condition
Author: Kenia Segura Abá
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

def main():
   # traits
   traits=("YPACETATE", "YPDCAFEIN40", "YPDHU", "YPETHANOL", "YPD14", "YPDCAFEIN50", "YPDKCL2M", "YPGALACTOSE", "YPD40", "YPDCHX05", "YPDLICL250MM", "YPGLYCEROL", "YPD42", "YPDCHX1", "YPDMV", "YPRIBOSE", "YPD6AU", "YPDCUSO410MM", "YPDNACL15M", "YPSORBITOL", "YPDANISO10", "YPDDMSO", "YPDNACL1M", "YPXYLOSE", "YPDANISO20", "YPDETOH", "YPDNYSTATIN", "YPDANISO50", "YPDFLUCONAZOLE", "YPDSDS", "YPDBENOMYL200", "YPDFORMAMIDE4", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDFORMAMIDE5")
   
   # Growth Conditions and Descriptions
   # growth temperature is 30C if not specified
   cond = {"YPDKCL2M":"YPD KCL 2M", "YPGALACTOSE":"YP Galactose 2%", "YPD40":"YPD 40C", \
   "YPDCHX05":"YPD Cycloheximide 0.5µg/ml", "YPDLICL250MM":"YPD LiCl 250mM", \
   "YPGLYCEROL":"YP Glycerol 2%", "YPD42":"YPD 42C", "YPDCHX1":"YPD Cycloheximide 1µg/ml", \
   "YPDMV":"YPD Methylviologen 20mM", "YPRIBOSE":"YP Ribose 2%", "YPD6AU":"YPD 6-Azauracile 600µg/ml", \
   "YPDCUSO410MM":"YPD CuSO4 10mM", "YPDNACL15M":"YPD NaCl 1.5M", "YPSORBITOL":"YP Sorbitol 2%", \
   "YPDANISO10":"YPD Anisomycin 10µg/ml", "YPDNACL1M":"YPD NaCl 1M", "YPXYLOSE":"YP Xylose 2%", \
   "YPDANISO20":"YPD Anisomycin 20µg/ml", "YPDETOH":"YPD Ethanol 15%", "YPDSDS":"YPD SDS 0.2%", \
   "YPDSODIUMMETAARSENITE":"YPD Sodium metaarsenite 2.5mM", "YPDNYSTATIN":"YPD Nystatin 10µg/ml", \
   "YPDFLUCONAZOLE":"YPD Fluconazole 20µg/ml", "YPACETATE":"YP Acetate 2%", "YPDCAFEIN40":"YPD Caffeine 40mM", \
   "YPDHU":"YPD Hydroxyurea 30mg/ml", "YPETHANOL":"YP Ethanol 2%", "YPD14":"YPD 14C", \
   "YPDCAFEIN50":"YPD Caffeine 50mM", "YPDDMSO":"YPD DMSO 6%", "YPDANISO50":"YPD Anisomycin 50µg/ml", \
   "YPDBENOMYL200":"YPD Benomyl 200µg/ml", "YPDFORMAMIDE4":"YPD Formamide 4%", \
   "YPDBENOMYL500":"YPD Benomyl 500µg/ml", "YPDFORMAMIDE5":"YPD Formamide 5%"}

   ##### First, we will be looking at the cross validation results #####
   # Loop through all files in a directory for each trait
   for trait in traits:
      # Create an empty list to contain data and file names
      df = []
      colnames = []
      for file in os.listdir("/mnt/scratch/seguraab/yeast_project/yeast_rrBLUP_results"):
         if file.startswith("R2_cv_results_rrBLUP_geno_Markers_"):
            if file.endswith(trait+".csv"):
               df.append(pd.read_csv(file, index_col=None, header=0))
               colnames.append(file)
      # Concatenate all data into one dataframe
      DF = pd.concat(df, axis=1, ignore_index=True)
      # Replace column names with file names
      DF.columns = colnames
      DF.columns = DF.columns.str.lstrip("1234567890_-abcdefghijklmnoqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") # remove any of these characters from the left side
      DF.columns = DF.columns.str.lstrip("p") # to prevent removal of # of features after "top"
      DF.columns = DF.columns.str.rstrip("_"+trait+".csv") # strip suffix from column names 
      # Check output: for col in DF.columns: print(col)
      # Sort the columns by column name
      DF.columns.map(type) # at this point, names are strings
      DF.columns = DF.columns.map(int) # convert names to integers
      DF = DF.sort_index(axis=1) # sort dataframe by column names
      # write dataframe into a file
      DF.to_csv("rrBLUP_FS_cv_R2_results_"+trait+".csv", header=True, index=False)
      """
      # Plot R2 vs # of Features
      av = DF.mean(axis=0) # average for each column
      se = DF.sem(axis=0) # standard error for each column
      sd = DF.std(axis=0) # standard deviation for each column
      fig = plt.figure(figsize=(30,5)) # initialize plot
      plt.plot(DF.columns, av) # plot lines
      plt.errorbar(DF.columns, av, yerr=se, fmt="-o") # add error bars w/points
      plt.xticks(DF.columns, rotation=45) # add all tick marks
      plt.title("CV R2 vs # of Features") # add title
      plt.savefig("cv_R2_vs_Features.png") # save plot
      plt.close
      """
      # Write summary statistics to a file
      stats = DF.describe()
      stats.to_csv("rrBLUP_FS_cv_R2_stats_"+trait+".csv", header=True, index=True)

   ###############################################
   ##### Finally, we will look at the test results #####
   # Loop through all files in a directory for each trait
   for trait in traits:
      # Create an empty list to contain data and file names
      df = []
      colnames = []
      for file in os.listdir("/mnt/scratch/seguraab/yeast_project"):
         if file.startswith("R2_test_results_exome_genoMarkers-"+trait):
            df.append(pd.read_csv(file, index_col=None, header=0))
            colnames.append(file)
      # Concatenate all data into one dataframe
      DF = pd.concat(df, axis=1, ignore_index=True)
      # Replace column names with file names
      DF.columns = colnames
      DF.columns = DF.columns.str.lstrip("1234567890_-abcdefghijklmnorstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") # remove any of these characters from the left side
      DF.columns = DF.columns.str.lstrip("p") # to prevent removal of # of features after "top"
      DF.columns = DF.columns.str.rstrip("_.txt.csv") # strip suffix from column names 
      # Check output: for col in DF.columns: print(col)
      # Sort the columns by column name
      DF.columns.map(type) # at this point, names are strings
      DF.columns = DF.columns.map(int) # convert names to integers
      DF = DF.sort_index(axis=1) # sort dataframe by column names
      # write dataframe into a file
      DF.to_csv("Feature_selection_test_"+trait+".csv", header=True, index=False)
      """
      # Plot R2 vs # of Features
      av = DF.mean(axis=0) # average for each column
      se = DF.sem(axis=0) # standard error for each column
      sd = DF.std(axis=0) # standard deviation for each column
      fig = plt.figure(figsize=(30,5)) # initialize plot
      plt.plot(DF.columns, av) # plot lines
      plt.errorbar(DF.columns, av, yerr=se, fmt="-o") # add error bars w/points
      plt.xticks(DF.columns, rotation=45) # add all tick marks
      plt.title("Test R2 vs # of Features") # add title
      plt.savefig("test_R2_vs_Features.png") # save plot
      plt.close
      """
      # Write summary statistics to a file
      stats = DF.describe()
      stats.to_csv("Feature_selection_stats_test_"+trait+".csv", header=True, index=True)

if __name__ == "__main__":
	main()
