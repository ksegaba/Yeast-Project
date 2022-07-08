#!/opt/software/Python/3.6.4-foss-2018a/bin/python
"""
Description: Automated feature selection for all 35 conditions to 
assess the contribution of additive genetic variance to the phenotypic variance
Input: rrBLUP R-squared (R2) test results for each set of selected markers for all 35 conditions
Output: Graph of phenotype variation vs each condition
Author: Kenia Segura Abá
"""
import os
import pandas as pd
import numpy as np
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
    "YPDANISO10":"YPD Anisomycin 50µg/ml", "YPDNACL1M":"YPD NaCl 1M", "YPXYLOSE":"YP Xylose 2%", \
    "YPDANISO20":"YPD Anisomycin 20µg/ml", "YPDETOH":"YPD Ethanol 15%", "YPDSDS":"YPD SDS 0.2%", \
    "YPDSODIUMMETAARSENITE":"YPD Sodium metaarsenite 2.5mM", "YPDNYSTATIN":"YPD Nystatin 10µg/ml", \
    "YPDFLUCONAZOLE":"YPD Fluconazole 20µg/ml", "YPACETATE":"YP Acetate 2%", "YPDCAFEIN40":"YPD Caffeine 40mM", \
    "YPDHU":"YPD Hydroxyurea 30mg/ml", "YPETHANOL":"YP Ethanol 2%", "YPD14":"YPD 14C", \
    "YPDCAFEIN50":"YPD Caffeine 50mM", "YPDDMSO":"YPD DMSO 6%", "YPDANISO50":"YPD Anisomycin 50µg/ml", \
    "YPDBENOMYL200":"YPD Benomyl 200µg/ml", "YPDFORMAMIDE4":"YPD Formamide 4%", \
    "YPDBENOMYL500":"YPD Benomyl 500µg/ml", "YPDFORMAMIDE5":"YPD Formamide 5%"}

    # Create an empty directory to contain summary statistics of R2 values from each trait
    summary_statistics = {}
    baseline_summary_statistics = {}

    # Loop through R2 test results files for each trait
    for trait in traits:
      # Collect data and file names for the trait
      data = [] # Test R2 results for each 5-fold cross validation repeat (10 times)
      colnames = [] # file names
      baseline_data = [] # baseline refers to using all 64456 markers
      baseline_colnames = []
      for file in os.listdir("/mnt/scratch/seguraab/yeast_project"):
        if file.startswith("R2_test_results_exome_genoMarkers-"+trait):
          data.append(pd.read_csv(file, index_col=None, header=0))
          colnames.append(file)
      # collect baseline test R2 data and file names when using all markers
      for file in os.listdir("/mnt/scratch/seguraab/yeast_project"):
        if file.startswith("R2_test_results_exome_genoMarkers-"+trait+"_top64456"): 
          baseline_data.append(pd.read_csv(file, index_col=None, header=0))
          baseline_colnames.append(file)
      # Concatenate all data into one dataframe
      DF = pd.concat(data, axis=1, ignore_index=True)
      DF_baseline = pd.concat(baseline_data, axis=1, ignore_index=True)
      # Replace column names with file names
      DF.columns = colnames
      DF.columns = DF.columns.str.lstrip("R2_test_results_exome_genoMarkers-") # strip prefix from column names 
      DF.columns = DF.columns.str.rstrip("_.txt.csv") # strip suffix from column names 
      # Check output: for col in DF.columns: print(col)
      DF_baseline.columns = baseline_colnames
      DF_baseline.columns = DF_baseline.columns.str.lstrip("R2_test_results_exome_genoMarkers-")
      DF_baseline.columns = DF_baseline.columns.str.rstrip("_.txt.csv")
      # Feature selection: Calculate summary statistics for the trait and obtain column containing the max median R2
      summary_stats = DF.describe()
      summary_statistics[summary_stats.loc["50%"].idxmax()] = summary_stats[summary_stats.loc["50%"].idxmax()]
      # summary stats for baseline test R2
      baseline_summary_stats = DF_baseline.describe()
      baseline_summary_statistics[baseline_summary_stats.loc["50%"].idxmax()] = baseline_summary_stats[baseline_summary_stats.loc["50%"].idxmax()]
      
    # write data to file
    summary_statistics = pd.DataFrame(summary_statistics)
    summary_statistics.to_csv("data_summary_test_R2_vs_Conditions.csv", header=True, index=True)
    baseline_summary_statistics = pd.DataFrame(baseline_summary_statistics)
    baseline_summary_statistics.to_csv("data_baseline_summary_test_R2_vs_Conditions.csv", header=True, index=True)

    """
    # Create bar plot
    fig, ax = plt.subplots(figsize=(15,15)) # initialize plot
    index = np.arange(35) # position of bars
    rects1 = ax.bar(index, mean_vals, width=0.5, yerr=std_vals, edgecolor="black", color="cyan", label="Test R2") # plot data after feature selection
    rects2 = ax.bar(index+0.5, mean_vals_baseline, width=0.5, yerr=std_vals_baseline, edgecolor="black", color="yellow", label="Baseline Test R2") # plot baseline data
    ax.tick_params(axis='x', which='major', labelsize=8) # change tick label font size
    ax.set_xlabel("Condition") # add x axis label
    ax.set_ylabel("rrBLUP Test R2", fontsize=26) # add y axis label
    ax.set_title("rrBLUP Test R2 in each condition", fontsize=26) # add title
    ax.legend() # add legend
    ax.legend(fontsize=26) # change legend font size
    ax.set_xticks(index+0.25) # x axis label position
    ax.set_xticklabels(traits, rotation=90) # add x axis tick labels
    plt.savefig("test_R2_vs_Conditions.png") # save plot
    """
if __name__ == "__main__":
	main()
