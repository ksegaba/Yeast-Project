#!/opt/software/Python/3.6.4-foss-2018a/bin/python
"""
Description: Combine h2 narrow-sense heritability results for each condition/trait into one data frame
Author: Kenia Segura Abá
"""
import os
import pandas as pd
import numpy as np

def main():
	# Traits
	traits=("YPACETATE", "YPDCAFEIN40", "YPDHU", "YPETHANOL", "YPD14", "YPDCAFEIN50", "YPDKCL2M", "YPGALACTOSE", "YPD40", "YPDCHX05", "YPDLICL250MM", "YPGLYCEROL", "YPD42", "YPDCHX1", "YPDMV", "YPRIBOSE", "YPD6AU", "YPDCUSO410MM", "YPDNACL15M", "YPSORBITOL", "YPDANISO10", "YPDDMSO", "YPDNACL1M", "YPXYLOSE", "YPDANISO20", "YPDETOH", "YPDNYSTATIN", "YPDANISO50", "YPDFLUCONAZOLE", "YPDSDS", "YPDBENOMYL200", "YPDFORMAMIDE4", "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDFORMAMIDE5")
	
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

	# Create empty lists to hold data
	h2 = [] # collect narrow-sense heritability values
	CI_upper = [] # collect upper confidence interval
	CI_lower = [] # collect lower confidence interval
	va = [] # collect additive variance
	ve = [] # collect residual/environmental variance
	loglik = [] # collect log-likelihood
	
	# Loop through h2_heritability_... files for each trait
	for trait in traits:
		for file in os.listdir("/mnt/home/seguraab/Shiu_Lab/Project"):
			if file.startswith("h2_heritability_"+trait):
				df = pd.read_csv(file, index_col=None, header=0)
				h2.append(df.loc[:,"h2"][0])
				CI_upper.append(df.loc[:,"conf.int1"][0])
				CI_lower.append(df.loc[:,"conf.int1"][1])
				va.append(df.loc[:,"va"][0])
				ve.append(df.loc[:,"ve"][0])
				loglik.append(df.loc[:,"loglik"][0])

	# Create a directory to hold data	
	data = {"Condition": traits, "h2":h2, "CI_upper":CI_upper, "CI_lower":CI_lower, "Va":va, "Ve":ve, "Log-likelihood":loglik}

	# Write data into a file
	data = pd.DataFrame(data)
	data.to_csv("h2_vs_Conditions.csv", header=True, index=True)