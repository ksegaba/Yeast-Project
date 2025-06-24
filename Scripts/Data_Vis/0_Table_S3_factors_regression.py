#!/usr/bin/env python3

import os
import pickle
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score
from glob import glob

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
### TABLE S3
################################################################################
### linear regression of single-env model performance on fitness-related factors 
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)  # fitness data
var = pheno.var(axis=0) ; var.name = "var"
med = pheno.median(axis=0) ; med.name = "median"

# how many environments have small variance?
cv = (pheno.std(axis=0) / pheno.mean(axis=0) * 100).sort_values() # coefficient of variation (%)
cv.describe()

h2 = pd.read_csv("Data/Peter_2018/Heritability_h2_H2_sommer_CORRECTED.csv") # trait heritabilities
h2.set_index("Conditions", inplace=True)
X = pd.concat([var, med, h2.h2], ignore_index=False, axis=1) # fitness factors
X_s = (X-X.mean())/X.std() # center and scale

for data_type in ["SNPs", "PAVs", "CNVs", "PCs"]:
	if data_type == "PCs": # model performance results
		Y = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")
	else:
		Y = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_%s_FS.txt" % data_type, sep="\t")
	
	Y.set_index("Y", inplace=True)
	
	## Regress model performance on all factors at once
	df = pd.concat([X_s, Y.r2_test], axis=1)
	mod = smf.ols(formula=f"r2_test ~ {' * '.join(df.columns[:3])}", data=df)
	res = mod.fit()
	yhats = res.predict()
	r2_scores = r2_score(df.r2_test, yhats, multioutput=None) # double check sm.ols R2
	
	with open(f"Scripts/Data_Vis/Section_3/Table_S3_factors_ols_{data_type}_results_CORRECTED.txt", "w") as out:
		out.write(res.summary().as_text())
		out.write(f"\nR-sq: {r2_scores}")
	# vars(res) # attributes
	
	pickle.dump(mod, open(f"Scripts/Data_Vis/Section_3/Table_S3_factors_ols_{data_type}_model_CORRECTED.pkl", 'wb')) # save the model
	
	yhats=pd.Series(yhats)
	yhats.index = df.index
	yhats.name = 'y_pred'
	pd.concat([Y.r2_test, yhats], ignore_index=False, axis=1).\
		to_csv(f"Scripts/Data_Vis/Section_3/Table_S3_factors_ols_{data_type}_preds_CORRECTED.csv")
	
	del df, Y, mod

# Regression model performance on the number of features
for data_type in ["SNPs", "PAVs", "CNVs", "PCs"]:
	if data_type == "PCs": # model performance results
		Y = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")
	else:
		Y = pd.read_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_%s_FS.txt" % data_type, sep="\t")
	
	Y.set_index("Y", inplace=True)
	
	mod = smf.ols(formula=f"r2_test ~ FeatureNum", data=Y)
	res = mod.fit()
	pickle.dump(mod, open(f"Scripts/Data_Vis/Section_3/Table_S3_featnum_ols_{data_type}_model_CORRECTED.pkl", 'wb'))
	with open(f"Scripts/Data_Vis/Section_3/Table_S3_featnum_ols_{data_type}_results_CORRECTED.txt", "w") as out:
		out.write(res.summary().as_text())

