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

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")

################################################################################
### TABLE S3
################################################################################
### linear regression of single-env model performance on fitness-related factors 
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)  # fitness data
var = pheno.var(axis=0) ; var.name = "var"
med = pheno.median(axis=0) ; med.name = "median"

h2 = pd.read_csv("Scripts/Data_Vis/Section_2/Heritability_h2_H2_sommer.csv") # trait heritabilities
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
	with open(f"Scripts/Data_Vis/Section_3/factors_ols_{data_type}_results.txt", "w") as out:
		out.write(res.summary().as_text())
		out.write(f"\nR-sq: {r2_scores}")
	# vars(res) # attributes
	pickle.dump(mod, open(f"Scripts/Data_Vis/Section_3/factors_ols_{data_type}_model.pkl", 'wb')) # save the model
	yhats=pd.Series(yhats)
	yhats.index = df.index
	yhats.name = 'y_pred'
	pd.concat([Y.r2_test, yhats], ignore_index=False, axis=1).\
		to_csv(f"Scripts/Data_Vis/Section_3/factors_ols_{data_type}_preds.csv")
	del df, Y, mod

################### GxE model performance linear regression ####################
# model performance regressed on degree of correlation among the correlated environments
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)  # fitness data
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)  # fitness data
pheno_test = pheno.loc[pheno.index.isin(test[0]),:]

## CV1 Multi-trait models
d = "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results"
to_regress = pd.DataFrame()
for env in pheno_test.columns:
	print(env)
	# Read in model predictions
	file = glob(f"{d}/multi-trait/cv1/GxE_cv1_pav_5_envs_{env}_*_preds_test_multitrait.txt")[0]
	preds = pd.read_csv(file, sep="\t")
	# Calculate test R2 values
	res = pd.concat([preds[["Rep", "Trait"]], preds.apply(lambda x: \
		r2_score(pheno_test[x.Trait], x.iloc[2:]), axis=1)], axis=1)
	res.rename(columns={0:"r2_test"}, inplace=True)
	res_summary = pd.DataFrame()
	res_summary["r2_test"] = res.groupby(["Trait"])["r2_test"].mean()
	res_summary["r2_test_sd"] = res.groupby(["Trait"])["r2_test"].std()
	res_summary["r2_test_se"] = res.groupby("Trait")["r2_test"].sem()
	# Get degree of correlation between the traits
	cor = pheno_test.loc[:,preds.Trait.unique()].corr(method="pearson")
	# keep = np.triu(np.ones(cor.shape), k=1).astype('bool') # the r values need to be sorted properly, hence, see below
	# to_regress = pd.concat([to_regress, pd.Series(res_summary.r2_test.to_list() \
	# 	+ cor.where(keep).stack().to_list())], axis=1)
	E1_cors = cor.loc[env,:].sort_values(ascending=False)[1:]
	E2_cors = cor.loc[E1_cors.index[0],:][E1_cors.index[1:]].sort_values(ascending=False)
	E3_cors = cor.loc[E1_cors.index[1],:][E1_cors.index[2:]].sort_values(ascending=False)
	to_regress = pd.concat([to_regress, pd.concat([res_summary.r2_test[[env]+list(E1_cors.index.values)],
		E1_cors.reset_index().iloc[:,1], E2_cors.reset_index().iloc[:,1],
		E3_cors.reset_index().iloc[:,1], pd.Series(cor.loc[E1_cors.index[2],\
		E1_cors.index[3]])], axis=0, ignore_index=True)], axis=1)
	
to_regress.columns = pheno_test.columns
to_regress.index = ["E1_r2_test", "E2_r2_test", "E3_r2_test",
	"E4_r2_test", "E5_r2_test", "E1_E2_cor", "E1_E3_cor", "E1_E4_cor",
	"E1_E5_cor", "E2_E3_cor", "E2_E4_cor", "E2_E5_cor", "E3_E4_cor",
	"E3_E5_cor", "E4_E5_cor"]
to_regress = to_regress.T
to_regress.to_csv("Scripts/Data_Vis/Section_2/RESULTS_RF_GxE_MULTI-TRAIT.txt", sep="\t")

# Factors regression on GxE model performance
mod = smf.ols(formula=f"E1_r2_test ~ {' + '.join(to_regress.columns[5:9])}", data=to_regress)
res = mod.fit()
yhats = res.predict()
r2_scores = r2_score(to_regress.E1_r2_test, yhats, multioutput=None) # double check sm.ols R2
with open(f"Scripts/Data_Vis/Section_3/factors_ols_multi-trait_results.txt", "w") as out:
	out.write(res.summary().as_text())
	out.write(f"\nR-sq: {r2_scores}")

# Plot E1_r2_test and E1_E3_cor since it has the significant p-value
sns.regplot(data=to_regress, x="E1_E3_cor", y="E1_r2_test")
plt.xlabel("Correlation of target env with second-most correlated env")
plt.ylabel("Multi-trait R2_test")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_3/factors_ols_multi-trait_plot.pdf")
plt.close()

