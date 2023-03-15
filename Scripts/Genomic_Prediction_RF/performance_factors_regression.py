"""
Linear regression of random forest model performances in each environment on
factors potentially influencing predictions based on fitness data.
"""
__author__ = "Kenia Segura Abá"

import os
import sys
import argparse
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.tree import DecisionTreeRegressor

os.chdir("/mnt/home/seguraab/Shiu_Lab/Project")  # Set the working directory


def warn(*args, **kwargs):
	pass


warnings.warn = warn


def kmeans_clustering(to_cluster, start, stop):
	distortions = []
	inertias = []
	mapping1 = {}
	mapping2 = {}
	for k in range(start,stop):
		clusterer = KMeans(n_clusters=k, random_state=10).fit(to_cluster)
		distortions.append(sum(
			np.min(cdist(
				to_cluster, clusterer.cluster_centers_, "euclidean"), 
				axis=1))/to_cluster.shape[0])
		inertias.append(clusterer.inertia_)
		mapping1[k] = sum(np.min(cdist(
			to_cluster, clusterer.cluster_centers_, "euclidean"), 
			axis=1)) / to_cluster.shape[0]
		mapping2[k] = clusterer.inertia_
		# Check clusters
		to_cluster = to_cluster.sort_values("median")
		colors = plt.cm.get_cmap("hsv", k)
		fig, ax1 = plt.subplots(1, 1)
		ax1.scatter(to_cluster.index, to_cluster['median'], c=clusterer.labels_, cmap=colors)
		ax1.set_ylabel("Median fitness")
		plt.xticks(rotation=50, ha='right', rotation_mode='anchor')
		plt.subplots_adjust(bottom=0.1)
		plt.savefig(f"Scripts/Genomic_Prediction_RF/fitness_severity_k{k}.pdf")
		plt.close()
		if k==10: # based on elbow plots
			cluster_labels = 1 - ((clusterer.labels_+1)/k) # get cluster labels to return
	# Elbow subplots
	fig, (ax1,ax2) = plt.subplots(1, 2)
	fig.set_size_inches(10, 5)
	ax1.plot(range(start,stop), distortions, 'bx-')
	ax1.set_xlabel('Values of K')
	ax1.set_ylabel('Distortion')
	ax1.set_title('The Elbow Method using Distortion')
	ax2.plot(range(start,stop), inertias, 'bx-')
	ax2.set_xlabel('Values of K')
	ax2.set_ylabel('Distortion')
	ax2.set_title('The Elbow Method using Inertia')
	plt.savefig(f"Scripts/Genomic_Prediction_RF/factors_kmeans_elbow.pdf")
	return to_cluster, cluster_labels


def process_data(df):
	""" Pre-process the fitness data to generate a factor table """
	cor = df.corr(method="pearson") # fitness environment correlations
	cor = cor.add_prefix(prefix="cor_")
	cov = df.cov() # fitness environment covariance
	cov = cov.add_prefix(prefix="cov_")
	var = df.var(axis=0) # variance in fitness per environment
	var.name = "fitness_variance"
	med = df.median(axis=0) # median fitness per environment
	med.name = "median_fitness"
	to_cluster = df.describe().transpose()
	to_cluster.insert(3, "median", med)
	to_cluster.insert(1, "var", var)
	# Assign environments to a level of severity by K-means clustering
	to_cluster, cluster_labels = kmeans_clustering(to_cluster, 2, 15)
	# make final feature table
	to_cluster.insert(10, "severity", cluster_labels) # level of stress severity
	to_cluster.drop("count", axis=1, inplace=True)
	X = pd.concat([cor, cov, to_cluster], axis=1) # feature table
	return X


def linear_mod(X, Y, n, prefix):
	""" Ordinary Least Squares Linear Regression """
	coefs = pd.DataFrame()
	yhats = pd.DataFrame()
	r2_scores = []
	for i in range(n):
		mod = LinearRegression().fit(X, Y)
		pickle.dump(mod, open(f"{prefix}_linear_model_factors_rep{n}.pkl", 'wb'))
		coefs = pd.concat([coefs, pd.Series(mod.coef_)])
		# calculate variation explained by all the factors
		yhat = np.dot(X, coefs) + mod.intercept_
		yhats = pd.concat([yhats, pd.Series(yhat)])
		r2_scores.append(r2_score(Y, yhat))
	return mod, coefs, yhats, r2_scores

def hyperparameter_tuning(parameters):
	
	return best_mod

def decisionTree_mod(X, Y):
	coefs = pd.DataFrame()
	yhats = pd.DataFrame()
	r2_scores = []
	for i in range(n):
		regressor = DecisionTreeRegressor(random_state=0)
	return mod, coefs


if __name__ == "__main__":
	# Argument parser
	parser = argparse.ArgumentParser(description="")
	req_group = parser.add_argument_group(title="Required Input")
	req_group.add_argument(
		"-Y", type=str, help="Path to CSV file of model performances")
	opt_group = parser.add_argument_group(title="Optional Input")
	opt_group.add_argument(
		"-sep", type=str, help="Delimeter for file Y", default=",")
	opt_group.add_argument(
		"-index_col", type=str, help="Name of column in Y to set as index", default="cond")
	opt_group.add_argument(
		"-y_col", type=str, help="Name of column in Y to set as dependent variable", default="r2_test")
	opt_group.add_argument(
		"-prefix", type=str, help="Prefix of file name to save models as. Prefix may include file path too.")

	if len(sys.argv) == 1:  # Help
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()  # Read arguments

	# Debugging arguments
	# Y="Results/RESULTS_RF_SNPs_FS.txt"
	# sep="\t"
	# index_col="cond"
	# y_col="r2_test"
	# prefix="Scripts/Genomic_Prediction_RF/snps_RF_FS"

	df = pd.read_csv("Data/Peter_2018/pheno.csv")  # Read in fitness data
	h2 = pd.read_csv("Results/Heritability_h2_H2_sommer.csv") # Read in trait heritabilities
	Y = pd.read_csv(args.Y, sep=args.sep)  # Read in model performance data

	h2.set_index("Conditions", inplace=True) # set environments as row index
	Y.set_index(args.index_col, inplace=True) # set environments as row index
	Y = Y.loc[:,args.y_col] # keep only dependent variable column
	X = process_data(df)  # make the feature table
	X = pd.concat([h2.h2, X], axis=1) # add heritabilities
	X.to_csv("Scripts/Genomic_Prediction_RF/")

	# Build models
	lmod, coefs = linear_mod(X, Y, args.prefix)  # fit the linear model
	dtmod, imps = decisionTree_mod(X, Y, args.prefix) # fit the decision tree model

	# Save model coefficients
	lm_df = pd.concat([pd.Series(X.columns), pd.Series(coefs)], axis=1)
	lm_df.columns = ["feature", "coef"]
	lm_df.sort_values("coef", ascending=False, inplace=True)
	lm_df.to_csv(f"{prefix}_linear_model_factor_coefs.csv", index=False)

	# write bullet points about why the different factors have higher coeffs, is there anything about the data 