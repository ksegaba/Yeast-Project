"""
Linear regression of random forest model performances in each environment on
factors potentially influencing predictions based on fitness data.
"""
__author__ = "Kenia Segura Abá"

import os
import sys
import argparse
import warnings
import pickle
import random
import time
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import explained_variance_score
from sklearn.metrics import r2_score

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
	plt.close()
	return to_cluster, cluster_labels


def process_data(df):
	""" Pre-process the fitness data to generate a factor table """
	# cor = df.corr(method="pearson") # fitness environment correlations
	# cor = cor.add_prefix(prefix="cor_")
	# cov = df.cov() # fitness environment covariance
	# cov = cov.add_prefix(prefix="cov_")
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
	# X = pd.concat([cor, cov, to_cluster], axis=1) # feature table
	# return X
	return to_cluster


def cumsum(evar):
	"""
	Get the cumulative sum of numerical items in a list
	Input: a list of values
	Returns: a list of cumulatively added values
	"""
	length=len(evar)
	cu_lst = [sum(evar[0:x]) for x in range(0, length+1)]
	return cu_lst


def pca(X):
	"""Principal Component Analysis for dimension reduction of feature table"""
	# scale and center the data
	X.boxplot(column=X.columns.to_list())
	plt.xticks(rotation=45, fontsize=6)
	plt.savefig("Scripts/Genomic_Prediction_RF/factor_box.pdf", height=5, width=15)
	plt.close()
	X_s = (X-X.mean())/X.std() # center and scale
	X_s.boxplot(column=X_s.columns.to_list())
	plt.xticks(rotation=45, fontsize=6)
	plt.savefig("Scripts/Genomic_Prediction_RF/factor_box_standardized.pdf", height=5, width=15)
	plt.close()
	# run PCA
	num_pc = min(X_s.shape[0]-1, X_s.shape[1]) # min(n-1, d)
	pca = PCA(n_components=num_pc)
	X_r = pca.fit(X_s).transform(X_s) # loading vectors? (PCs, eigenvectors of the square covariance matrix)
	evar = pca.explained_variance_ratio_.tolist() # explained variance of each PC
	cu_evar = cumsum(evar) # cumulative sum of explained variance
	# determine number of principal components to keep using scree plot
	plt.plot(evar, marker='o')
	plt.plot(cu_evar, marker='o')
	plt.xlabel("Principal Component")
	plt.ylabel("Proportion of Variance Explained")
	plt.title("Scree Plot")
	plt.savefig("Scripts/Genomic_Prediction_RF/scree_plot.pdf")
	plt.close()
	# Kaiser-Guttman criterion for choosing principal components (PCs)
	ev = list(abs(pca.singular_values_)) # eigenvalues (explained variance of each PC)
	nPCs = len([x for x in ev if x>=1])
	X_r = pd.DataFrame(X_r)
	X_r_sub = X_r.iloc[:,0:nPCs]
	X_r_sub.index = X.index
	# determine which features are linearly combined in each principal component
	# Each PC is a linear combination of the variables:
	# equation: PC = w1x1 + w2x2 + ... + wnxn
	# The weight each variable has are the loadings of the PC (eigenvector)
	V = pd.DataFrame(pca.components_, columns=X.columns) # singular vectors (multiply to X to get X_r)
	V = V.T
	V.columns = [f"PC{i}" for i in range(1, num_pc+1)]
	V.sort_values(by="PC1")
	V.to_csv("Scripts/Genomic_Prediction_RF/singular_vectors.csv")
	plt.scatter(X_r[0], X_r[1]) # plot the first two PCs
	plt.xlabel("PC 1")
	plt.ylabel("PC 2")
	for i in range(V.shape[0]): # plot the loading vectors
		plt.arrow(0, 0, V.iloc[i,0], V.iloc[i,1], color="r", alpha=0.5)
		# plt.text(V.iloc[i,0]*1.15, V.iloc[i,1]*1.15, V.index[i], color="b",
		# 	ha="center", va="center")
	plt.savefig("Scripts/Genomic_Prediction_RF/biplot.pdf")
	# plt.close()
	# plt.heatmap(V.iloc[:,:5])
	# plt.savefig("Scripts/Genomic_Prediction_RF/loadings_heatmap.pdf")
	return (X_r_sub, evar[:nPCs])


def Linear_mod(X, y, prefix, method="all"):
	""" Ordinary Least Squares Linear Regression """
	if method=="all":
		X_s = (X-X.mean())/X.std() # scale and center feature table
		X2 = sm.add_constant(X_s)
		mod = sm.OLS(np.array(y).reshape(-1,1),
			np.array(X2).reshape(-1,X2.shape[1]))
		res = mod.fit()
		yhats = res.predict()
		r2_scores = r2_score(y, yhats) # calculate variation explained by all the factors
		with open(f"{prefix}_ols_env_results.txt", "w") as out:
			out.write(res.summary().as_text())
			out.write(f"\nR-sq: {r2_scores}") # formula might be a bit different from sm.ols
		# vars(res) # attributes
		pickle.dump(mod, open(f"{prefix}_ols_env.pkl", 'wb')) # save the model
		coefs = pd.Series(res.params[1:])
		coefs = pd.concat([coefs, pd.Series(res.pvalues[1:]), pd.Series(X.columns)], axis=1)
		coefs.columns = ["coef", "p-value", "factor"]
		coefs.index.name = "Factor"
		coefs.sort_values("coef", ascending=False, inplace=True)
		coefs.to_csv(f"{prefix}_ols_env_coefs.csv")
		yhats=pd.Series(yhats)
		yhats.index = y.index
		yhats.name = 'y_pred'
		pd.concat([y, yhats], axis=1).to_csv(f"{prefix}_ols_env_preds.csv")
	if method=="byCol":
		coefs = pd.DataFrame()
		yhats = pd.DataFrame()
		r2_scores = []
		for col in X.columns:
			X2 = sm.add_constant(X[col])
			mod = sm.OLS(np.array(y).reshape(-1, y.shape[1]),
				np.array(X2).reshape(-1, X2.shape[1]))
			res = mod.fit()
			print(res.summary())
			coefs[col] = pd.Series(mod.coef_[0][0])
			# calculate variation explained by all the factors
			# yhat = np.dot(np.array(X_r).reshape(-1,X_r.shape[1]),res.params[1:]) + res.params[1]
			yhat = res.predict()
			yhats[col] = pd.Series(yhat.ravel())
			r2_scores.append(r2_score(y, yhat))
			print("R-sq:", r2_score(y, yhat))
		coefs.to_csv(f"{prefix}_ols_factors_coefs.csv")
		yhats.to_csv(f"{prefix}_ols_factors_preds.csv")
		r2_scores = pd.Series(r2_scores)
		r2_scores.index = X.index
		r2_scores.to_csv(f"{prefix}_ols_factors_results.csv")
	if method=="pc":
		# X_r, evar = pca(X.T)
		# factors_pc = copy.deepcopy(X_r)
		# factors_pc['sum'] = 0
		# for i in range(X_r.shape[0]):
		# 	factors_pc.iloc[i, factors_pc.shape[1]-1] = factors_pc.iloc[i,:].abs().sum() # row sums of loadings across eigenvectors
		# factors_pc.sort_values(by="sum", ascending=False)
		# factors_pc.to_csv(f"{prefix}_ols_env_pca_matrix.csv") # save sorted PCA matrix
		X_r, evar = pca(X) # re-run PCA, we need to have 35 rows, not 81
		X2 = sm.add_constant(X_r)
		X_r.to_csv(f"{prefix}_ols_env_pca_matrix.csv") # save sorted PCA matrix
		mod = sm.OLS(np.array(y).reshape(-1,1),
			np.array(X2).reshape(-1,X2.shape[1]))
		res = mod.fit()
		yhats = res.predict()
		r2_scores = r2_score(y, yhats) # calculate variation explained by all the factors
		with open(f"{prefix}_ols_env_pca_results.txt", "w") as out:
			out.write(res.summary().as_text())
			out.write(f"\nR-sq: {r2_scores}") # formula might be a bit different from sm.ols
		# vars(res) # attributes
		pickle.dump(mod, open(f"{prefix}_ols_env_pca.pkl", 'wb')) # save the model
		coefs = pd.Series(res.params[1:])
		coefs = pd.concat([coefs, pd.Series(res.pvalues[1:]), pd.Series(evar)], axis=1)
		coefs.columns = ["coef", "p-value", "explained_variance"]
		coefs.index.name = "PC"
		coefs.sort_values("coef", ascending=False, inplace=True)
		coefs.to_csv(f"{prefix}_ols_env_pca_coefs.csv")
		yhats=pd.Series(yhats)
		yhats.index = y.index
		yhats.name = 'y_pred'
		pd.concat([y, yhats], axis=1).to_csv(f"{prefix}_ols_env_pca_preds.csv")


def hyperparameter_tuning(reg, params, nGS, X_train, y_train, prefix):
	"""Hyperparameter tuning using grid search and cross-validation"""
	print("==========PERFORMING GRID SEARCH ON HYPERPARAMETER SPACE ==========")
	start_time = time.time()
	# Dataframe to save gridsearch results to    
	gs_results = pd.DataFrame(
		columns=["mean_test_score", "params"])
	for rep in range(nGS):
		print(f"Round {rep+1} of {nGS}")
		gs = GridSearchCV(estimator=reg,
			param_grid=params,
			scoring="r2",
			cv = 5,
			n_jobs=-1,
			verbose=2)
		# fitted_model = gs.fit(X_train,y_train) # fit the model
		# return fitted_model.best_params_, fitted_model.best_estimator_)
		gs.fit(X_train, y_train) # fit to training data
		rep_results = pd.DataFrame(gs.cv_results_)
		gs_results = pd.concat([gs_results, 
			rep_results[["params", "mean_test_score"]]])
	
	# Break params into seperate columns
	gs_results2 = pd.concat([gs_results.drop(["params"], axis=1),
			gs_results["params"].apply(pd.Series)], axis=1)
	
	# Find the mean scores for each parameter combination across nGS reps
	param_names = list(gs_results2)[1:] # parameters tested
	gs_results_mean = gs_results2.groupby(param_names).mean() # mean of mean test scores 
	gs_results_mean = gs_results_mean.sort_values("mean_test_score",
		0, ascending=False) # sort
	top_params = gs_results_mean.index[0] # best parameter combination
	print(gs_results_mean.head()) # print gridsearchcv results
	print(f"Parameter sweep time: {(time.time()- start_time)} seconds") # elapsed time
	
	# Save grid search results
	outName = open(f"{prefix}_GridSearch.csv", "w")
	outName.write(f"# {(time.time() - start_time)} sec\n")
	gs_results_mean.to_csv(outName)
	outName.close()
	return top_params

def DecisionTree_mod(X, y, n, nGS, prefix):
	"""Decision Tree Regression model training with cross-validation"""
	# stratified train-test split
	X_train, X_test, y_train, y_test = train_test_split(
		X, y, test_size=0.2, random_state=42)
	
	# tune hyperparameters
	reg = DecisionTreeRegressor(random_state=random.seed(123))
	params = {"max_depth":[1, 2, 3, 5],
		"max_features":[0.1, 0.5, "auto", "sqrt", "log2", None]}
	best_params = hyperparameter_tuning(reg, params, nGS, X_train, y_train, prefix)
	print("Best parameters: ", best_params)
	max_depth, max_features = best_params
	
	# train the model
	start_time = time.time()
	results_cv = []
	results_test = []
	imp = pd.DataFrame()
	for i in range(n):
		print(f"Running {i+1} of {n}")
		
		reg = DecisionTreeRegressor(max_depth=max_depth, 
			max_features=max_features, random_state=i)
		
		# train best estimator
		cv_pred = cross_val_predict(
			reg, X, y, cv=4, n_jobs=-1) # predictions
		
		# Performance statistics on validation set
		mse_val = mean_squared_error(y_train, cv_pred)
		rmse_val = np.sqrt(mean_squared_error(y_train, cv_pred))
		evs_val = explained_variance_score(y_train, cv_pred)
		r2_val = r2_score(y_train, cv_pred)
		cor_val = np.corrcoef(np.array(y_train), cv_pred)
		print("Val MSE: %f" % (mse_val))
		print("Val RMSE: %f" % (rmse_val))
		print("Val R-sq: %f" % (r2_val))
		print("Val PCC: %f" % (cor_val[0, 1]))
		result_val = [mse_val, rmse_val, evs_val, r2_val, cor_val[0, 1]]
		results_cv.append(result_val)

		# Evaluate the model on the test set
		reg.fit(X_train, y_train)
		y_pred = reg.predict(X_test)
		imp[i] = reg.feature_importances_
		# yhats[i] = y_pred # save predictions

		# Performance statistics on the test set
		mse = mean_squared_error(y_test, y_pred)
		rmse = np.sqrt(mean_squared_error(y_test, y_pred))
		evs = explained_variance_score(y_test, y_pred)
		r2 = r2_score(y_test, y_pred)
		cor = np.corrcoef(np.array(y_test), y_pred)
		print("Test MSE: %f" % (mse))
		print("Test RMSE: %f" % (rmse))
		print("Test R-sq: %f" % (r2))
		print("Test PCC: %f" % (cor[0, 1]))
		result_test = [mse, rmse, evs, r2, cor[0, 1]]
		results_test.append(result_test)
	
	print(f"Training time: {(time.time() - start_time)} seconds")
	return (results_cv, results_test, imp)



if __name__ == "__main__":
	# Argument parser
	parser = argparse.ArgumentParser(description="")
	req_group = parser.add_argument_group(title="Required Input")
	req_group.add_argument(
		"-Y", type=str, help="Path to CSV file of model performances")
	req_group.add_argument(
		"-alg", type=str, help="Algorithm to run (lin or dt)")
	req_group.add_argument(
		"-prefix", type=str, help="Prefix of file name to save models as. Prefix may include file path too.")
	opt_group = parser.add_argument_group(title="Optional Input")
	opt_group.add_argument(
		"-sep", type=str, help="Delimeter for file Y", default=",")
	opt_group.add_argument(
		"-index_col", type=str, help="Name of column in Y to set as index", default="cond")
	opt_group.add_argument(
		"-y_col", type=str, help="Name of column in Y to set as dependent variable", default="r2_test")
	opt_group.add_argument(
		"-method", type=str,
		help="Method for linear model ('all' features, 'byCol' one feature at a time, 'pc' run on PCA matrix of features",
		default="all")
	opt_group.add_argument(
		"-nGS", type=int, help="Number of times to repeat grid search for DecisionTree", default=10)
	opt_group.add_argument(
		"-n", type=int, help="Number of times to repeat DecisionTree model training", default=10)
	

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
	y = pd.read_csv(args.Y, sep=args.sep)  # Read in model performance data

	h2.set_index("Conditions", inplace=True) # set environments as row index
	y.set_index(args.index_col, inplace=True) # set environments as row index
	y = y.loc[:,args.y_col] # keep only dependent variable column
	X = process_data(df)  # make the feature table
	X = pd.concat([h2.h2, X], axis=1) # add heritabilities
	X.to_csv("Scripts/Genomic_Prediction_RF/factors_features.csv")# 11 factors
	# X.to_csv("Scripts/Genomic_Prediction_RF/factors_features_for_pca.csv") # 81 factors

	# Build models
	if args.alg == "lin":
		Linear_mod(X, y, args.prefix, args.method)  # fit the linear model

	if args.alg == "dt":
		dtmod, imps = DecisionTree_mod(X, y, args.n, args.nGS, args.prefix) # fit the decision tree model (this model performs really poorly)

	# determine which variables are important
	# X_r, evar = pca(X)
	# factors_pc = copy.deepcopy(X_r)
	# factors_pc['sum'] = 0
	# for i in range(81):
	# 	factors_pc.iloc[i,26] = factors_pc.iloc[i,0:4].abs().sum() # sum of loadings in the first 4 eigenvectors
	# factors_pc['sum'].sort_values(ascending=False) # features that captures the most variance based on PCs
	# factors_pc['cor'] = 0
	# Xcorr = X.corr().stack().reset_index() # did not scale features
	# left = pd.DataFrame(factors_pc['sum'])
	# right = pd.DataFrame(Xcorr)
	# df = pd.merge(left, right, left_index=True, right_on='level_0')
	# df = df.loc[df[0]<1] # remove correlations of 1
	# df2 = df.loc[df[0]>.8] # highly correlated features
	# df2.level_0.unique()
	# df2.loc[df2.level_0=='75%'] # see which features correlate with it
	# df2.loc[df2.level_0=='cor_YPDCHX05']
	# df2.loc[df2.level_0=='cor_YPDNACL15M']
	# Xcorr.to_csv('Scripts/Genomic_Prediction_RF/corr_factors.csv') # save correlation matrix
	# X_s = (X-X.mean())/X.std() # scale and center feature table
	# xy = pd.concat([y,X_s], axis=1)
	# xy.corr().r2_test.sort_values() # see which features correlate with r2_test
	'''In the end I just stuck with the 29 PCs. I was worried that I needed 
	to do more stringent feature selection because non of the coefficients
	had significant p-values. Shinhan thinks it's ok, bc we want the estimated
	coefs to be close to the true coefs (accept null hypothesis). Thus, the 
	higher the p-value, the closer the estimates are to the true values.'''
