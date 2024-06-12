"""
Generate feature selection subsets for a given XGBoost model
"""

import sys
import argparse
import pickle
import pandas as pd
import datatable as dt
import numpy as np


def parse_arguments():
	"""Argument parser"""
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-m", "--model", type=str, help="path to XGBoost model rep_0 file (to get feature names)",
		required=True, default=""
	)
	parser.add_argument(
		"-i", "--imp", type=str, help="path to XGBoost feature importance file",
		required=False, default=""
	)
	parser.add_argument(
		"-start", type=int, help="Number of features to start with",
		required=True
	)
	parser.add_argument(
		"-stop", type=int, help="Number of features to end with",
		required=True
	)
	parser.add_argument(
		"-o", "--save", type=str, help="path and/or prefix to save feature subset as",
		required=True
	)
	parser.add_argument(
		"-step", type=int, help="Step size of features to select", default=None
	)
	parser.add_argument(
		"-base", type=int, help="Base of exponential number of features to select, e.g., base^start to base^stop",
		default=None
	)
	return parser.parse_args()

if __name__ == "__main__":
	# Arguments help
	if len(sys.argv)==1:
		print("python feature_selection.py [-m] [-start] [-stop] [-step] [-o] [-x optional] [-i optional]")
		sys.exit()
	
	# Arguments
	args = parse_arguments()
	model_file = args.model
	imp_file = args.imp
	start = args.start
	stop = args.stop
	step = args.step
	base = args.base
	save = args.save

	# Load model
	mod = pickle.load(open(model_file, "rb"))

	if imp_file != "":
		# Read in feature importances
		imp = dt.fread(imp_file).to_pandas()
		imp.index = mod.feature_names_in_ # set feature names as index
		imp = imp.iloc[:,1:]
		mean_imp = imp.mean(axis=1) # take the mean across importances per feature
		mean_imp.sort_values(ascending=False, inplace=True)
		og_featnum = len(mean_imp)
		mean_imp = mean_imp[mean_imp.values != 0.0] # drop non-important features
		imp = pd.DataFrame(mean_imp.copy(deep=True))
	else:
		# Read in feature importances
		imp = mod.feature_importances_
		imp = pd.DataFrame(imp, index=mod.feature_names_in_)
		imp.sort_values(by=0, ascending=False, inplace=True) # sort features
		og_featnum = len(imp)
		imp = imp.loc[imp.iloc[:,0]!=0.0,:] # drop non-important features

	
	# Select subsets of features
	if stop > len(imp):
		print(f"After removing unimportant features, only {len(imp)} out of {og_featnum} remained.")
		print(f"Setting stop argument to {len(imp)}")
		stop = len(imp)

	if step != None:
		for i in range(start, stop+step, step):
			if i == start:
				continue
			else:
				subset = imp.iloc[start:i,:].index.values
				print(f"saving {start} to {i} features to file {save}_top_{len(subset)}")
				np.savetxt(f"{save}_top_{len(subset)}", subset, fmt="%s")
	
	if base != None:
		for i in range(start, stop+1):
			if i == start-1:
				continue
			else:
				subset = imp.iloc[start-1:base**i,:].index.values
				print(f"saving {start-1} to {base**i} features to file {save}_top_{len(subset)}")
				np.savetxt(f"{save}_top_{len(subset)}", subset, fmt="%s")