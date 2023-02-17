#!/usr/bin/env python3
"""
Multi-Output Random Forest Regressor for GxE analysis

    "-X", help="path to feature data file", required=True)
    req_group.add_argument(
        "-Xindex", help="name of index column in X", required=True, default="ID")
    req_group.add_argument(
        "-y", help="path to label data file", required=True)
    req_group.add_argument(
        "-test", help="path to file of test set instances", required=True)
    req_group.add_argument(
        "-save", help="prefix for output files", required=True)

Output:

Command:
X=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv
y=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv
test=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt
python Multi-Output_RF.py -X $X -Xindex ID -y $y -test $test -save RF_GxE
"""
__author__ = "Kenia Segura Abá"

import warnings
import sys
import os
import argparse
import datetime
import datatable as dt
import pandas as pd
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score
from numpy import corrcoef
from numpy import std
import pickle

def warn(*args, **kwargs):
    pass


warnings.warn = warn


def train_test_split(X, y, test):
    """Split X and y into training and testing sets"""
    X_train = X.loc[~X.index.isin(test.iloc[:,0]),:]
    X_test = X.loc[X.index.isin(test.iloc[:,0]),:]
    y_train = y.loc[~y.index.isin(test.iloc[:,0]),:]
    y_test = y.loc[y.index.isin(test.iloc[:,0]),:]
    return (X_train, X_test, y_train, y_test)


def Multi_Output_RF_reg(X_train, y_train):
    """
    Multi-Output Random Forest Regressor with 
    hyperparameter tuning and cross-validation
    """
    # Model
    model = MultiOutputRegressor(RandomForestRegressor())

    # Hyperparameter tuning
    hyperparameters = dict(
        estimator__n_estimators=[50, 100, 300, 500, 700],
        estimator__max_depth=[3, 5, 10],
        estimator__min_samples_split=[2, 4, 7],
        estimator__min_samples_leaf=[1, 2, 5, 10],
        estimator__max_features=[1, "auto", "sqrt"],
        estimator__max_leaf_nodes=[5, 20, 50, 100, 300],
        estimator__min_impurity_decrease=[0, 0.1, 0.2])

    gs = GridSearchCV(
        estimator=model, 
        param_grid=hyperparameters, 
        cv=5,
        scoring="r2",
        verbose=2,
        n_jobs=-1)

    tuning = gs.fit(X_train, y_train) # fit to training data
    tuned_model = tuning.best_estimator_ # best model
    best_scores = tuning.best_score_ # best model scores

    # Cross-validation
    # for j in range(0, 10): # repeat cv 10 times
    #     print("Running %i of 10" % (j + 1))
    #     cv = StratifiedKFold(n_splits=5)
    #     cv_preds = cross_val_predict(
    #         model, X_train, y_train, cv=cv, n_jobs=-1)

    regr = MultiOutputRegressor(RandomForestRegressor(random_state=123)).fit(X_train, y_train)
    regr.predict(X[[0]])

    return(tuned_model, best_scores)#, cv_preds)

if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="")
    # Required input
    req_group = parser.add_argument_group(title="Required Input")
    req_group.add_argument(
        "-X", help="path to feature data file", required=True)
    req_group.add_argument(
        "-Xindex", help="name of index column in X", required=True, default="ID")
    req_group.add_argument(
        "-y", help="path to label data file", required=True)
    req_group.add_argument(
        "-test", help="path to file of test set instances", required=True)
    req_group.add_argument(
        "-save", help="prefix for output files", required=True)
    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args() # Read arguments

    # Read in data
    X = dt.fread(args.X) # features
    X = X.to_pandas()
    X.set_index(args.Xindex, drop=True, inplace=True)
    
    y = pd.read_csv(args.y, sep=args.sep) # labels
    
    test = pd.read_csv(args.test) # test instances

    # Debugging arguments
    # X = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv")
    # X = X.to_pandas()
    # X.set_index('ID', drop=True, inplace=True)
    # y = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", index_col=0)
    # test = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", header=None)

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test)

    # Training
    #model, cv_preds = Multi_Output_RF_reg
    model, cv_scores = Multi_Output_RF_reg

    # Evaluation
    y_pred = model.predict(X_test)

    # Results

    # Save model to file
    with open(f"{args.save}_model.pkl", "wb") as f:
        pickle.dump(model, f)