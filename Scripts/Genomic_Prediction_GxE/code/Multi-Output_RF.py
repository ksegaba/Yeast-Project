#!/usr/bin/env python3
"""
Multi-Output Random Forest Regressor for GxE analysis

Arguments:
    -X (str): 
    -y (str):
    -label1 (str): 
    -label2 (str):
    -test (str):
    -save (str):
    -feat (str):
    -Xindex (str):
    -nGS (int):
    -nTrain (int):

Output:
    [1] _GridSearch.csv
    [2] _fitted_model.pkl
    [3] _results_cv.csv
    [4] _results_test.csv
    [5] _imp.csv
    [6] _preds_cv.csv
    [7] _preds_test.csv
"""
__author__ = "Kenia Segura Abá"

import warnings
import sys
import os
import argparse
import datetime
import time
import pickle
import datatable as dt
import pandas as pd
import numpy as np
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import explained_variance_score
from numpy import corrcoef


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


def GridSearch(hyperparameters, nGS, X_train, y_train, save):
    """
    Perform a repeated hyperparameter sweep using GridSearchCV with 5 folds.
    """
    start_time = time.time()
    
    # Model
    model = MultiOutputRegressor(RandomForestRegressor(random_state=123))
    
    # Dataframe to save gridsearch results to    
    gs_results = pd.DataFrame(
        columns=["mean_test_score", "params"])
    for rep in range(nGS):
        print(f"Round {rep+1} of {nGS}")
        gs = GridSearchCV(
            estimator=model, 
            param_grid=hyperparameters, 
            cv=5,
            scoring="r2",
            verbose=2,
            n_jobs=-1)
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
    outName = open(f"{save}_GridSearch.csv", "w")
    outName.write(f"# {(time.time() - start_time)} sec\n")
    gs_results_mean.to_csv(outName)
    outName.close()
    return top_params


def Run_MultiOutputRF_reg(X, y, test, nGS, nTrain, save):
    """
    Perform a Multi-Output Random Forest Regression with 
    hyperparameter tuning and 5-fold cross-validation to train the model.
    """
    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test)

    # Hyperparameter tuning
    hyperparameters = dict(
        estimator__max_depth=[3, 5, 10],
        estimator__max_features=[0.1, 0.5, "sqrt", "log2"],
        estimator__n_estimators=[100, 300, 500, 700, 1000])
    top_params = GridSearch(hyperparameters, nGS, X_train, y_train, save)

    print("==========PERFORMING 5-FOLD CV TO TRAIN THE BEST MODEL ==========")
    start_time = time.time()
    results_cv = {} # cross-validation results
    results_test = {} # test set results
    imp = {} # feature importance scores
    cv_preds = {}
    test_preds = {}
    
    # Train the model using 5-fold cross-validation
    max_depth, max_features, n_estimators = top_params # hyperparameters
    
    for rep in range(nTrain):
        print(f"Running {rep+1} of {nTrain}")
        rep_name = f"rep_{str(rep + 1)}"
        imp[rep] = {}
        
        reg = MultiOutputRegressor(RandomForestRegressor(
            max_depth=int(max_depth), 
            max_features=max_features, 
            n_estimators=int(n_estimators),
            random_state=rep)) # model to train
        
        cv_pred = cross_val_predict(
            estimator=reg, X=X_train, y=y_train, cv=5, n_jobs=5, verbose=2)
        
        # Get cross validation performance statistics and predicted values
        cv_pred_df = pd.DataFrame(data=cv_pred, index=X_train.index, columns=y_train.columns)
        cv_preds[rep] = cv_pred_df.T.to_dict()
        
        r2 = r2_score(y_train, cv_pred, multioutput="raw_values")
        mse = mean_squared_error(y_train, cv_pred, multioutput="raw_values")
        evs = explained_variance_score(y_train, cv_pred, multioutput="raw_values")
        cor = []
        for trait in range(cv_pred.shape[1]):
            cor.append(corrcoef(np.array(y_train.iloc[:,trait]),
                np.array(cv_pred_df.iloc[:,trait]))[0, 1])
        
        result = {}
        for trait in range(y_train.shape[1]):
            result[y_train.columns[trait]] = np.array(
                [r2[trait], mse[trait], evs[trait], cor[trait]])
        results_cv[rep] = result

        # Fit the model and evaluate it
        reg.fit(X_train, y_train) # fit the model
        test_pred = reg.predict(X_test) # evaluate the model
        for i in range(len(reg.estimators_)):
            imp[rep][y_test.columns[i]] = reg.estimators_[i].feature_importances_

        # Save the fitted model 
        with open(f"{save}_fitted_model.pkl", "wb") as f:
            pickle.dump(reg, f)

        # Get performance statistics on the test set and predicted values
        test_pred_df = pd.DataFrame(
            data=test_pred, index=X_test.index, columns=y_test.columns)
        test_preds[rep] = test_pred_df.T.to_dict()
        
        r2_test = r2_score(y_test, test_pred, multioutput="raw_values")
        mse_test = mean_squared_error(y_test, test_pred, multioutput="raw_values")
        evs_test = explained_variance_score(y_test, test_pred, multioutput="raw_values")
        cor_test = []
        for trait in range(test_pred.shape[1]):
            cor_test.append(corrcoef(np.array(y_test.iloc[:,trait]),
                np.array(test_pred_df.iloc[:,trait]))[0, 1])
        
        result_test = {}
        for trait in range(y_test.shape[1]):
            result_test[y_test.columns[trait]] = np.array(
                [r2_test[trait], mse_test[trait], evs_test[trait], cor_test[trait]])
        results_test[rep] = result_test

    print(f"Training time: {(time.time() - start_time)} seconds")

    return (results_cv, results_test, imp, cv_preds, test_preds)


def nestedDict_to_df(nested_dict):
    new_dict = {}
    for outerKey, innerDict in nested_dict.items():
        for innerKey, values in innerDict.items():
            new_dict[(outerKey, innerKey)] = values
    return (pd.DataFrame(new_dict).T)


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="")
    # Required input
    req_group = parser.add_argument_group(title="Required Input")
    req_group.add_argument(
        "-X", type=str, help="path to feature data file", required=True)
    req_group.add_argument(
        "-y", type=str, help="path to label data file", required=True)
    req_group.add_argument(
        "-test", type=str, help="path to file of test set instances", required=True)
    req_group.add_argument(
        "-save", type=str, help="prefix for output files (may include path)", required=True)
    # Optional input
    opt_group = parser.add_argument_group(title="Optional Input")
    opt_group.add_argument(
        "-feat", type=str, help="", default="all")
    opt_group.add_argument(
        "-Xindex", type=str, help="name of index column in X", default="ID")
    opt_group.add_argument(
        "-labels", type=list, nargs='+', help="Comma delimited column names of label data file", default="all")
    opt_group.add_argument(
        "-nGS", type=int, help="Number of times to repeat grid search", default=10)
    opt_group.add_argument(
        "-nTrain", type=int, help="Number of times to repeat model training", default=100)

    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args() # Read arguments

    # Debugging arguments
    # X = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv")
    # X = X.to_pandas()
    # X.set_index('ID', drop=True, inplace=True)
    # y = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", index_col=0)
    # test = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", header=None)
    # feat = pd.read_csv("/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs/feat_rf_YPDCUSO410MM_top_1000", header=None)
    # labels = ["YPDCUSO410MM", "YPDBENOMYL500"]
    # save = "multi-output-RF-YPDCUSO410MM"

    # Read in data
    X = dt.fread(args.X) # features
    X = X.to_pandas()
    X.set_index(args.Xindex, drop=True, inplace=True)
    y = pd.read_csv(args.y, index_col=0) # labels
    test = pd.read_csv(args.test) # test instances

    # Subset data
    if args.feat != "all":
        feat = pd.read_csv(args.feat)
        X = X.loc[:,feat[0]]
    
    if args.labels != "all":
        print(args.labels)
        labels = []
        for i in range(len(args.labels)):
            strng=""
            for char in args.labels[i]:
                if char != ",":
                    strng += char
            labels.append(strng)
        y = y[labels]

    # Hyperparameter tuning, model training, and evaluation
    results_cv, results_test, imp, cv_preds, test_preds = Run_MultiOutputRF_reg(
        X, y, test, args.nGS, args.nTrain, args.save)

    # Save results to files
    results_cv_df = nestedDict_to_df(results_cv)
    results_cv_df.columns = ["r2", "mse", "evs", "pcc"]
    results_cv_df.index.rename(['rep', 'trait'], inplace=True)
    results_cv_df.to_csv(f"{args.save}_results_cv.csv")
    results_test_df = nestedDict_to_df(results_test)
    results_test_df.columns = ["r2", "mse", "evs", "pcc"]
    results_test_df.index.rename(['rep', 'trait'], inplace=True)
    results_test_df.to_csv(f"{args.save}_results_test.csv")
    imp_df = nestedDict_to_df(imp)
    imp_df.columns = X.columns
    imp_df.index.rename(['rep', 'trait'], inplace=True)
    imp_df.to_csv(f"{args.save}_imp.csv")
    cv_preds_df = nestedDict_to_df(cv_preds)
    cv_preds_df.index.rename(['rep', 'instance'], inplace=True)
    cv_preds_df.to_csv(f"{args.save}_preds_cv.csv")
    test_preds_df = nestedDict_to_df(test_preds)
    test_preds_df.index.rename(['rep', 'instance'], inplace=True)
    test_preds_df.to_csv(f"{args.save}_preds_test.csv")
