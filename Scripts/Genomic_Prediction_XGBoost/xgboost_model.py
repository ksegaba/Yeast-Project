"""
XGBoost for regression on SNP, ORF, and CNO data from Peter et al. 2018

# Required Inputs
    -X      Path to feature matrix
    -Y      Path to label matrix
    -test   Path to list of test instances file
    -trait  Column name of target trait in Y matrix
    -save   Path to folder to save output files to
    # Optional
    -feat   Path to list of training features file (default is all)
    -type   Feature types (e.g. SNP (default), ORF, CNO)
    -fold   k folds for Cross-Validation (default is 5)
    -n      Number of CV repetitions (default is 10)

# Outputs (prefixed with <type>_<trait>)
    _GridSearch.csv     Grid Search CV R-sq scores
    _lm_test.pdf        Regression plot of predicted and actual test labels
    _model.save         XGBoost model
    _imp.csv            Feature importance scores
    _top20.pdf          Plot of top 20 features' importance scores
    _cv_results.csv     Cross-validation results (various metrics)
    _test_results.csv   Evaluation results (various metrics)
    RESULTS_xgboost.txt Aggregated results (various metrics)    
"""
__author__ = "Kenia Segura Abá"

from configparser import ExtendedInterpolation
import sys, os, argparse, time
from datetime import datetime
import random
from operator import index
import datatable as dt
import pandas as pd
import numpy as np
import pickle
import xgboost as xgb
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.metrics import explained_variance_score
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import seaborn as sns

# Need to add args
# Need trait arg bc the for loop may cause it to crash on the command line
# So I need a slurm script to keep track 

def tune_model(xg_reg, X_train, y_train, data_type, trait, fold):
    """ Hyperparameter sweep using grid search with cross-validation """
    start = time.time()

    # Hyperparameter tuning
    parameters = {"gamma": [0,0.01], # min split loss
                "eta":[0.3, 0.1, 0.03], # learning rate
                "max_depth":[3, 6, 7], # tree depth
                "subsample": [0.3, 1], # instances per tree
                "colsample_bytree": [0.3, 0.5, 1], # features per tree
                "n_estimators": [100, 300, 500]} # sample training instances
    gs = GridSearchCV(estimator=xg_reg,
                param_grid=parameters,
                cv=fold,           # cross validation folds
                scoring='r2',   # find model with the best R-squared
                verbose=2,
                n_jobs=10)
    fitted_model = gs.fit(X_train, y_train) # fit the model

    run_time = time.time() - start
    print("GridSearch Run Time: %f" % (run_time))

    # Write results to file
    gs_results = pd.DataFrame(gs.cv_results_)
    gs_results.to_csv(f"{args.save}/{data_type}_{trait}_GridSearch.csv")

    print("Best parameters: ", fitted_model.best_params_)

    return (fitted_model.best_params_, gs.best_estimator_)


def xgb_reg(trait, fold, n, data_type):
    """ Train XGBoost Regression Model """
    print(trait)
    y = Y[trait]

    # Train-test split
    X_train = X.loc[~X.index.isin(test[0])]
    X_test = X.loc[X.index.isin(test[0])]
    y_train = y.loc[~y.index.isin(test[0])]
    y_test = y.loc[y.index.isin(test[0])]

    # Convert trainig set to Dmatrix data structure 
    data_dmatrix = xgb.DMatrix(data=X_train,label=y_train)

    # Instantiate XGBoost regressor object
    xg_reg = xgb.XGBRegressor(random_state=random.seed(123))

    # Hyperparameter tuning
    best_params, best_model = tune_model(
        xg_reg, X_train, y_train, data_type, trait, fold)

    ################## Training with Cross-Validation ##################
    results_cv = [] # hold performance metrics of cv reps
    results_test = [] # hold performance metrics on test set

    # Training with Cross-validation
    # cv = RepeatedKFold(n_splits=fold, n_repeats=n, random_state=123)
    # cv_scores = xgb.cv(
    #     params=best_params,
    #     dtrain=data_dmatrix,
    #     num_boost_round=10,
    #     stratified=True,
    #     folds=cv,
    #     metrics="rmse", # there is no r2_score equivalent
    #     seed=123)

    for j in range(0, n): # repeat cv 10 times
        print(f"Running {j+1} of {n}")
        
        cv = StratifiedKFold(n_splits=fold) # stratified samples
        cv_pred = cross_val_predict(
            best_model, X_train, y_train, cv=fold, n_jobs=-1) # predictions

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
    #best_model.fit(X_train, y_train)
    y_pred = best_model.predict(X_test)

    # Performance on the test set
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

    # Plot linear regression of actual and predicted test values
    sns.regplot(x=y_pred, y=y_test, fit_reg=True, ci=95, seed=123, color="black")
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title(f"{trait}")
    plt.savefig(f"{args.save}/{data_type}_{trait}_lm_test.pdf", format="pdf")

    # Save the fitted model to a file
    filename = f"{args.save}/{data_type}_{trait}_model.save"
    pickle.dump(best_model, open(filename, 'wb'))

    # Save feature importance scores to file
    importances = pd.DataFrame(best_model.feature_importances_)
    importances.to_csv(f"{args.save}/{data_type}_{trait}_imp.csv")

    # Plot feature importances
    xgb.plot_importance(
        best_model, grid="off", max_num_features=20, 
        title=f"{trait} Feature Importances", xlabel="Weight")
    plt.savefig(f"{args.save}/{data_type}_{trait}_top20.pdf", format="pdf")
    
    return (results_cv, results_test)


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(
        description="XGBoost Regression on SNP and ORF data")
    
    # Required input
    req_group = parser.add_argument_group(title="Required Input")
    req_group.add_argument(
        "-X", help="path to feature data file", required=True)
    req_group.add_argument(
        "-Y", help="path to label data file", required=True)
    req_group.add_argument(
        "-test", help="path to file of test set instances", required=True)
    req_group.add_argument(
        "-trait", help="name of trait column in y dataframe", required=True)
    req_group.add_argument(
        "-save", help="path to save output files", required=True)
    
    # Optional input
    req_group.add_argument(
        "-feat", help="Path to list of training features file (default is all)", default="all")
    req_group.add_argument(
        "-type", help="data type of X matrix (e.g. SNP, ORF, CNO)", default="SNP")
    req_group.add_argument(
        "-fold", help="k number of cross-validation folds", default=5)
    req_group.add_argument(
        "-n", help="number of cross-validation repetitions", default=10)
    
    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args() # Read arguments

    # Read in data
    X = dt.fread(args.X)
    X = X.to_pandas()
    X.set_index(X.columns[0], inplace=True)
    Y = pd.read_csv(args.Y, index_col=0)
    test = pd.read_csv(args.test, sep="\t", header=None)

    # Filter out features not in feat file
    if args.feat != "all":
        print("Using subset of %s features" % args.feat)
        with open(args.feat) as f:
            features = f.read().strip().splitlines()
            features = ["ID"] + features
        X = X.reindex(columns=features)
        print(X.shape)

    # Train the model
    start = time.time()
    results_cv, results_test = xgb_reg(
        args.trait, int(args.fold), int(args.n), args.type)
    run_time = time.time() - start
    print("Training Run Time: %f" % (run_time))

    # Save results to file
    results_cv = pd.DataFrame(
        results_cv, 
        columns=["MSE_val", "RMSE_val", "EVS_val", "R2_val", "PCC_val"])
    results_cv.to_csv(f"{args.save}/{args.type}_{args.trait}_cv_results.csv")
    results_test = pd.DataFrame(
        results_test, 
        columns=["MSE_test", "RMSE_test", "EVS_test", "R2_test", "PCC_test"])
    results_test.to_csv(
        f"{args.save}/{args.type}_{args.trait}_test_results.csv")

    # Aggregate results and save to file
    if not os.path.isfile(f"{args.save}/RESULTS_xgboost.txt"):
        out = open(f"{args.save}/RESULTS_xgboost.txt", "a")
        out.write("Date\tRunTime\tData\tTrait\tMSE_val\tMSE_val_sd\
            \tRMSE_val\tRMSE_val_sd\tEVS_val\tEVS_val_sd\tR2_val\
            \tR2_val_sd\tPCC_val\tPCC_val_sd\tMSE_test\tMSE_test_sd\
            \tRMSE_test\tRMSE_test_sd\tEVS_test\tEVS_test_sd\tR2_test\
            \tR2_test_sd\tPCC_test\tPCC_test_sd\n")
        out.close()

    out = open(f"{args.save}/RESULTS_xgboost.txt", "a")
    out.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\t\
        {run_time}\t{args.type}\t{args.trait}\t\
        {np.mean(results_cv.MSE_val)}\t{np.std(results_cv.MSE_val)}\t\
        {np.mean(results_cv.RMSE_val)}\t{np.std(results_cv.RMSE_val)}\t\
        {np.mean(results_cv.EVS_val)}\t{np.std(results_cv.EVS_val)}\t\
        {np.mean(results_cv.R2_val)}\t{np.std(results_cv.R2_val)}\t\
        {np.mean(results_cv.PCC_val)}\t{np.std(results_cv.PCC_val)}\t\
        {np.mean(results_test.MSE_test)}\t{np.std(results_test.MSE_test)}\t\
        {np.mean(results_test.RMSE_test)}\t{np.std(results_test.RMSE_test)}\t\
        {np.mean(results_test.EVS_test)}\t{np.std(results_test.EVS_test)}\t\
        {np.mean(results_test.R2_test)}\t{np.std(results_test.R2_test)}\t\
        {np.mean(results_test.PCC_test)}\t{np.std(results_test.PCC_test)}")
