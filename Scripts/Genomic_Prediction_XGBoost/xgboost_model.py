"""
XGBoost for regression

Kenia Segura Aba
"""
import sys, os, time
from datetime import datetime
import random
from operator import index
import datatable as dt
import pandas as pd
import numpy as np
import pickle
import xgboost as xgb
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import GridSearchCV
random.seed(123)

# Need to add args
# Need trait arg bc the for loop may cause it to crash on the command line
# So I need a slurm script to keep track 

def hyperparameter_tuning(X_train, y_train):
    # Parameters to tune
    param_tuning = {
        'learning_rate': [0.01, 0.1],
        'max_depth': [3, 5, 7, 10],
        'min_child_weight': [1, 3, 5],
        'subsample': [0.5, 0.7],
        'colsample_bytree': [0.5, 0.7],
        'n_estimators' : [100, 200, 500],
        'objective': ['reg:squarederror']
    }
    xgb_model = xgb.XGBRegressor() # model
    # Grid Search
    gsearch = GridSearchCV(estimator = xgb_model,
                           param_grid = param_tuning,
                           scoring = 'neg_mean_squared_error',
                           cv = 5,
                           n_jobs = -1,
                           verbose = 1)
    gsearch.fit(X_train,y_train) # fit the model
    return gsearch.best_params_

def main():
    # Read in data
    #X = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv")
    #X = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv") # ORF copy number
    X = dt.fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv") # ORF presence/absence
    X = X.to_pandas()
    X.set_index(X.columns[0], inplace=True)
    Y = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", index_col=0)
    test = pd.read_csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", sep="\t", header=None)

    data_type = "ORF" # SNP  CNO
    #tmp=["YPDSDS", "YPRIBOSE"]#"YPDFLUCONAZOLE", "YPDMV"] #"YPDANISO10", "YPDETOH"]#, , 
    tmp=["YPDFORMAMIDE5"]
    for trait in tmp:#Y.columns:
        print(trait)
        y = Y[trait]

        # Convert to Dmatrix data structure 
        data_dmatrix = xgb.DMatrix(data=X,label=y)

        # Create train and test set
        X_train = X.loc[~X.index.isin(test[0])]
        X_test = X.loc[X.index.isin(test[0])]
        y_train = y.loc[~y.index.isin(test[0])]
        y_test = y.loc[y.index.isin(test[0])]

        # Parameter sweep using grid search
        start = time.time()
        #params = hyperparameter_tuning(X_train, y_train)

        # Instantiate XGBoost regressor object with the best parameters
        xg_reg = xgb.XGBRegressor(objective='reg:squarederror', colsample_bytree=0.3, learning_rate=0.1,
                                max_depth=5, alpha=10, n_estimators=500, random_state=random.seed(123))
        
        # Hyperparameter tuning
        """xg_reg = xgb.XGBRegressor(random_state=random.seed(123))
        parameters = {'nthread':[4], 
                    'objective':['reg:squarederror'],
                    'learning_rate': [.03, 0.05, .07, .1], #so called `eta` value
                    'max_depth': [3, 5, 7, 9],
                    'min_child_weight': [4],
                    'silent': [1],
                    'subsample': [0.7],
                    'colsample_bytree': [0.3, 0.5, 0.7],
                    'alpha': [5, 7, 10],
                    'n_estimators': [200, 500]}
        xgb_gs = GridSearchCV(
                        xg_reg,
                        parameters,
                        cv=5,           # cross validation folds
                        verbose=1, 
                        scoring='r2',   # find model with the best R-squared
                        n_jobs=8)       # number of concurrent jobs, you need to
        
        # Fit the regressor to the training set and evaluate on the test set
        fitted_model = xgb_gs.fit(X_train,y_train)
        best_model = xgb_gs.best_estimator_ # best model
        y_pred = best_model.predict(X_test)

        # Save the model
        filename = "%s_model_xgb_gridsearch.save"%trait
        pickle.dump(xgb_gs.best_estimator_, open(filename, 'wb'))
        """
        # Cross-validation (this is done during gridsearch)
        preds_val = cross_val_predict(xg_reg, X_train, y_train, cv=10, n_jobs=-1)
        #scores = cross_val_score(xg_reg, X_train, y_train, cv=5, n_jobs=-1)
        
        # Cross-validation method #2
        #params = {"objective":"reg:squarederror",'colsample_bytree': 0.3,'learning_rate': 0.1, 'max_depth': 5, 'alpha': 10}
        #cv_results = xgb.cv(dtrain=data_dmatrix, params=params, nfold=5,
        #            num_boost_round=500,early_stopping_rounds=10,metrics="r2",as_pandas=True,seed=123)

        # Fit the regressor to the training set and evaluate on the test set
        xg_reg.fit(X_train,y_train)
        y_pred = xg_reg.predict(X_test)

        # Save the model
        #filename = "%s_model_xgb_not_tuned.save"%trait
        #filename = "%s_cno_model_xgb_not_tuned.save"%trait
        filename = "%s_orf_model_xgb_not_tuned.save"%trait
        pickle.dump(xg_reg, open(filename, 'wb'))

        # Model performance
        rmse_val = np.sqrt(mean_squared_error(y_train, preds_val))
        r2_val = r2_score(y_train, preds_val)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        r2 = r2_score(y_test, y_pred)
        print("RMSE: %f" % (rmse))
        print("R-sq: %f" % (r2))
        run_time = time.time() - start
        print("Run Time: %f" % (run_time))
            
        # Save results to file
        if not os.path.isfile("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_XGBoost/RESULTS_xgboost.txt"):
            out = open("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_XGBoost/RESULTS_xgboost.txt", "a")
            out.write("Date\tRunTime\tData\tTrait\tn\tRMSE_val\tR2_val\tRMSE_test\tR2_test\n")
            out.close()

        out = open("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_XGBoost/RESULTS_xgboost.txt", "a")
        out.write("%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\n"%(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),run_time,data_type,trait,rmse_val,r2_val,rmse,r2))

if __name__ == "__main__":
    main()