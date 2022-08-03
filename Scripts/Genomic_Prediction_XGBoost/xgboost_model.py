"""
XGBoost for regression

Kenia Segura Aba
"""
from configparser import ExtendedInterpolation
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
from sklearn.metrics import explained_variance_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import GridSearchCV
random.seed(123)

# Need to add args
# Need trait arg bc the for loop may cause it to crash on the command line
# So I need a slurm script to keep track 

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
    #tmp=["YPDSDS", "YPRIBOSE"]#"YPDFLUCONAZOLE", "YPDMV"] #"YPDANISO10", "YPDETOH"]#, , "YPDFORMAMIDE5"
    tmp=['YPDBENOMYL500']
    for trait in tmp:#Y.columns: #
        print(trait)
        y = Y[trait]

        # Convert to Dmatrix data structure 
        data_dmatrix = xgb.DMatrix(data=X,label=y)

        # Create train and test set
        X_train = X.loc[~X.index.isin(test[0])]
        X_test = X.loc[X.index.isin(test[0])]
        y_train = y.loc[~y.index.isin(test[0])]
        y_test = y.loc[y.index.isin(test[0])]

        ########## Parameter sweep using grid search ##########
        start = time.time()
        
        gs_results = pd.DataFrame(columns=['mean_test_score', 'params'])
        
        """GS_REPS = 10
        for j in range(GS_REPS):
            print("Round %s of %s" % (j + 1, GS_REPS))
            
            # Instantiate XGBoost regressor object
            xg_reg = xgb.XGBRegressor()#random_state=random.seed(123))
        
            # Hyperparameter tuning
            parameters = {"subsample":[0.5, 0.75, 1],
                        "colsample_bytree":[0.5, 0.75, 1],
                        "max_depth":[3, 10],
                        "min_child_weight":[1, 5],
                        "learning_rate":[0.3, 0.1, 0.03],
                        "n_estimators":[100, 500]}
            gs = GridSearchCV(estimator=xg_reg,
                        param_grid=parameters,
                        cv=5,           # cross validation folds
                        scoring='r2',   # find model with the best R-squared
                        verbose=1)
            fitted_model = gs.fit(X_train, y_train) # fit the model
            
            # Add results to dataframe
            j_results = pd.DataFrame(gs.cv_results_)
            gs_results = pd.concat([gs_results, j_results[['params',
                'mean_test_score']]])
            
        # Break params into seperate columns and save to file
        gs_results2 = pd.concat([gs_results.drop(['params'], axis=1),
            gs_results['params'].apply(pd.Series)], axis=1)
        gs_results.to_csv("%s_%s_GridSearchFULL.txt"%(data_type, trait))

        print("Best parameters: ", fitted_model.best_params_)
        """#######################################################
        # Run XGBoost model
        results_val = []
        results_test = []
        n = 1
        for j in range(0, n):
            print("Running %i of %i" % (j + 1, args.n))
            #best_model = gs.best_estimator_
			best_model =  xgb.XGBRegressor()

            # Cross-validation
            cv_pred = cross_val_predict(estimator=best_model, X=X_train, y=y_train,
                                        cv=5, n_jobs=-1)
            
            # Performance statistics from cross-validation
            mse_val = mean_squared_error(y_train, cv_pred)
            rmse_val = np.sqrt(mean_squared_error(y_train, cv_pred))
            evs_val = explained_variance_score(y_train, cv_pred)
            r2_val = r2_score(y_train, cv_pred)
            cor_val = np.corrcoef(np.array(y_train), cv_pred)
            print("Val RMSE: %f" % (rmse))
            print("Val R-sq: %f" % (r2))
            print("Val PCC: %f" % (cor_val))
            result_val = [mse_val, rmse_val, evs_val, r2_val, cor_val[0, 1]]
            results_val.append(result_val)

            # Fit the model
            best_model.fit(X_train, y_train)

            # Apply the model to the test set
            test_pred = best_model.predict(X_test)

            # Performance on the test set
            mse = mean_squared_error(y_test, test_pred)
            rmse = np.sqrt(mean_squared_error(y_test, test_pred))
            evs = explained_variance_score(y_test, test_pred)
            r2 = r2_score(y_test, test_pred)
            cor = np.corrcoef(np.array(y_test), test_pred)
            print("Test RMSE: %f" % (rmse))
            print("Test R-sq: %f" % (r2))
            print("Test PCC: %f" % (cor))
            result_test = [mse, rmse, evs, r2, cor[0, 1]]
            results_test.append(result_test)

            run_time = time.time() - start
            print("Run Time: %f" % (run_time))

        # Save the fitted model
        filename = "%s_%s_model.save"%(data_type, trait)
        pickle.dump(best_model, open(filename, 'wb'))

        # Save results to file
        if not os.path.isfile("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/ \
            Genomic_Prediction_XGBoost/RESULTS_xgboost.txt"):
            out = open("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/ \
                Genomic_Prediction_XGBoost/RESULTS_xgboost.txt", "a")
            out.write("Date\tRunTime\tData\tTrait\tMSE_val\tRMSE_val\tEVS_val \
                \tR2_val\tPCC_val\tMSE_test\tRMSE_test\tEVS_test\tR2_test \
                \tPCC_test\n")
            out.close()

        out = open("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/ \
            Genomic_Prediction_XGBoost/RESULTS_xgboost.txt", "a")
        out.write("%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%
        (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),run_time,data_type,trait,
            mse_val,rmse_val,evs_val,r2_val,cor_val,mse,rmse,evs,r2,cor))

        # Feature importances
        #importances = best_model.feature_importances_
if __name__ == "__main__":
    main()
