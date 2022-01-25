# Description:
# Author: Kenia Segura Aba
# Reference: Jason Brownlee, Deep Learning with Time Series Forecasting, Machine Learning Mastery, Available from https://machinelearningmastery.com/machine-learning-with-python/, accessed August 3rd, 2021.

import sys, os, argparse, time
from datetime import datetime
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow import keras
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold

tf.__version__
keras.__version__

# read in data
pheno = pd.read_csv("~/Shiu_Lab/Project/pheno.csv", index_col=0)
geno = pd.read_csv("~/Shiu_Lab/Project/geno.csv", index_col=0)
X = geno.copy()
y = pheno.YPACETATE

# define testing and training sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=125)

# build sequential model    
def seq_model():
    model = keras.models.Sequential()
    model.add(keras.layers.Dense(64456, input_dim=64456, kernel_initializer="normal", activation="sigmoid"))
    model.add(keras.layers.Dense(20000, kernel_initializer="normal", activation="sigmoid"))
    model.add(keras.layers.Dense(20000, kernel_initializer="normal", activation="sigmoid"))
    model.add(keras.layers.Dense(10000, kernel_initializer="normal", activation="sigmoid"))
    model.add(keras.layers.Dense(6000, kernel_initializer="normal", activation="sigmoid"))
    model.add(keras.layers.Dense(2000, kernel_initializer="normal", activation="sigmoid"))
    model.add(keras.layers.Dense(1, kernel_initializer="normal", activation="sigmoid"))
    model.compile(loss="mean_squared_error", optimizer="adam", metrics=["mse"])
    return model

# evaluate the model using 5-fold cross validation
estimator = keras.wrappers.scikit_learn.KerasRegressor(build_fn=seq_model, epochs=10, batch_size=5)#, verbose=0)
kfold = KFold(n_splits=2)
results = cross_val_score(estimator, X_train, y_train, cv=kfold)
print("Results: %.2f (%.2f) MSE" % (results.mean(), results.std()))

# fit the sequential model on the training data
#model.fit(X_train, y_train, epochs=10)

# evaluate the sequential model
#_, accuracy = model.evaluate(X_test, y_test)

# make predictions on test data
#predictions = model.predict(X_test)
#rounded = [round(x[0]) for x in predictions]
