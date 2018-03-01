## Machine learning testing script

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import posanal as usrpos
import random

data_init = pd.read_csv('results_exp1-5.csv')
# ---- Data Pre-processing ---- #
dataset = data_init
X = dataset.iloc[:, :-1].values
y = dataset.iloc[:,6].values

print dataset
print '\nX = '
print X
print '\nY = '
print y

    # ---- Split into training and testing ---- #

from sklearn.cross_validation import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    test_size = 0.2,
                                                    random_state = 0)
    # ---- Least Squares Method ---- #

from sklearn.linear_model import LinearRegression
regressor = LinearRegression()
regressor.fit(X_train, y_train)

y_pred = regressor.predict(X_test)

print 'y_pred = \n', y_pred
print '\ny_test = \n', y_test

    # ---- Plot this test ---- #

plt.plot(y_test, (y_test-y_test), c='red', label = 'Test Data')
plt.scatter(y_test, (y_test-y_pred)/y_test, c='blue', label = 'Predicted Data')
plt.xlabel('$R_g$ (test)',size='large')
plt.ylabel('$R_g$ (test, predicted)',size='large')
plt.legend()
plt.show()

# Apply model

# data_ML

# Use model to generate results for all kap and lam

# plot everything

# data_actual
