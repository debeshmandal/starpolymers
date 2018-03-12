## Machine learning testing script

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import posanal as usrpos
import random


plot_ML = True

data_init = pd.read_csv('results_exp1-5.csv')
# ---- Data Pre-processing ---- #
dataset = data_init
X = dataset.iloc[:, :-1].values
y = dataset.iloc[:,6].values

#print dataset
#print '\nX = '
#print X
#print '\nY = '
#print y

    # ---- Split into training and testing ---- #
results_RF = pd.DataFrame()
results_SLR = pd.DataFrame()
test_size = np.arange(0.05,1,0.05)
print test_size

for i in range(19):
    ts=float(i+1)/20
    from sklearn.cross_validation import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size = ts,
                                                        random_state = 10)
        # ---- Regression Model ---- #

    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import LinearRegression
    regressor = RandomForestRegressor()
    regressor.fit(X_train, y_train)
    regressor2 = LinearRegression()
    regressor2.fit(X_train, y_train)

    y_pred = regressor.predict(X_test)
    y_pred2 = regressor2.predict(X_test)

    # print 'y_pred = \n', y_pred
    # print '\ny_test = \n', y_test

        # ---- Plot this test ---- #

    if plot_ML == True:
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(y_test, y_test, 'k--', label = 'Test Data')
        ax.scatter(y_test, y_pred, c='blue', label = 'Random Forest Prediction', s=10)
        ax.scatter(y_test, y_pred2, c='red', label = 'Linear Regression Prediction',s=10)
        ax.set_xlabel('$R_g$ (test)',size='large')
        ax.set_ylabel('$R_g$ (test, predicted)',size='large')
        ax.legend()
        ax.set_title('Test Size = {}'.format(ts))
        fname = str('prediction_{}.png'.format(ts))
        plt.savefig(fname)

        # ---- Evaluate results ---- #

    def error(pred, test):
        N = len(test)
        err_df = np.ndarray([N,1])
        for i in range(N):
            err_df[i][0] = abs((pred[i]-test[i]))/test[i]

        err_df=pd.DataFrame(err_df)
        mean_err = err_df.mean(axis=0)
        std_err = err_df.std(axis=0)
        errors = np.array([mean_err,std_err])
        #print ts, errors, '\n'
        return errors

    err = 100*error(y_pred, y_test).transpose()
    err2 = 100*error(y_pred2, y_test).transpose()
    temp_RF = pd.DataFrame(err)
    temp_SLR = pd.DataFrame(err2)
    results_RF = results_RF.append(temp_RF,ignore_index=True)
    results_SLR = results_SLR.append(temp_SLR,ignore_index=True)

print results_RF
print results_SLR
fig = plt.figure()
ax = fig.gca()
ax.errorbar(test_size, results_RF[0],yerr=results_RF[1],
             c='red',
             label='Random Forest Error',
             fmt='o:',
             capthick=2)

ax.errorbar(test_size, results_SLR[0],yerr=results_SLR[1],
             c='blue',
             label='Linear Regression Error',
             fmt='s:',
             capthick=2)
ax.xlabel('Test Size')
ax.ylabel('Error (%)')
ax.legend()
plt.show()
