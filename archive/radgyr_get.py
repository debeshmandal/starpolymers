# module to read the LAMMPS fix ave/time gyration output
# First 3 lines are to be ignored
# From here we can do a few things
# The most useful would be to take an average of the Rg values
# over the time series

# since the format of the rest of the file is on alternating lines

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def gyrstats(kap,lam):
    fname = 'results/star'+str(kap)+'_'+str(lam)+'/gyr.out'
    x=pd.read_csv(fname, delimiter=' ', skiprows=3)
    x=x['1']
    x=x[x!=1.0]
    x=x.reset_index(drop=True)
    mean = x.mean()
    std = x.std()
    return mean, std

results_ML = pd.DataFrame()
results_mean = pd.DataFrame(index=list(range(3,12)),columns=list(range(3,12)))

for kap in range(3,12):
    for lam in range(3,12):
        results_mean[kap][lam] = gyrstats(kap,lam)[0]
        temp = np.ndarray([1,4])
        temp[0][0] = int(kap)
        temp[0][1] = int(lam)
        temp[0][2] = gyrstats(kap,lam)[0]
        temp[0][3] = gyrstats(kap,lam)[1]
        temp_df = pd.DataFrame(temp)
        results_ML = results_ML.append(temp_df, ignore_index=True)

results_ML.to_csv('rg_ML.csv', index=False)
results_mean.to_csv('rg_mean.csv', index=False)

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
a = np.arange(3,12)
b = np.arange(3,12)
A, B = np.meshgrid(a,b)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel(r'Number of Arms ($\kappa$)', size='large')
ax.set_ylabel(r'Length of Arms ($\lambda$)', size='large')
ax.set_zlabel(r'$R_g$', size='large')
surf = ax.plot_surface(A,B,results_mean,cmap=cm.coolwarm,linewidth=0)
fig.colorbar(surf, shrink=0.5, aspect =5)
plt.show()
