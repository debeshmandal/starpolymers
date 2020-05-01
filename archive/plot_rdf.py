# plot g(r)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

kap = 10
lam = [11]

for j in lam:
    fname = 'rdfdata.csv'
    data=pd.read_csv(fname, header=None)
    print data
plt.scatter(data[0],data[1],label=('$\lambda=$'+str(j)))
plt.legend()
plt.show()
