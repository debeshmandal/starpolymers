# plot rg_mean.csv

import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('rg_mean.csv')
axesrange = [3,4,5,6,7,8,9,10,11]
data = data.set_index([axesrange])
data = data.transpose()
print(data)
for i in axesrange:
    if i % 2 == 0:
        plt.scatter(axesrange, data[i], label=('$\lambda = $'+str(i)))
plt.xlabel('$\kappa$')
plt.ylabel('$R_g$')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102), loc=3, ncol=4, mode = "expand",borderaxespad=0.)
plt.show()
