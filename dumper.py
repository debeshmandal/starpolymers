# analysis of lamps dump files

import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances

def com(df):
    """
    returns a list of arrays [x_mean, y_mean, z_mean, mol]
    """
    result = []
    n_mols = df['mol'].max()
    for i in range(n_mols):
        data = df[df['mol']==i+1]
        temp = np.ndarray([1,4])
        temp[0][0] = data['x'].mean()
        temp[0][1] = data['y'].mean()
        temp[0][2] = data['z'].mean()
        temp[0][3] = data['mol'].mean()
        result.append(temp)
    return result

def distances(r1, r2):
    dif = r1-r2
    dif = np.square(dif)
    sum_dif =dif[:,0]+dif[:,1]+dif[:,2]
    dist = np.sqrt(sum_dif)
    return dist
    

class DumpReader():

    def __init__(self, ID):
        self.ID = str(ID)

    def read(self, step, kind='positions-short'):
        
        fname = 'results/{0}/dump.{0}.{1}'.format(self.ID, step)
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                if i == 8:
                    columns = line
        columns = columns.lstrip('ITEMS: ATOMS ')
        columns = columns.split()
        positions = pd.read_csv(fname, skiprows=9, delimiter = ' ',
                                header=None)
        delete = len(positions.columns.tolist())-1
        positions = positions.drop(delete, axis=1)
        positions = positions.set_axis(columns, axis='columns', inplace=False)

        if kind == 'all':
            data = positions

        elif kind == 'positions-short':
            data = positions[['x', 'y', 'z', 'mol']]

        elif kind == 'positions-long':
            data = positions[['id', 'x', 'y', 'z', 'mol', 'q']]
            
        return data

    def com(self, molecule, steps):
        
        result = pd.DataFrame()
        for i in steps:
            data = self.read(i)
            centres = com(data)
            index = int(molecule)-1
            temp = centres[index]
            temp_df = pd.DataFrame(temp)
            temp_df['Step'] = i
            temp_df = temp_df.drop(3, axis=1)
            result = result.append(temp_df, ignore_index=True)
        return result        

    def complex_distance(self, steps):

        d1 = self.com('1', steps)[[0,1,2]].values
        d2 = self.com('2', steps)[[0,1,2]].values
        data = pd.DataFrame(distances(d1,d2))
        data['Step'] = steps
        return data
            
        
        
    
        
