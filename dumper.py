# analysis of lamps dump files

import pandas as pd
import numpy as np

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

def distances(r1, r2, box=50):
    dist = np.abs(r1-r2)
    dist = np.where(dist>0.5*box, dist-box, dist)
    dist = np.power(dist, 2)
    dist = np.sum(dist)
    dist = np.power(dist, 0.5)
    return dist
    

class DumpReader():

    def __init__(self, ID):
        self.ID = str(ID)
        self.root = 'results/'
        self.path = '{}/{}'.format(self.root, self.ID)

    def change_path(self, path, kind='path'):
        if kind=='path':
            self.path = path
        elif kind=='root':
            self.root = path

    def read(self, step, kind='positions-short'):
        
        fname = '{0}/dump.{1}.{2}'.format(self.path, self.ID, step)
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
        positions = positions.sort_values('id').reset_index(drop=True)

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

        
            
        
        
    
        
