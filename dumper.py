# analysis of lamps dump files

import pandas as pd
import numpy as np

def distances(r1, r2, box=50):
    dist = np.abs(r1-r2)
    dist = np.where(dist>0.5*box, dist-box, dist)
    dist = np.power(dist, 2)
    dist = np.sum(dist)
    dist = np.power(dist, 0.5)
    return dist

def delta(molecule, salt):
    cds = []
    pos = molecule[['x','y','z']].values
    #print 'molecule is: \n', pos
    cis = salt[['x','y','z']].values
    #print 'salt is: \n', cis
    for i in range(len(pos)):
        dist = []
        for j in range(len(salt)):
            dist.append(distances(pos[i], cis[j]))
        cds.append(min(dist))
            
    return cds

class DumpReader():

    def __init__(self, ID, box=50):
        self.ID = str(ID)
        self.root = 'results'
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

    def com(self, molecule):
        fname = self.path + '/com{}.out'.format(molecule)
        data = pd.read_csv(fname,
                           comment='#',
                           delim_whitespace = True,
                           header=None)
        com = pd.DataFrame()
        com['ts'] = data[data[1]==3][0].values
        com['x'] = data[data[0]==1][1].values
        com['y'] = data[data[0]==2][1].values
        com['z'] = data[data[0]==3][1].values
        
        return com
        

    def complex_distance(self, steps):

        d1 = self.com(1)[[0,1,2]].values
        d2 = self.com(2)[[0,1,2]].values
        data = pd.DataFrame(distances(d1,d2))
        data['Step'] = steps
        return data

    def counterion_distance(self, step, bound=20):

        data = self.read(step, kind='positions-long')
        pairs = [['star','neg'],['dna','pos']]
        dfs = {'star': data[data['mol']==1],
               'dna': data[data['mol']==2],
               'pos': data[(data['q']>0) & (data['mol']==3)],
               'neg': data[(data['q']<0) & (data['mol']==3)]}

        
        coms = {'star': self.com(1)[self.com(1)['ts']==step],
               'dna': self.com(2)[self.com(2)['ts']==step]}


        # get all salt ions from within bound from relevant com

        # rescale relevant dfs by centre of mass
        distances = dict()

        distances['star'] = delta(dfs['star'], dfs['neg'])
        #print distances['star']
        distances['dna'] = delta(dfs['dna'], dfs['pos'])
        #print distances['dna']

        # concat both above dataframes

        result = pd.concat([pd.DataFrame(distances['star']),
                            pd.DataFrame(distances['dna'])])
        return np.mean(result.values)
        #result.std()

        # take mean and std

        # return mean and std

        
        

        
            
        
        
    
        
