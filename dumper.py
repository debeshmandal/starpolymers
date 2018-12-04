# analysis of lamps dump files

import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, euclidean
from numpy.linalg import norm
import math

def distances(r1, r2, box=50):
    dist = np.abs(r1-r2)
    for i in range(3):
        dist[i] = dist[i] - box*round(dist[i]/box)
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

def com(array):
    return np.mean(array, axis=0)

class DumpReader():

    def __init__(self, ID, box=50):
        self.ID = str(ID)
        self.root = 'results'
        self.path = '{}/{}'.format(self.root, self.ID)
        self.box = box

    def change_path(self, path, kind='path'):
        if kind=='path':
            self.path = path
        elif kind=='root':
            self.root = path

    def read(self, step, kind='positions-short', unwrap=False):
        
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
        if unwrap == True:
            try:
                positions = positions.rename(columns={'xu': 'x',
                                                      'yu': 'y',
                                                      'zu': 'z'})
                positions['x'] = positions['x'].values % self.box
                positions['y'] = positions['y'].values % self.box
                positions['z'] = positions['z'].values % self.box
            except:
                None

        if kind == 'all':
            data = positions

        elif kind == 'positions-short':
            data = positions[['id',
                              'x',
                              'y',
                              'z',
                              'mol']]


        elif kind == 'positions-long':
            data = positions[['id',
                              'x',
                              'y',
                              'z',
                              'mol', 'q']]
           
        return data

    def read_com(self, molecule):
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

        
        #coms = {'star': self.com(1)[self.com(1)['ts']==step],
        #       'dna': self.com(2)[self.com(2)['ts']==step]}


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

    def counterion_number(self, step, bound=2):
        data = self.read(step, kind='positions-long')
        dfs = {'star': data[data['mol']==1].reset_index(drop=True),
               'dna': data[data['mol']==2].reset_index(drop=True),
               'pos': data[(data['q']>0) & (data['mol']==3)].reset_index(drop=True),
               'neg': data[(data['q']<0) & (data['mol']==3)].reset_index(drop=True)}
        count = 0

        dist = np.transpose(distance_matrix(dfs['star'][['x','y','z']].values,
                                    dfs['neg'][['x','y','z']].values))

        count += len(dist[dist[:,0]<=bound])

        dist = np.transpose(distance_matrix(dfs['dna'][['x','y','z']].values,
                                    dfs['pos'][['x','y','z']].values))

        count += len(dist[dist[:,0]<=bound])
        
        #for i in range(len(dfs['pos'])):
        #    ion=dfs['pos'][['x', 'y', 'z']].loc[[i]].values
        #    for j in range(len(dfs['dna'])):
        #        atom = dfs['dna'][['x', 'y', 'z']].loc[[j]].values
        #        
        #        if distances(ion, atom)<bound:
        #            count +=1
        #            continue

        

                    
        #for i in range(len(dfs['neg'])):
        #    ion = dfs['neg'][['x', 'y', 'z']].loc[[i]].values
        #    for j in range(len(dfs['star'])):
        #        atom = dfs['star'][['x', 'y', 'z']].loc[[j]].values
        #        #count_2 = 0
        #        if distances(ion, atom)<bound:
        #            count +=1
        #            continue

        #count = count_1 + count_2

        return count

    def xyz(self, steps, out=True, folder=None):
        for step in steps:
            data = self.read(step)
            data = data[['mol', 'x', 'y', 'z']]
            if folder != None:
                fname = '{0}/{1}.{2}.xyz'.format(folder, self.ID, step)
            else:
                fname = '{0}.{1}.xyz'.format(self.ID, step)
            with open(fname, 'w') as f:
                f.write(str(data.shape[0]-1))
                f.write('\n')
                data.to_csv(f, header=False, index=None, sep=' ')

    def dep_pqr(self, config, steps, out=True, folder=None, radius=2.0,
            molecules=[1,2]):
        
        for step in steps:
            data = self.read(step, kind='positions-long')
            r = radius
            # write HETATMs
            master = data.rename(columns={'id': 'residueNumber',
                                          'x': 'X',
                                          'y': 'Y',
                                          'z': 'Z',
                                          'mol': 'molName',
                                          'q': 'charge'})
            master['recordName'] = ['HETATM']*len(master)
            master['radius'] = [radius] * len(master)
            atomNames = ['C'] * len(master)
            master['atomName'] = atomNames
            #master['molName'] = ['ION'] * len(master)
            master = master[master['molName'] != 3]
            master = master[['recordName', 'residueNumber',
                             'atomName', 'molName',
                             'X',
                             'Y',
                             'Z',
                             'charge',
                             'radius']]
            
            
            # write CONECTs

            conects = pd.DataFrame()

            # read from BONDS to ANGLES

            # generate list of atom_IDs - this becomes ATOM_ID

            # for all atoms (get total number of atoms)

                # ATOM 1 = each atom - create line

                # search column 1 - if atom appears in this ATOM $(counter+2)
                # equals number in column 2, counter += 1

                # search column 2 - if atom appears in this, ATOM $(counter+2)
                # equals number in column 1, counter += 1, if counter == 3,
                # new line
            
            if folder != None:
                fname = '{0}/{1}.{2}.pqr'.format(folder, self.ID, step)
            else:
                fname = '{0}.{1}.pqr'.format(self.ID, step)

            # with open as f:
            master = master.round(2)
            master.to_csv(fname, header=False, index=None, sep=' ')
                # conects.to_csv
                # END
            print master

    def gyration(self, step, salt=3):
        """

        Returns the radius of gyration of all molecules except the salt
        molecule.

        Arguments:

        step: integer value for timestep
        salt=3: which molecule to exclude

        """
        
        # read file ensuring that ix iy and iz are present
        try:
            data = self.read(step, kind='all')[['mol','x','y','z',
                                                'ix', 'iy', 'iz']]
        except KeyError:
            print "Periodic Boundary image (i.e. ix, iy, iz) not present"
        
        # remove salt atoms
        data = data[data['mol']!=salt]

        # for atoms that are not in i*==0, move them by box*i*
        data['xs'] = data['x'].values + self.box*data['ix'].values
        data['ys'] = data['y'].values + self.box*data['iy'].values
        data['zs'] = data['z'].values + self.box*data['iz'].values
        data = data[['xs', 'ys', 'zs']]

        dist = 0
        for i in range(len(data)):
            for j in range(len(data)):
                dist += distances(data.values[i], data.values[j]) ** 2
            rg2 = dist/(2*len(data)**2)
            rg = math.sqrt(rg2)
        return rg
        
    
        
