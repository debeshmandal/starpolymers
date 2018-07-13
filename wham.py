# Conducting Free Energy Calculations on
# Umbrella sampled simulations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from math import exp, log
from pandas import DataFrame as df


kBT = 1.2
beta = 1.0/kBT
bf = exp(-beta) # e^-kBT - same as exp(beta)


def write_to_csv(path):
    data = pd.read_csv('{}/xi.hist'.format(path), header=None, delim_whitespace=True,
                       skiprows=4).rename(columns={
                           0: 'bin',
                           1: 'xi',
                           2: 'counts',
                           3: 'p_bias'})
        
    data.to_csv('{}/xi.csv'.format(path), index=False)

def compute_to_csv(wham_obj, centre, p=0.02):
    data = wham_obj.get_data(centre)
    data = data[data['p_bias']>p].reset_index(drop=True)
    data['u_bias'] = bias(data['xi'].values, centre, K=wham_obj.K)['u_bias']
    
    data['F_bias'] = free_energy(data)['F_bias']
    data.to_csv('{}_{}/xi_mod.csv'.format(wham_obj.path, centre), index=False)

def bias(array, centre, K=2):

    # compute the bias force field
    
    coeff = 0.5 * K
    delta = array - centre
    result = coeff * np.square(delta)
    result = np.exp(-result) * bf
    
    return result

def free_energy(array):

    print array
    result = np.log(array)
    result = -kBT * result
    
    return result
    

def wham_F(p, w):

    """

    Returns F_i which is a 1D array with length that is equal to the number of centres

    """
    x = np.multiply(p,w)
    F_i = np.sum(x, axis=0)
    return F_i

def wham_p(p_b, F, w):

    """

    Returns p_u which is a vector with length the number xis

    """

    p_u = np.ndarray([len(p_b), len(F)])
    for i in range(len(F)):
        weight = w[:, i] * F[i]
        p_u[:,i] = np.multiply(p_b[:, i], weight)

    return p_u

def boltzmann(x, bf):
    return math.exp()

class WHAM():
    def __init__(self, ID, centres, K=10, root='results'):
        self.ID = ID
        self.centres = centres # e.g. [1, 2, 3, 4...]
        self.K = K
        self.root = root
        self.path = '{}/{}'.format(self.root, ID)
        self.xis = pd.read_csv('{}_{}/xi.hist'.format(self.path, centres[0]), header=None, delim_whitespace=True,
                               skiprows=4).rename(columns={
                                   0: 'bin',
                                   1: 'xi',
                                   2: 'counts',
                                   3: 'p_bias'})[['xi']] 
        self.master = {'p_unbias': pd.DataFrame(columns=self.centres+['xi']),
                       'p_bias': pd.DataFrame(columns=self.centres+['xi']),
                       'master':pd.DataFrame(columns=['xi', 'p']),
                       'exp(-bF)': pd.DataFrame(columns=['centre', 'exp(-bF)']),
                       'exp(bw)': pd.DataFrame(columns=self.centres+['xi'])}
        self.master['p_bias']['xi'] = self.xis
        self.master['p_unbias']['xi'] = self.xis
        self.master['master']['xi'] = self.xis
        self.master['exp(-bF)']['centre'] = self.centres
        self.master['exp(bw)']['xi'] = self.xis
      

    # initialise a master dataframe that will be updated each iteration

    def change_path(self, path, keyword='path'):
        if keyword == 'path':
            self.path = path
        if keyword == 'root':
            self.root = path
            self.path = '{}/{}'.format(self.root, self.ID)

    def initialise_csv_files(self):
        """

        Get all data and write to .csv files on disk to
        reduce memory usage
        
        """
        for i in self.centres:
            path = '{}_{}'.format(self.path, i)
            write_to_csv(path)

    def add_to_data(self, centre=None, do_all=True):

        # create for loop

        # for each centre, calculate the biased free energy
        # also calculate the biasing potential
        # Calculate F and the biasing potential as a function of xi
        # then rewrite to csv

        if do_all == True:
            for centre in self.centres:
                compute_to_csv(self, centre)               
        else:
            if centre != None:
                compute_to_csv(self, centre)
            else:
                return "Error, set centre to a value using it as a keyword argument"

    def get_data(self, centre):

        """

        Read data from .csv files for each xi
        
        """

    # get data for a given value of xi
    
        fname = '{}_{}/xi.csv'.format(self.path, centre)
        data = pd.read_csv(fname)
        return data
        

    def initialise_master(self):
        """

        Read the first 

        """

        # get p_bias from files
        for i in self.centres:
            data = self.get_data(i)
            self.master['p_bias'][i]=data['p_bias']
            self.master['exp(bw)'][i] = bias(self.xis.values, i, self.K)

        # initialise F_i - set all values to 1 and return exp(F)*bf

        F = np.ones([len(self.master['exp(-bF)']),1])
        self.master['exp(-bF)']['exp(-bF)'] = np.exp(F) * bf
        

    def iterate(self):

        """

        Perform the WHAM iterations until complete
        
        """

        # iterate previously defined functions and return difference

        p_b = self.master['p_bias'][self.centres].values
        F = self.master['exp(-bF)']['exp(-bF)'].values # exponentials
        w = self.master['exp(bw)'][self.centres].values

        ## run WHAM equations

        p_u = wham_p(p_b, F, w)
        new_F = wham_F(p_u, w)

        

        ## calculate convergence condition

        p = np.sum(p_u, axis=1)
        
        self.master['p_unbias'][self.centres] = p_u
        self.master['master']['p'] = p
        self.master['exp(-bF)']['exp(-bF)'] = new_F


    def full_run(self, max_iterations=100, conv=0.01, lower_bound=0.001):

        self.initialise_csv_files()
        self.initialise_master()
        self.iterate()
        master = self.master['master']
        p_old = master['p'].values
        self.iterate()
        p_new = master['p'].values
        convergence = np.mean(p_new - p_old)
        counter = 0
        
        while convergence > conv:
            while counter < max_iterations:
                p_old = p_new
                self.iterate()
                p_new = master['p']
                convergence = np.mean(p_new - p_old)

        self.master['master'] = self.master['master'][self.master['master']['p']>lower_bound]
        p = self.master['master']['p'].values
        self.master['master']['PMF'] = free_energy(p)
            
    def merge(self, merger='p_bias', plot='off'):

        results = pd.DataFrame()
        results['xi'] = self.xis
        merge_list = []
        for i in self.centres:
            name = '{}_{}'.format(merger, i)
            data = pd.read_csv('{}_{}/xi.hist'.format(self.path, i), header=None, delim_whitespace=True,
                               skiprows=4).rename(columns={
                                   0: 'bin',
                                   1: 'xi',
                                   2: 'counts',
                                   3: 'p_bias'})[['xi', merger]].rename(columns={
                                       merger: name})
            
            results = pd.merge(results, data, on='xi')
            merge_list.append(name)
        results['{}_merged'.format(merger)] = results[merge_list].sum(axis=1)
        results = results[['xi','{}_merged'.format(merger)]]
        if plot == 'on':
            plt.bar(results.iloc[:,0], results.iloc[:,1], width=0.35)
            plt.xlabel(r'$\xi$ [$\sigma$]')
            plt.ylabel('{}_merged'.format(merger))
        return results
