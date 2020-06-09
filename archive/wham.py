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

def pprint(a):
    print(('\n', pd.DataFrame(a), '\n'))

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

    result = np.log(array)
    result = -kBT * result
    
    return result
    

def wham_F(p, w):

    """

    Returns F_i which is a 1D array with length that is equal to the number of centres

    """
    F_i = list()
    for i in range(w.shape[1]):
        F_i.append(np.sum(np.multiply(w[:,i],p), axis=0))
    F_i = np.array(F_i)
    return F_i

def wham_p(counts, F, w):

    """

    Returns p which is a vector with length the number xis

    """
    
    
    numer = np.sum(counts, axis=1)
    N = np.sum(counts,axis=0)
    wf = np.multiply(1.0/F, 1.0/w)
    denom = np.sum(np.multiply(1.0/N,1.0/wf), axis=1)
    pprint(denom)
    

    p = denom*numer
    #pprint(p)
    
    return p

class WHAM():
    def __init__(self, ID, centres, K=2, root='results', max_counts = 5000):
        self.ID = ID
        self.centres = centres # e.g. [1, 2, 3, 4...]
        self.K = K
        self.max_counts = max_counts
        self.root = root
        self.path = '{}/{}'.format(self.root, ID)
        self.xis = pd.read_csv('{}_{}/xi.hist'.format(self.path, centres[0]), header=None, delim_whitespace=True,
                               skiprows=4).rename(columns={
                                   0: 'bin',
                                   1: 'xi',
                                   2: 'counts',
                                   3: 'p_bias'})[['xi']]
        self.master = {'counts': pd.DataFrame(columns=self.centres),
                       'master':pd.DataFrame(columns=['xi', 'p']),
                       'exp(-bF)': pd.DataFrame(columns=['centre', 'F', 'exp(-bF)']),
                       'exp(bw)': pd.DataFrame(columns=self.centres+['xi'])}
        self.master['counts']['xi'] = self.xis['xi']
        self.master['master']['xi'] = self.xis['xi']
        self.master['exp(-bF)']['centre'] = self.centres
        self.master['exp(bw)']['xi'] = self.xis['xi']
      

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
            self.master['counts'][i]=data['counts']
            self.master['exp(bw)'][i] = bias(self.xis.values, i, self.K)

        # initialise F_i - set all values to 1 and return exp(F)*bf

        F = np.ones([len(self.master['exp(-bF)']),1])
        self.master['exp(-bF)']['F'] = F
        self.master['exp(-bF)']['exp(-bF)'] = np.exp(F) * bf

    def raw_histogram(self):
        fig, ax = plt.subplots()
        ax.set_prop_cycle('color',
                          plt.cm.viridis(np.linspace(0,1,len(self.centres))))
        for centre in self.centres:
            plt.plot(self.master['counts']['xi'],
                     self.master['counts'][centre], alpha=0.5)
        plt.plot(self.master['counts']['xi'],
                 self.master['counts'][self.centres].sum(axis=1),
                 'k:')
        plt.ylabel('counts')
        plt.xlabel(r'$xi$ [$sigma$]')
        plt.show()
        

    def iterate(self):

        """

        Perform the WHAM iterations until complete
        
        """

        # iterate previously defined functions and return difference

        counts = self.master['counts'][self.centres].values
        F = self.master['exp(-bF)']['exp(-bF)'].values # exponentials
        w = self.master['exp(bw)'][self.centres].values

        ## run WHAM equations

        p = wham_p(counts, F, w)
        new_F = wham_F(p, w)        

        ## calculate convergence condition        
        
        self.master['master']['p'] = p
        self.master['exp(-bF)']['exp(-bF)'] = new_F


    def full_run(self, max_iterations=100, conv=0.001, lower_bound=0.001):

        self.initialise_csv_files()
        self.initialise_master()
        print(((self.master['counts'][self.centres].sum(axis=1)/self.master['counts'][self.centres].sum(axis=1).sum()).sum()))
        F_old = np.copy(self.master['exp(-bF)']['exp(-bF)'].values)
        self.iterate()
        F_new = self.master['exp(-bF)']['exp(-bF)'].values
        convergence = abs(np.mean(F_new - F_old))
        counter = 0

        while (counter < max_iterations) & (abs(convergence) > conv):
            print(('\nconvergence = ', convergence))
            print(('\ntotal probability = ',self.master['master']['p'].sum()))
            F_old = np.copy(self.master['exp(-bF)']['exp(-bF)'].values)
            self.iterate()
            F_new = self.master['exp(-bF)']['exp(-bF)'].values
            convergence = abs(np.mean(F_new - F_old))
            counter+=1
                

        self.master['master'] = self.master['master'][self.master['master']['p']>lower_bound]
        p = self.master['master']['p'].values
        self.master['master']['PMF'] = free_energy(p)
            
    def merge(self, merger='counts', plot='off'):

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
            plt.xlabel('xi')
            plt.ylabel('{}_merged'.format(merger))
        return results
