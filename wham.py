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
    data.to_csv('{}_{}/xi.csv'.format(wham_obj.path, centre), index=False)

def bias(array, centre, K=2):

    # initialise result DF

    result = df()
    result['xi'] = array

    # compute the bias force field
    
    coeff = 0.5 * K
    delta = array - centre
    result['u_bias'] = coeff * np.square(delta)
    
    return result

def free_energy(dataframe):
    result = df()
    result['xi'] = dataframe['xi']
    
    result['F_bias'] = (-kBT * np.array([math.log(num) for num in dataframe['p_bias'].values]))
    
    return result
    

def wham_F(A, w, shifts):

    ## for each xi:

    ## return the sum of a shift(centre) * 1/weight(xi, centre) * A(centre)

    w = np.mat(np.power(w, -1))
    A = np.mat(A)
    F = np.matmul(w, A) * shifts    
    
    return F_xis # this should be an array

def wham_A(F, w):
    
    ## for each centre:

    ## return the sum of [weight(xi, centre)*Free_energy(xi)] for all xis

    w = np.mat(w)
    F = np.mat(F)
    A = np.matmul(w, F)
    
    return A

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
        self.master = {'A': pd.DataFrame(columns=['centre', 'A','exp(-bA)']),
                       'F': pd.DataFrame(columns=['xi', 'F', 'exp(-bF)']),
                       'exp_u': np.zeros([len(self.xis), len(centres)])}
        self.shifts = np.ones([1, len(centres)])
      

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
        

    def initialise_master(self, K=2):
        """

        Read the first 

        """

        # calculate the weights

        for centre in self.centres:
            self.master['exp_u'][:,self.centres.index(centre)] = np.exp(bias(self.xis, centre, K)['u_bias'])

        # set the free energies to zero
        # calculate the exp(-bf) column
        
        F = np.zeros([len(self.xis), 1])
        self.master['F']['xi'] = self.xis
        self.master['F']['F'] = F[:,0]
        self.master['F']['exp(-bF)'] = np.exp(F) * bf

        # set Aj to zero
        # calculate the exp(-bA) column
        
        A = np.zeros([len(self.centres), 1])
        self.master['A']['centre'] = self.centres
        self.master['A']['A'] = A[:,0] 
        self.master['A']['exp(-bA)'] = np.exp(A) * bf
        return self.master

        

    def iterate(self, weights, shifts):

        """

        Perform the WHAM iterations until complete
        
        """

        # iterate previously defined functions until a condition is met

        bA = self.master['A']['exp(-bA)'].values # exponentials for all
        bF = self.master['F']['exp(-bF)'].values # exponentials for all

        ## run WHAM equations

        new_bA = wham_A(bF, weights)
        new_bF = wham_F(bA, weights, shifts)

        ## calculate convergence condition

        convergence = new_bF-bF
        
        self.master['A']['exp(-bA)'] = new_bA
        self.master['F']['exp(-bF)'] = new_bF

        return convergence


    def full_run(self, max_iterations=100, conv=0.01):
        
        self.master = self.initialise_master()
        u = self.master['exp_u']
        c = self.shifts

        # iterate once

        convergence = self.iterate(u, self.shifts)
        
        for i in range(max_iterations):
            
            if np.mean(convergence['F']) > conv:
                result = self.iterate(u, self.shifts)
                convergence = self.master['F'] - result['F']
                self.master['F'] = result['F']
                
                # print np.mean(convergence['F'])
            else:
                return [np.mean(convergence['F']), i]
            
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
            plt.xlabel('xi')
            plt.ylabel('{}_merged'.format(merger))
        return results
