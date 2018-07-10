# Conducting Free Energy Calculations on
# Umbrella sampled simulations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from math import exp, log
from pandas import DataFrame as df


kBT = 1.2
bf = exp(-1.0/kBT) # e^-kBT - same as exp(beta)


def write_to_csv(path):
    data = pd.read_csv('{}/xi.hist'.format(path), header=None, delim_whitespace=True,
                       skiprows=4).rename(columns={
                           0: 'bin',
                           1: 'xi',
                           2: 'counts',
                           3: 'p'})
        
    data.to_csv('{}/xi.csv'.format(path), index=False)

def compute_to_csv(wham_obj, centre, p=0.02):
    data = wham_obj.get_data(centre)
    data = data[data['p']>p].reset_index(drop=True)
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
    
    result['F_bias'] = (-kBT * np.array([math.log(num) for num in dataframe['p'].values]))
    
    return result
    

def wham_F(A, xis, centres):

    ## log(c(z)) - log(n) - log(sum{exp[w_j(z) + BAj})

    # 
    #

    # sum over all j for each xi value, return two column array
    # where first column is xi (colvar) and second column is F

    # need an array of Aj values

    
    
    
    return F_xis # this should be an array

def wham_A(F, centres, xis):
    
    ## -BAj = ln sum_z{e^[-BF(z) + w_j(z)]}

    # sum over all xis for each j value, return two column array
    # where first column is j (center) and second column is F

    # J is the same as range(len(xis))
    
    summing_term = exp(-1/kBT * wham_F(xis) + bias(centres, xis))
    BAj = - log()


    # return a DataFrame with two columns
    return A_centers

def boltzmann(x, bf):
    return math.exp()

class WHAM():
    def __init__(self, ID, centres, K=10):
        self.master = pd.DataFrame()
        self.ID = ID
        self.root = 'results'
        self.path = '{}/{}'.format(self.root, ID)
        self.hists = list()
        self.master = {'A': pd.DataFrame().rename(columns={0:'j', 1:'A'}),
                       'F': pd.DataFrame().rename(columns={0:'xi', 1: 'F'})}
        self.columns = list()
        self.biases = list()
        self.centres = centres # e.g. [1, 2, 3, 4...]
        self.xis = data = pd.read_csv('{}_{}/xi.hist'.format(self.path, centres[0]), header=None, delim_whitespace=True,
                       skiprows=4).rename(columns={
                           0: 'bin',
                           1: 'xi',
                           2: 'counts',
                           3: 'p'})[['xi']]
        self.K = K
      

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
        return False

    def iterate(self):

        """

        Perform the WHAM iterations until complete
        
        """

        # iterate previously defined functions until a condition is met

        A = self.master['A']
        F = self.master['F'] 

        ## run equation 1

        A = wham_A(A)
        F = wham_F(F)

        ## run equation 2

        ## calculate convergence condition

        F_old = F
        F_new = F
        convergence = F_old-F_new


        ## see if it is low enough

        if converging(convergence) != True:

            self.iterate()
            
        else:

        ## if not then repeat

        ## if it is then return final dataframe

            return self.master

    def full_run(self, max_iterations=100, conv=0.01):
        self.master = self.initialise()
        
        for i in range(max_iterations):
            
            if np.mean(convergence['F']) > conv:
                result = self.iterate()
                convergence = self.master['F'] - result['F']
                self.master['F'] = result['F']
                
                # print np.mean(convergence['F'])
            else:
                return [np.mean(convergence['F']), i]
        
