#nc.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pmf
from pmf import PMF
from design import *
from rdf import RDF
from starpolymers.dumper import DumpReader

def _get_min(PMF, bound=5):
    return PMF.min(bound=bound) # returns range of xi as tuple

def _label_generator(variables, parameters, units=None):
    dictionary = dict()
    for i in range(len(parameters)):
        label = ''
        params = parameters[i]
        for j in range(len(variables)):
            label += variables[j]
            label += '='
            label += str(params[j])
            if units != None:
                    label += units[j]
            if j+1 != len(variables):
                label += ', '
        dictionary[i+1] = label
    
    return dictionary

def _manning_radius(f, mol1, mol2, q2, L=100):

    atoms = DumpReader('none', 
                       box=50,
                       fname=f).read(use_fname=True,
                                     kind='positions-long')

    molecule = atoms[atoms['mol']==mol1][['x','y','z']].values
    counterions = atoms[atoms['mol']==mol2]
    counterions = counterions[counterions['q']==q2][['x','y','z']].values
    myrdf = rdf.RDF(molecule, counterions, L=L, rmax=10.0, bin_width=0.5)
    radius = myrdf.condensed_radius
    return radius

def _nc(ts, radius, mol1, mol2, q2, L=100, fdump='dump.{}.lammpstrj'):
    results = []
    for i in ts:
        fname = fdump.format(i)
        atoms = DumpReader('none', 
                           box=50,
                           fname=f).read(use_fname=True,
                                         kind='positions-long')

        molecule = atoms[atoms['mol']==mol1][['x','y','z']].values
        counterions = atoms[atoms['mol']==mol2]
        counterions = counterions[counterions['q']==q2][['x','y','z']].values
        results.append(rdf.RDF(molecule, counterions, L=100, rmax=10.0, bin_width=0.5, condensed_radius=radius).NC)

    data = pd.DataFrame()
    data['ts'] = ts
    data[r'$N_C$'] = results
    return data
    

def _get_nc(runs, timesteps, radius, mol1, mol2, q2, root=None, 
            f_traj='out.colvars.traj'):
    """

    From a set of N runs starting at 1 and ending at N, return
    a DataFrame with columns=['xi', 'mean', 'std'], note that
    the gyr will be scaled by the last 5 values (i.e. towards the 
    Colvars targetCenters argument)

    Parameters:
    -----

    runs : int (positive)
           Number of different sampling runs

    fname : string
            name of gyr containing file

    
    Returns
    -----
     
    data : pd.DataFrame
           return DataFrame with columns=['xi', 'mean', 'std']

    """

    xi = pd.Series()
    nc = pd.Series()
    if root != None:
        root = '{}/'.format(root)
    for run in range(1, runs+1):
        fname = 'dump.{}.lammpstrj'
        fin = '{}{}/{}'.format(root, run, fname)
       
        try:

            temp = _nc(timesteps, radius, mol1, mol2, q2, L=L, 
                       fdump=fin)
            PMF_traj = '{}{}/{}'.format(root, run, f_traj)
            try:
                traj = pmf._get_traj(PMF_traj) # has columns ['ts', 'xi']
            except:
                None
            # merge temp and traj
            merged = pd.merge(temp, traj, on='ts')         
            
            xi=xi.append(merged['xi'])
            nc=nc.append(merged[r'$N_C$'])
            
        except IOError:
            print "Warning: run {} did not work!".format(run)
    
    data = pd.DataFrame() # should end with long dataframe with
                          # columns ['xi', 'NC']
    
    data['xi'] = xi
    data[r'$N_C$']=nc
    data = data.sort_values(by='xi')
    data = data[data['xi']<75]
    data = data.reset_index(drop=True)
    
    return data

class NC():
    """

    Class containing gyr data in DataFrame form with method to write
    to csv file

    Parameters
    -----

    fname : string
            The name of the file (without extension) to be used when
            writing the file using the write() method

    runs : int (positive)
           Number of different sampling runs

    Methods
    -----

    write() : writes file to '{self.fname}.csv'

    """
    def __init__(self, manning_f, PMF, mol1, mol2, q2, fname='gyr', runs=10, root=None, bound=5):
        if root != None:
            self.fname = '{}/{}.csv'.format(root, fname)
        else:
            self.fname = '{}.csv'.format(fname)
        self.condensed_radius = _get_manning_radius(manning_f)
        self.nc = _get_nc(runs, root=root)
    
    def write(self):
        _write_nc(self.nc, self.fname)
        return

class NC_AVERAGE():
    def __init__(self):
        self.fname = None