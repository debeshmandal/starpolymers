# aim is to read, write and plot gyr files and DeltaG
# together

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pmf
from pmf import PMF

markers = ['o','d','*','P','h','x','+']

def _get_min(PMF, bound=1):
    return PMF.min(bound=bound) # returns range of xi as tuple

def _collate(gyr_list, variables, parameters):
    data = pd.DataFrame()
    for i in range(len(variables)):
        data[variables[i]]=np.array(parameters)[:,i]
    gyr = []
    std = []
    for i in range(len(gyr_list)):
        gyr.append(gyr_list[i].complex['mean'])
        std.append(gyr_list[i].complex['std'])
    data['gyr'] = gyr
    data['std'] = std
    return data

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

def _get_gyr(runs, fname='gyr.out', root=None, 
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
    gyr = pd.Series
    if root != None:
        root = '{}/'.format(root)
    for run in range(1, runs+1):
        fin = '{}{}/{}'.format(root, run, fname)
       
        try:

            temp = pd.read_csv(fin, 
                               delim_whitespace=True,
                               comment='#',
                               header=None)
            
            temp = temp.rename(columns={0:'ts', 
                                        1:'star',
                                        2:'siRNA',
                                        3:'gyr'})

            PMF_traj = '{}{}/{}'.format(root, run, f_traj)
            traj = pmf._get_traj(PMF_traj) # has columns ['ts', 'xi']

            # merge temp and traj
            merged = pd.merge(temp, traj, on='ts')[['xi', 'gyr']]
                     
            
            
            xi=xi.append(merged['xi'])
            gyr=gyr.append(merged['gyr'])
            
        except IOError:
            print "Warning: run {} did not work!".format(run)
    
    data = pd.DataFrame() # should end with long dataframe with
                          # columns ['xi', 'gyr']
    
    data['xi'] = xi
    data['gyr'] = gyr
    print 'data:'
    print data
    return data

def _get_complex_size(GYR, PMF, bound=1):
    xi = _get_min(PMF, bound=bound)
    print GYR.gyr
    data = GYR.gyr['xi'].isin(xi)
    print data
    mean = data['gyr'].mean()
    std = data['std'].std()
    return [mean, std]

def _write_gyr(dataframe, fname):
    dataframe.to_csv(fname, index=False)
    return

def _write_dg():
    return

def _plot_gyr(gyr, ax):
    gyr = gyr.gyr
    ax.errorbar(gyr['xi'], gyr['mean'], yerr=gyr['std'],
                capsize=2, fmt='kx', markersize=5, elinewidth=1)
    return

def _plot_gyr_list(GYR_LIST, ax):
    labels = GYR_LIST.labels
    for i in range(1, GYR_LIST.N+1):
        gyr = GYR_LIST.GYR
        ax.errorbar(gyr['ts'], gyr['{}_mean'.format(i)],
                    yerr=gyr['{}_std'.format(i)], label=labels[i],
                    capsize=2, fmt='{}'.format(markers[i-1]),
                    markersize=5, elinewidth=1)
    return

class GYR():
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
    def __init__(self, PMF, fname='gyr', runs=10, root=None, bound=1):
        if root != None:
            self.fname = '{}/{}.csv'.format(root, fname)
        else:
            self.fname = '{}.csv'.format(fname)
        self.gyr = _get_gyr(runs, root=root)
        self.complex = _get_complex_size(self, PMF, bound=bound)
    
    def write(self):
        _write_gyr(self.gyr, self.fname)
        return    

class GYR_LIST():
    def __init__(self, gyr_list, variables=[''], parameters=[[]],
                 units=None, fname='gyrs'):

        self.N = len(gyr_list)

        # self.dg should be a dataframe with
        # columns=['variables[0], ..., variables[N], dg, err]
        self.gyr = _collate(gyr_list, variables, parameters)

        # self.labels is given by e.g. {1: 'variable = parameter[0]',...
        # or {1: 'var[0] = param[0][0], var[1] = param[0][1]',...
        self.labels = _label_generator(variables, parameters)
       

    def plot(self):#, subax_dimensions):
        fig, ax = plt.subplots()
        _plot_gyr_list(self, ax)

        # make subax and plot dG
        #subax = _make_inset(ax, subax_dimensions)
        #_plot_dg(self, subax)

        plt.legend()

        plt.show()
