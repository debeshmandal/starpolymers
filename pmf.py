# pmf function

# aim is to read, write and plot PMF files and DeltaG
# together

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

markers = ['o','d','*','P','h','x','+']

def _scale(series, n=5):
    factor = series.tail(n).mean()
    array = series.values - factor    
    return array

def _collate(pmf_list):
    data = pd.DataFrame()
    for i in (range(len(pmf_list))):
        pmf = pmf_list[i].pmf
        if i == 0:
            data['xi'] = pmf['xi']
        data['{}_mean'.format(i+1)] = pmf['mean']
        data['{}_std'.format(i+1)] = pmf['std']
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
            label += units[j]
            if j+1 != len(variables):
                label += ', '
    dictionary[i+1] = label
    return dictionary

def _get_pmf(runs, fname='out.harmonic1.ti.pmf', root=None):
    """

    From a set of N runs starting at 1 and ending at N, return
    a DataFrame with columns=['xi', 'mean', 'std'], note that
    the PMF will be scaled by the last 5 values (i.e. towards the 
    Colvars targetCenters argument)

    Parameters:
    -----

    runs : int (positive)
           Number of different sampling runs

    fname : string
            name of PMF containing file

    
    Returns
    -----
     
    data : pd.DataFrame
           return DataFrame with columns=['xi', 'mean', 'std']

    """
    data = pd.DataFrame()
    pmfs = pd.DataFrame()
    xi = True
    if root != None:
        root = '{}/'.format(root)
    for run in range(1, runs+1):
        fin = '{}{}/{}'.format(root, run, fname)
        
        try:
            temp = pd.read_csv(fin, 
                               delim_whitespace=True,
                               comment='#',
                               header=None)

            pmfs[run] = temp.values[:,1]
            if xi:
                data['xi'] = temp.values[:,0]
                xi = False
        except IOError:
            print "Warning: run {} did not work!".format(run)
    
    data['mean'] = np.mean(pmfs.values, axis=1)
    data['mean'] = _scale(data['mean'])
    data['std'] = np.std(pmfs.values, axis=1)
    return data

def _write_pmf(dataframe, fname):
    dataframe.to_csv(fname, index=False)
    return

def _write_dg():
    return

def _plot_pmf(PMF, ax):
    pmf = PMF.pmf
    ax.errorbar(pmf['xi'], pmf['mean'], yerr=pmf['std'],
                capsize=2, fmt='kx', markersize=5, elinewidth=1)
    return

def _plot_pmf_list(PMF_LIST, ax):
    labels = PMF_LIST.labels
    for i in range(1, PMF_LIST.N+1):
        pmf = PMF_LIST[i].pmf
        ax.errorbar(pmf['xi'], pmf['{}_mean'.format(i)],
                    yerr=pmf['{}_std'.format(i)], label=labels[i],
                    capsize=2, fmt='{}'.format(markers[i-1]),
                    markersize=5, elinewidth=1)
    return

def _plot_dg(ax):
    return

class PMF():
    """

    Class containing pmf data in DataFrame form with method to write
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
    def __init__(self, fname='pmf', runs=10, root=None):
        if root != None:
            self.fname = '{}/{}.csv'.format(root, fname)
        else:
            self.fname = '{}.csv'.format(fname)
        self.pmf = _get_pmf(runs, root=root)
    
    def write(self):
        _write_pmf(self.pmf, self.fname)
        return

    def plot(self, fname='pmf', show=True):
        fig, ax = plt.subplots()
        _plot_pmf(self, ax)
        ax.set_xlabel(r'$\xi$ [$\sigma$]', fontsize='large')
        ax.set_ylabel(r'$\mathcal{W}$ [$k_BT$]', fontsize='large')
        ax.tick_params(labelsize='large')
        fig.savefig('{}.pdf'.format(fname))
        if show:
            plt.show()

class PMF_LIST():
    def __init__(self, pmf_list, variables=[''], parameters=[[]],
                 fname='pmfs'):

        self.N = len(pmf_list)

        # set self.PMF as a 2N+1-column DataFrame
        # ['xi', '1_mean', '1_std', ..., 'N_mean', N_std']
        self.PMF = _collate(pmf_list)

        # self.dg should be a dataframe with
        # columns=['variables[0], ..., variables[N], dg, err]
        self.dg = _get_dg(self.PMF, variables, parameters)

        # self.labels is given by e.g. {1: 'variable = parameter[0]',...
        # or {1: 'var[0] = param[0][0], var[1] = param[0][1]',...
        self.labels = _label_generator(variables, parameters)

    def plot(self):#, subax_dimensions):
        fig, ax = plt.subplots()
        _plot_pmf_list(self, ax)

        # make subax and plot dG
        #subax = _make_inset(ax, subax_dimensions)
        #_plot_dg(self, subax)

        plt.show()

class PMF_ARRAY():
    def __init__(self, series_list, fname='dg'):
        self.DG = _get_series(series_list)
        self.fname = fname

    def write(self):
        _write_dg(self.DG)
        return

    def plot_dg(self):
        _plot_dg(self, ax)
