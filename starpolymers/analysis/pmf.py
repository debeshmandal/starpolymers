# aim is to read, write and plot PMF files and DeltaG
# together

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from starpolymers.tools.design import base, cdict, markers

def _scale(series, n=5):
    factor = series.tail(n).mean()
    array = series.values - factor    
    return array

def _collate(pmf_list):
    data = pd.DataFrame()
    for i in (list(range(len(pmf_list)))):
        pmf = pmf_list[i].pmf
        if i == 0:
            data['xi'] = pmf['xi']
        data['{}_mean'.format(i+1)] = pmf['mean']
        data['{}_std'.format(i+1)] = pmf['std']
    return data

def _integrate(pmf):
    data = pmf.pmf.values
    dr = data[1, 0] - data[0, 0]
    result = data[:,0]**2 * data[:,1] * dr
    return np.sum(result)

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

def _get_traj(fname):
    try:
    
        data = pd.read_csv(fname, delim_whitespace=True,
                       header=None, comment='#')
        data = data.rename(columns={0: 'ts', 1: 'xi'})
    except IOError:
        print(('Warning: No traj for {}!'.format(fname)))
        return 'Error'
    return data

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
            print(("Warning: run {} did not work!".format(run)))
    
    data['mean'] = np.mean(pmfs.values, axis=1)
    data['mean'] = _scale(data['mean'])
    data['std'] = np.std(pmfs.values, axis=1)
    return data

def _dg(pmf, std):
    MAX = pmf.max()
    MIN = pmf.min()
    difference = MAX - MIN
    if pmf.idxmax() <= pmf.idxmin():
        sign = 1
    else:
        sign = -1
    value = sign * difference
    error = std[pmf.idxmax()] + std[pmf.idxmin()]
    return [value, error]

def _get_dg(PMF_LIST, variables, parameters):
    # self.dg should be a dataframe with
    # columns=['variables[0], ..., variables[N], dg, err]

    dataframe = pd.DataFrame()
    params = np.array(parameters)
    for i in range(len(variables)):
        dataframe[variables[i]] = params[:,i]
    
    dg = []
    std = []
    
    for i in range(PMF_LIST.N):
        temp = _dg(PMF_LIST.PMF['{}_mean'.format(i+1)],
                   PMF_LIST.PMF['{}_std'.format(i+1)])
        dg.append(temp[0])
        std.append(temp[1])

    dataframe['dg'] = dg
    dataframe['std'] = std
    
    return dataframe

def _get_integral(PMF_LIST, variables, parameters):

    dataframe = pd.DataFrame()
    params = np.array(parameters)
    for i in range(len(variables)):
        dataframe[variables[i]] = params[:,i]
 
    I = []
    for pmf in PMF_LIST:
        I.append(pmf.integral)
    dataframe['integral'] = I
    return dataframe

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

def _colour(PMF_LIST, order_by=0):
    colours = []
    params = PMF_LIST.dg.values[:, order_by]
    group, counts = np.unique(params, return_counts=True)
    #print group, counts
    #print cdict 
    if len(group)/len(params) == 1:
        return False
    for i in range(len(group)):
        c = base[i]
        for j in range(counts[i]):
            #print cdict[c][j]
            colours.append(cdict[c][j])
    #print colours
    return colours
           

def _plot_pmf_list(PMF_LIST, ax):
    labels = PMF_LIST.labels
    colours = _colour(PMF_LIST)
    if colours == False:
        ax.set_prop_cycle(color=plt.cm.viridis(np.linspace(0.1, 0.9,
                                               len(labels))))
    for i in range(1, PMF_LIST.N+1):
        pmf = PMF_LIST.PMF            
        if colours == False:
 
            ax.errorbar(pmf['xi'], pmf['{}_mean'.format(i)],
                        yerr=pmf['{}_std'.format(i)], 
                        label=labels[i], capsize=2, 
                        fmt='{}'.format(markers[i-1]),
                        markersize=5, elinewidth=1)

        else:
            
            ax.errorbar(pmf['xi'], pmf['{}_mean'.format(i)],
                        yerr=pmf['{}_std'.format(i)], 
                        label=labels[i], capsize=2, 
                        fmt='{}'.format(markers[i-1]),
                        mfc=colours[i-1], mec=colours[i-1],
                        ecolor=colours[i-1], markersize=5, 
                        elinewidth=1)

    ax.set_xlabel(r'$\xi$ [$\sigma$]')
    ax.set_ylabel(r'$\mathcal{W}$ [$k_BT$]')

    return

def _plot_dg(self, ax, var=0, x_axis=1):
    params = self.dg.values[:,var]
    group, counts = np.unique(params, return_index=True)
    if len(group) == len(params):
        ax.errorbar(self.dg.values[:,var], self.dg.values[:,-2],
                    yerr=self.dg.values[:,-1], fmt='kx:')
    else:
        for i in range(len(group)):
            temp = pd.DataFrame(self.dg.values[:, [var, x_axis,-2,-1]])
            temp = temp.rename(columns={0: 'group',
                                        1: 'X',
                                        2: 'Y',
                                        3: 'err'})
            data = temp[temp['group']==group[i]]
            X = data['X'].values
            Y = data['Y'].values
            err  = data['err'].values
            lab = group[i]
            colour = base[i]
            marker = markers[i]
            ax.errorbar(X, Y, yerr=err, label=lab, 
                        capsize=2, fmt='{}:'.format(marker),
                        mfc=colour, mec=colour,
                        ecolor=colour, markersize=5,
                        elinewidth=1, color=colour)
    ax.set_xlabel(self.dg.columns[x_axis])
    ax.set_ylabel(r'$\Delta \mathcal{W}$ [$k_BT$]')
    ax.tick_params(direction='in')
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
        self.integral = _integrate(self)
    
    def write(self):
        _write_pmf(self.pmf, self.fname)
        return

    def plot(self, fout='pmf', show=True, legend_on=True,
             ax=None):
        if ax == None:
            fig, ax = plt.subplots()
        _plot_pmf(self, ax)
        ax.set_xlabel(r'$\xi$ [$\sigma$]', fontsize='large')
        ax.set_ylabel(r'$\mathcal{W}$ [$k_BT$]', fontsize='large')
        ax.tick_params(labelsize='large')
        if fout != None:
            fig.savefig('{}.pdf'.format(fout))
        if show:
            plt.show()

    def min(self, bound=1):
        """
        
        returns a tuple of (xi_min, xi_max) which
        provides the range where the free energy is a minimum

        """
        pmf = self.pmf
        pmf_min = pmf['mean'].min() - bound
        pmf_max = pmf['mean'].min() + bound

        data = pmf[pmf['mean']<pmf_max]
        data = data[data['mean']>pmf_min]
        
        xi_min = data['xi'].min()
        xi_max = data['xi'].max()

        return (xi_min, xi_max)
    

class PMF_LIST():
    def __init__(self, pmf_list, variables=[''], parameters=[[]],
                 units=None, fname='pmfs'):

        self.N = len(pmf_list)

        # set self.PMF as a 2N+1-column DataFrame
        # ['xi', '1_mean', '1_std', ..., 'N_mean', N_std']
        self.PMF = _collate(pmf_list)

        # self.dg should be a dataframe with
        # columns=['variables[0], ..., variables[N], dg, err]
        self.dg = _get_dg(self, variables, parameters)

        self.integral = _get_integral(pmf_list, variables, parameters)

        # self.labels is given by e.g. {1: 'variable = parameter[0]',...
        # or {1: 'var[0] = param[0][0], var[1] = param[0][1]',...
        self.labels = _label_generator(variables, parameters)
       

    def plot(self, fout='pmf.pdf', legend_cols=2, 
             sub=[0.25, 0.675, 0.2, 0.2], sub_var=0, 
             sub_x=1, ax=None, legend_on=True, 
             show=False):
            # subax_dimensions):
        if ax == None:
            fig, ax = plt.subplots()
        _plot_pmf_list(self, ax)
        if legend_on==True:
            plt.legend(frameon=False, ncol=legend_cols)

        # make subax and plot dG
        subax = plt.axes(sub)
        _plot_dg(self, subax, var=sub_var, x_axis=sub_x)
        if fout != None:
            plt.savefig(fout)
        if show:
            plt.show()

    def write(self, fout='pmfs.csv'):
        self.PMF.to_csv(fout, index=False)
