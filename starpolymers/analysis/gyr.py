"""
Objects to calculate the radius of gyration, using
information from a potential of mean force (PMF) generated
by Colvars
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from . import pmf
from .pmf import PMF
from starpolymers.tools.design import base, cdict, markers

def _collate(gyr_list):
    data = dict()
    for i in (list(range(len(gyr_list)))):
        gyr = gyr_list[i].gyr        
        data[i+1] = gyr
    return data

def _colour(GYR_LIST, order_by=0):
    colours = []
    params = GYR_LIST.COMPLEX.values[:, order_by]
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

def _get_min(PMF, bound=5):
    return PMF.min(bound=bound) # returns range of xi as tuple

def _collate_complex(gyr_list, variables, parameters):
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
    gyr = pd.Series()
    star = pd.Series()
    siRNA = pd.Series()
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
            try:
                traj = pmf._get_traj(PMF_traj) # has columns ['ts', 'xi']
            except:
                None
            # merge temp and traj
            merged = pd.merge(temp, traj, on='ts')
                     
            
            
            xi=xi.append(merged['xi'])
            star = star.append(merged['star'])
            siRNA = siRNA.append(merged['siRNA'])
            gyr=gyr.append(merged['gyr'])
            
        except IOError:
            print(("Warning: run {} did not work!".format(run)))
    
    data = pd.DataFrame() # should end with long dataframe with
                          # columns ['xi', 'gyr']
    
    data['xi'] = xi
    data['gyr'] = gyr
    data['siRNA'] = siRNA
    data['star'] = star
    data = data.sort_values(by='xi')
    data = data[data['xi']<75]
    data = data.reset_index(drop=True)
    
    return data

def _get_size(GYR, PMF, bound=5, mol='gyr'):
    xi = _get_min(PMF, bound=bound) 
    data = GYR.gyr
    data = data[data['xi'] > xi[0]]
    data = data[data['xi'] < xi[1]+1]
    mean = data[mol].mean()
    std = data[mol].std()
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

def _plot_xi(GYR_LIST, ax, mol='star'):
    labels = GYR_LIST.labels
    colours = _colour(GYR_LIST)
    if colours == False:
        ax.set_prop_cycle(color=plt.cm.viridis(np.linspace(0.1, 0.9,
                                                           len(labels))))

    
    
    for i in range(1, GYR_LIST.N+1):
        
        gyr = GYR_LIST.GYR[i]
        
        data = pd.DataFrame()
        data['xi'] = gyr['xi']
        data['rg'] = gyr['{}'.format(mol)]
        
        win = len(data['xi'])/100
        roller = data.rolling(win).mean()
        xi = roller['xi']
        rg = roller['rg']
        
        if colours == False:
           ax.plot(xi, rg, label=labels[i])
                    
        else:
           ax.plot(xi, rg, color=colours[i-1],
                   label=labels[i], marker=markers[i],
                   markevery=500)
                   
    ax.set_xlabel(r'$\xi$ [$\sigma$]')
    ax.set_ylabel(r'$\langle R_g \rangle$({}) [$\sigma$]'.format(mol))
    return

def _plot_complex(GYR_LIST, ax, var=0, x_axis=1):
    labels = GYR_LIST.labels
    #colours = _colour(GYR_LIST)
    data = GYR_LIST.COMPLEX.values
    params = data[:, var]
    group, counts = np.unique(params, return_index=True)
    if len(group) == len(params):

        ax.errorbar(data[:,0], data[:,1], yerr=data[:,2],
                    capsize=2, fmt='kx', markersize=5, elinewidth=1)

    else:

        for i in range(len(group)):
            temp = pd.DataFrame(data[:, [var, x_axis, -2, -1]])
            temp = temp.rename(columns={0:'group',
                                        1:'X',
                                        2:'Y',
                                        3:'err'})
            temp = temp[temp['group']==group[i]]
            X = temp['X'].values
            Y = temp['Y'].values
            err = temp['err'].values
            lab = '{}={}'.format(GYR_LIST.COMPLEX.columns[var], group[i])
            colour = base[i]
            marker = markers[i]
            ax.errorbar(X, Y, yerr=err, label=lab,
                        capsize=2, fmt='{}:'.format(marker),
                        mfc=colour, mec=colour,
                        ecolor=colour, markersize=5,
                        elinewidth=1, color=colour)
    
    ax.set_xlabel(GYR_LIST.COMPLEX.columns[x_axis], fontsize='large')
    ax.set_ylabel(r'$\langle R_g \rangle _{complex}$ [$\sigma$]',
                  fontsize='large')
    ax.tick_params(direction='in', labelsize='large')
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
    def __init__(self, PMF, fname='gyr', runs=10, root=None, bound=5):
        if root != None:
            self.fname = '{}/{}.csv'.format(root, fname)
        else:
            self.fname = '{}.csv'.format(fname)
        self.gyr = _get_gyr(runs, root=root)
        self.complex = {'mean':_get_size(self, PMF, bound=bound)[0],
                        'std':_get_size(self, PMF, bound=bound)[1]}
    
    def write(self):
        _write_gyr(self.gyr, self.fname)
        return    

class GYR_LIST():
    def __init__(self, gyr_list, variables=[''], parameters=[[]],
                 units=None, fname='gyrs'):

        self.N = len(gyr_list)
        self.GYR = _collate(gyr_list)
        

        # self.dg should be a dataframe with
        # columns=['variables[0], ..., variables[N], dg, err]
        self.COMPLEX = _collate_complex(gyr_list, variables, parameters)

        # self.labels is given by e.g. {1: 'variable = parameter[0]',...
        # or {1: 'var[0] = param[0][0], var[1] = param[0][1]',...
        self.labels = _label_generator(variables, parameters)

    def plot_xi(self, fout='gyr.pdf', mol='star', ax=None, 
                legend_on=True, show=False):
        if ax == None:
            fig, ax = plt.subplots()
        _plot_xi(self, ax, mol=mol)

        # make subax and plot dG
        #subax = _make_inset(ax, subax_dimensions)
        #_plot_dg(self, subax)
        
        if legend_on:
            plt.legend(frameon=False, fontsize='large')
        if fout != None:
            plt.savefig(fout)
        if show:
            plt.show()

    def plot_complex(self, fout='gyr.pdf', legend_cols=2,
                     var=0, x_axis=1, ax=None, legend_on=True,
                     show=False):
        if ax == None:
            fig, ax = plt.subplots()
        _plot_complex(self, ax, var=var, x_axis=x_axis)
        if legend_on:
            plt.legend(frameon=False, fontsize='large')
        if fout != None:
            plt.savefig(fout)
        if show:
            plt.show()
