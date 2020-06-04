import numpy as np
import scipy.spatial.distance as D
import pandas as pd
from math import pi
import matplotlib.pyplot as plt

def _volume(bins):
    shell_width = bins[1] - bins[0]
    return 4 * pi * bins ** 2 * shell_width

def _get_distance_vector(first, second, rmin, rmax):
    array = D.cdist(first, second)
    vector = np.ravel(array)
    vector = vector[vector>=rmin]
    vector = vector[vector<rmax]
    return vector

def _bins(rmin, rmax, binwidth):
    return np.arange(rmin, rmax, binwidth)

def _get_histogram(vector, bins):
    return np.histogram(vector, bins=bins)

def _get_rdf(histogram_object, bulk_density):
    bins = histogram_object[1] + 0.5
    bins = bins[:-1]
    values = histogram_object[0] / _volume(bins)
    data = pd.DataFrame()
    data['r'] = bins
    data['g(r)'] = values / bulk_density
    return data

def _get_condensed_radius(df):
    return df['r'][df['g(r)'].idxmax]

def _plot_rdf(rdf, ax):
    ax.plot(rdf.values[:,0], rdf.values[:,1],'-')
    ax.set_xlabel('r', fontsize='large')
    ax.set_ylabel('g(r)', fontsize='large')

def _get_average(RDF_LIST):
    counter = 0
    dataframe = pd.DataFrame()
    temp = pd.DataFrame()
    for rdf in RDF_LIST:
        counter +=1
        if counter == 1:
            dataframe['r'] = rdf.rdf['r']
        temp[counter] = rdf.rdf['g(r)'].values
    
    dataframe['g(r)'] = np.mean(temp.values, axis=1)

    return dataframe

def _get_NC(first, second, cutoff):
    count = 0
    for atom in second:
        distances = D.cdist(first, atom.reshape(1,3))
        if not len(distances[distances<=cutoff]):
            continue
        else:
            count+=1
    return count

class RDF():
    def __init__(self, atoms1, atoms2, bin_width=1.0, rmin=0.0, rmax=20.0, L=70, condensed_radius=None):
        self.distances = _get_distance_vector(atoms1, atoms2, rmin, rmax)
        self.histogram = _get_histogram(self.distances, _bins(rmin, rmax, bin_width))
        self.density = float(len(atoms2)) / L ** 3
        self.rdf = _get_rdf(self.histogram, self.density)
        if condensed_radius != None:
            self.condensed_radius = condensed_radius
        else:
            self.condensed_radius = _get_condensed_radius(self.rdf)
        self.NC = _get_NC(atoms1, atoms2, self.condensed_radius)
        
    
    def plot(self, show=True):
        fig, ax = plt.subplots()
        _plot_rdf(self.rdf, ax)
        if show:
            plt.show()

class RDF_AVERAGE():
    def __init__(self, RDF_LIST):
        self.rdf = _get_average(RDF_LIST)

    def plot(self, show=True):
        fig, ax = plt.subplots()
        _plot_rdf(self.rdf, ax)
        if show:
            plt.show()

def _collate(RDF_LIST):
    counter = 1
    results = pd.DataFrame()
    for i in RDF_LIST:
        if counter == 1:
            results['r'] = i.rdf['r']
        results['{}_g(r)'.format(counter)] = i.rdf['g(r)']
        counter+=1
    return results

def _plot_list(RDF_LIST_OBJ, ax):
    for i in range(1, RDF_LIST_OBJ.N+1):
        ax.plot(RDF_LIST_OBJ.rdf['r'], RDF_LIST_OBJ.rdf['{}_g(r)'.format(i)], '-', label=i)

class RDF_LIST():
    def __init__(self, RDF_LIST):
        self.N = len(RDF_LIST)
        self.rdf = _collate(RDF_LIST)
    
    def plot(self, show=True):
        fig, ax = plt.subplots()
        _plot_list(self, ax)
        ax.set_xlabel('r', fontsize='large')
        ax.set_ylabel('g(r)',fontsize='large')
        ax.legend()
        if show:
            plt.show()