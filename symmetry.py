# Symmetry analysis using PCA

import pandas as pd
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def product(vector):
    result = 1
    for i in range(len(vector)):
        result = result * vector[i]
    return result

def read_data(fname):
        """

        Reads data from LAMMPS dump file

        """
        raw_data = pd.read_csv(fname, delimiter=' ', skiprows=9, header=None)
        positions = raw_data[raw_data[9]!=3].iloc[:,2:5].values
        return positions

class SymmetryAnalyser:
    def __init__(self, filename):
        self.f = filename
        self.positions = read_data(filename)    

    def plot_complex(self):
        """

        Plots complexes onto 3D axis

        """

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        positions = self.positions
        ax.scatter(positions[:,0], positions[:,1], positions[:,2], label='{}'.format(self.f))
        plt.show()

    def statistics(self):
        positions = self.positions
        mean = np.array([np.mean(positions[:,0]),
                         np.mean(positions[:,1]),
                         np.mean(positions[:,2])])

        covariance = np.cov(positions, rowvar = False)
        eigenvalues = la.eig(covariance)[0]
        eigenvectors = la.eig(covariance)[1]
        scaled_eigenvalues = eigenvalues / eigenvalues.max()
        symmetry_index = product(scaled_eigenvalues)
        return [mean, covariance, eigenvalues, eigenvectors, scaled_eigenvalues, symmetry_index]
