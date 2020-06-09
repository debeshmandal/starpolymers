# Symmetry analysis using PCA

import pandas as pd
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from starpolymers.readers.dumper import DumpReader as dr

def product(vector):
    result = 1
    for i in range(len(vector)):
        result = result * vector[i]
    return result

def read_data(fname, exclude=3):
        """

        Reads data from LAMMPS dump file

        """
        raw_data = pd.read_csv(fname, delimiter=' ', skiprows=9, header=None)
        positions = raw_data[raw_data[9]!=exclude].iloc[:,2:5].values
        return positions

def scale(positions, box=50):
    
    # pick one particle

    ID = len(positions)/2

    # centre this at (0,0,0)

    scale_by = -positions[ID,:]

    # for all particles
    # scale by same amount

    scaled = positions+scale_by

    # measure distance in x,y,z direction

    for i in range(len(scaled)):
        for j in range(3):
            scaled[i,j] = scaled[i,j] - box*round(scaled[i,j]/box)

    return scaled

class SymmetryAnalyser:
    def __init__(self, filename, exclude=3):
        self.f = filename
        self.positions = read_data(filename, exclude)

    def scale_complex(self):
        self.positions = scale(self.positions)

    def plot_complex(self, eig=False):
        """

        Plots complexes onto 3D axis

        """

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        positions = self.positions
        ax.scatter(positions[:,0], positions[:,1], positions[:,2], label='{}'.format(self.f), alpha = 0.1)
        if eig==True:
            A = self.statistics()
            C = A[0]
            basis = np.mat(np.zeros([3,3]))
            for i in range(3):
                basis[i] = C + (np.array(A[2])[i] * np.array(A[3])[:,i])
                ax.plot([C[0],basis[i,0]],[C[1],basis[i,1]],[C[2],basis[i,2]],c='black', alpha=0.5)
        plt.show()

    def statistics(self, pretty=False, symmetry_only=False):
        positions = self.positions
        mean = np.array([np.mean(positions[:,0]),
                         np.mean(positions[:,1]),
                         np.mean(positions[:,2])])

        covariance = np.cov(positions, rowvar = False)
        eigenvalues = la.eig(covariance)[0]
        eigenvectors = la.eig(covariance)[1]
        scaled_eigenvalues = eigenvalues / eigenvalues.max()
        symmetry_index = product(scaled_eigenvalues)
        if pretty == True:
            print(('\nThe mean (centre-of-mass) is \n', mean))
            print(('\nThe covariance matrix is \n', covariance))
            print(('\nThe eigenvalues are \n', eigenvalues))
            print(('\nThe eigenvectors are \n', eigenvectors))
            print(('\nThe scaled eigenvalues are \n', scaled_eigenvalues))
            print(('\nThe symmetry_index is {0:.4e}'.format(symmetry_index)))
        if symmetry_only == True:
            return symmetry_index
        else:
            return [mean, covariance, eigenvalues, eigenvectors, scaled_eigenvalues, symmetry_index]

    def transform(self, plot=False):
        positions = np.mat(self.positions)
        transformer = np.mat(self.statistics()[3])
        new = positions * transformer
        if plot == True:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(new[:,0],new[:,1],new[:,2])
            plt.show()
        return new
        

        
