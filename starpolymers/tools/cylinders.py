# Convert set of 3D points to a cylindrical mesh
# to represent the surface

import numpy as np
import matplotlib.pyplot as plt
from math import pi
import pandas as pd

identity_matrix = np.array([[1,0,0],[0,1,0],[0,0,1]])

def cylindrical2cartesian(array):
    new = pd.DataFrame()
    r = array[:,0]
    new['x'] = r * np.cos(array[:,1])
    new['y'] = r * np.sin(array[:,1])
    new['z'] = array[:,2]
    return new.values

def init_cylinder(radius, N):
    height = np.linspace(0., 1., num=int(np.sqrt(N)))
    angle = np.linspace(0., 2.*pi, num=int(np.sqrt(N)))
    array = np.ndarray((len(height) * len(angle), 3))
    count = 0
    for i in range(len(height)):
        for j in range(len(angle)):
            array[count,0]=radius
            array[count,1]=angle[j]
            array[count,2]=height[i]
            count+=1
    points = cylindrical2cartesian(array)
    return points

def unit(vector):
    return vector/np.linalg.norm(vector)

def rotate(vector_1, vector_2):
    vec_1=unit(vector_1)
    vec_2=unit(vector_2)
    cross = np.cross(vec_1, vec_2)
    dot = np.dot(vec_1, vec_2)
    skew_matrix = np.array([[0        , -cross[2], cross[1]],
                            [cross[2] , 0        , -cross[0]],
                            [-cross[1], cross[0] , 0       ]])
    
    rotation = identity_matrix + skew_matrix + (1./(1.+dot))*np.dot(skew_matrix, skew_matrix)
    return rotation   


def transform(init, atom_1, atom_2):

    # rotate
    target = atom_2 - atom_1
    unit_vector = np.array([0,0,1])
    R = rotate(unit_vector, target)

    # translate
    T = atom_1
    
    # scale length
    L = np.array([[1,0,0,0],
                  [0,1,0,0],
                  [0,0,np.linalg.norm(target),0],
                  [0,0,0,1]])

    transformation_matrix = np.array([[R[0,0], R[0,1],   R[0,2]  ,  T[0]],
                                      [R[1,0], R[1,1],   R[1,2]  ,  T[1]],
                                      [R[2,0], R[2,1],   R[2,2],    T[2]],
                                      [0     , 0     ,       0   ,   1  ]])

    transformation_matrix = np.matmul(transformation_matrix, L)

    temp_points = np.ones((init.shape[0], init.shape[1]+1))
    temp_points[:,:-1] = init
    points = np.matmul(transformation_matrix, temp_points.transpose())
    points = points.transpose()[:,:-1]
    
    return points
    

def cylinder(atom_1, atom_2, radius, N):
    init = init_cylinder(radius, N)
    final = transform(init, atom_1, atom_2)
    return final


class Cylinder():

    def __init__(self, pair, radius, N):
        """
        Constructs a cylinder between a pair of atoms

        Input:

        pair : np.array([[x1,y1,z1],[x2,y2,z2]]) - positions of atoms
        radius : float - radius of cylinder
        N : int - number of points overall

        Attributes:

        points : np.ndarray((N, 3)) evenly distributed grid of points
        """
        atom_1 = pair[0]
        atom_2 = pair[1]
        self.points = cylinder(atom_1, atom_2, radius, N)


class Molecule():
    def __init__(self, array, sigma=1.1, N=100):
        self.atoms = array
        self.radius = 2. * sigma
        points = np.array([[0,0,0]])
        for i in range(len(array)-1):
            pair = array[i:i+2, :]
            grid = Cylinder(pair, self.radius, N).points
            points = np.concatenate((points, grid))
        self.points = points[1:,:]

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.points[:, 0],
                self.points[:, 1],
                self.points[:, 2], 'bo', markersize=1)
        ax.plot(self.atoms[:, 0],
                self.atoms[:, 1],
                self.atoms[:, 2], 'ko-', markersize=10)
        plt.show()