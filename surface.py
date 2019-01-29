# surface.py

from starpolymers.packages import *
from starpolymers.dumper import DumpReader as dr
from multiprocessing import Pool
from scipy.spatial import distance
from mpl_toolkits.mplot3d import Axes3D
from itertools import permutations

def lj(r, eps, sig, rc):
    rr = 1/r ** 2
    sr6 = float(sig)/r ** 6
    potential = 4 * eps * (sr6 - sr6*rr)
    return potential                 

class Molecule():
    def __init__(self, atoms):
        self.atoms = atoms # list

    def __call__(self):
        for atom in self.atoms:
            print atom.position

    def settings(self, dictionary):
        self.bound = dictionary['bound']
        self.rc = dictionary['rc']
        self.sigma = dictionary['sigma']
        self.epsilon = dictionary['epsilon']

    def initialise_space(self):
        atom_array = []
        for atom in self.atoms:
            atom_array.append(atom.position)
        atom_array = np.array(atom_array)
        length = np.max(distance.pdist(atom_array))
        return length

    def recentre(self):
        scaler = -self.atoms[len(self.atoms)/2].position
        for atom in self.atoms:
            atom.position = atom.position + scaler                

    def calculate_surface(self, grid):
        self.surface = grid.grid
        subgrids = np.vsplit(self.surface, len(self.surface))
        outputs = pd.DataFrame()
        for i in range(len(self.atoms)):
            output = []
            print 'atom number {} is running'.format(i)
            for j in subgrids:
                self.atoms[i].settings(self.bound,
                                       self.epsilon,
                                       self.sigma,
                                       self.rc)
                output.append(self.atoms[i].lj(j))
            outputs[i] = output
        self.surface = self.surface[outputs.any(axis=1)]
        
class Grid():
    def __init__(self, gridpoints, box=50):
        self.box = box
        self.coarse = gridpoints
        x = np.linspace(-self.box, self.box, self.coarse)
        self.grid = np.array(list(permutations(x,3)))
        

class Atom():
    def __init__(self, position):
        self.position = position

    def settings(self, bound, epsilon, sigma, rc):
        self.bound = bound
        self.eps = epsilon
        self.sigma = sigma
        self.cutoff = rc

    def lj(self, test_particle):
        r = distance.euclidean(self.position, test_particle)
        if r > self.cutoff:
            return False
        pot = lj(r, self.eps, self.sigma, self.cutoff) 
        if pot > self.bound[0] and pot < self.bound[1]:            
            return True
        else:
            return False

class Converter():
    def __init__(self, path, ID, step, box=50, exclude=[3]):
        self.path = path
        self.ID = ID
        self.step = step
        self.box = box
        self.exclude = exclude

    def create_molecule(self):
        # read dump file
        # ignore atoms that belong to exclude molecule

        reader = dr(self.ID, box=self.box)
        reader.change_path(self.path)
        data = reader.read(self.step)
        for i in self.exclude:
            data = data[data['mol']!=i]
        data = data[['x', 'y', 'z']].values

        # for each row
        # parse each row into a position

        atom_list = []

        for i in range(len(data)):
            atom = Atom(data[i])
            atom_list.append(atom)

        # now we should have a list of 1x3 arrays

        molecule = Molecule(atom_list)

        return molecule

def run(path, ID, step, box=50, res=0.1,
        settings={'bound':0.0,
                  'rc':2.5,
                  'epsilon':1.0,
                  'sigma':2.5}, exclude=3):
    molecule = Converter(path, ID, step, box=box,
                         exclude=exclude).create_molecule()
    molecule.settings(settings)
    # get longest distance between molecules
    # recentre on middle atom
    # set box size to this r_max
    grid = Grid(res)
    molecule.calculate_surface(grid)
    return molecule

def plot_surface(molecule):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(molecule.surface[:,0],
               molecule.surface[:,1],
               molecule.surface[:,2],
               'o', s=1)
    for atom in molecule.atoms:
        ax.scatter(atom.position[0],
                   atom.position[1],
                   atom.position[2],
                   'o',c='black', s=20,
                   alpha=0.4)
    plt.show()



    
