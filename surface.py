# surface.py

from starpolymers.packages import *
from starpolymers.dumper import DumpReader as dr
from multiprocessing import Pool
from scipy.spatial import distance
from mpl_toolkits.mplot3d import Axes3D
from itertools import permutations

def lj(r, eps, sig, rc):
    rr = 1/r ** 2
    sr6 = float(sigma)/r ** 6
    potential = 4 * eps * (sr6 - sr6*rr)
    return potential                 

class Molecule():
    def __init__(self, atoms):
        self.atoms = atoms

    def settings(self, dictionary):
        self.bound = dictionary['bound']
        self.rc = dictionary['rc']
        self.sigma = dictionary['sigma']
        self.epsilon = dictionary['epsilon']

    def surface(self, grid, pool_size=5):
        self.surface = grid.grid
        subgrids = np.vsplit(self.surface, 1)
        pool = Pool(pool_size)
        outputs = pd.Dataframe()
        print 'test'
        for i in range(len(atoms)):
            print 'atom number {} is running'.format(i)
            atoms[i].settings(self.bound,
                              self.epsilon,
                              self.sigma,
                              self.rc)
            output[i] = pool.map(atoms[i].lj, subgrids)
            self.surface = self.surface[output[i]]
        
class Grid():
    def __init__(self, resolution, fine_graining=10, box=50):
        self.box = box
        self.coarse = resolution
        self.fine = float(resolution)*10

    def initialise(self):
        n = self.box*2/self.coarse
        x = np.linspace(-self.box, self.box, int(n))
        self.grid = np.array(list(permutations(x,3)))
        print self.grid.shape

class Atom():
    def __init__(self, position):
        self.position = position

    def settings(bound, epsilon, sigma, rc):
        self.bound = bound
        self.eps = epsilon
        self.sigma = sigma
        self.cutoff = rc

    def lj(self, test_particle):
        r = distance.euclidean(self.position, test_particle)
        if r > self.cutoff:
            return False
        if lj(r, self.eps, self.sigma, self.cutoff) > self.bound:
            return True
        else:
            return False

class Converter():
    def __init__(self, path, ID, step, box=50, exclude=3):
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
        data = data[data['mol']!=self.exclude]
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
        settings={'bound':1.0,
                  'rc':6.0,
                  'epsilon':1.0,
                  'sigma':2.5}, exclude=3):
    molecule = Converter(path, ID, step, box=box,
                         exclude=exclude).create_molecule()
    molecule.settings(settings)
    grid = Grid(res)
    grid.initialise()
    return molecule.surface(grid)

def plot_surface(surface):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(surface[:,0], surface[:,1], surface[:,2],
               'o', s=1)
    plt.show()
    
