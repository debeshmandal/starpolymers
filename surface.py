from starpolymers.dumper import DumpReader as dr
from starpolymers.tools.cylinders import Molecule
from starpolymers.packages import *

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

def _get_min(PMF, bound=5):
    return PMF.min(bound=bound)

def _surf():
    results = []
    results_nosalt = []
    for i in range(timesteps):
                results.append(Surface().potential)
                results_nosalt.append(Surface().potential_nosalt)
    data = pd.DataFrame()
    data['ts'] = ts
    data['pot'] = results
    data['pot_nosalt'] = results_nosalt
    return

def _get_surf(runs, timesteps, root=None, 
              f_traj='out.colvars.traj', radius=2.0, 
              bjerrum=2, T=1.2, dna=2):
    xi = pd.Series()
    surf = pd.Series()
    surf_nosalt = pd.Series()
    if root != None:
        root = '{}/'.format(root)
           
    for run in range(1, runs+1):
          
        fname = 'dump.{}.lammpstrj'
        fin = '{}{}/{}'.format(root, run, fname)
       
        try:
            os.system('cd {}/{} && tar xzf dumps.tar.gz'.format(root, run))

            _surf()
            
            os.system('cd {}/{} && rm *.lammpstrj'.format(root, run))
           
            PMF_traj = '{}{}/{}'.format(root, run, f_traj)
            try:
                traj = pmf._get_traj(PMF_traj) # has columns ['ts', 'xi']
            except:
                None
            # merge temp and traj
            merged = pd.merge(temp, traj, on='ts')         
            
            xi=xi.append(merged['xi'])
            surf=surf.append(merged['pot'])
            surf_nosalt = surf_nosalt.append(merged['pot_nosalt'])
            
        except IOError:
            print "Warning: run {} did not work!".format(run)
    
    data = pd.DataFrame() # should end with long dataframe with
                          # columns ['xi', 'pot']
    
    data['xi'] = xi
    data['pot']=surf
    data['pot_nosalt'] = surf_nosalt
    data = data.sort_values(by='xi')
    data = data[data['xi']<75]
    data = data.reset_index(drop=True)
   
    return data

def _get_dna(dataframe, mol):
    array = dataframe[dataframe['mol']==mol][['x','y','z']].values
    return array

def _potential(surface, atoms, bjerrum, T):
    inverse = atoms[:,-1].transpose() * (1./cdist(surface, atoms[:,:-1]))
    return pd.DataFrame(inverse).sum(axis=1)
    
def _get_potential(Surface, bjerrum, T, salt=True):
    gridpoints = Surface.surface
    if salt:
        potential = _potential(Surface.surface, Surface.atoms, bjerrum, T)
    else:
        potential = _potential(Surface.surface, Surface.atoms_nosalt, bjerrum, T)
    results = pd.DataFrame()
    results['x'] = gridpoints[:,0]
    results['y'] = gridpoints[:,1]
    results['z'] = gridpoints[:,2]
    results['pot'] = potential
    return results

def _plot_complex(SURF_LIST, ax, var=0, x_axis=1):
    labels = SURF_LIST.labels
    #colours = _colour(SURF_LIST)
    data = SURF_LIST.COMPLEX.values
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
            lab = '{}={}'.format(SURF_LIST.COMPLEX.columns[var], group[i])
            colour = base[i]
            marker = markers[i]
            ax.errorbar(X, Y, yerr=err, label=lab,
                        capsize=2, fmt='{}:'.format(marker),
                        mfc=colour, mec=colour,
                        ecolor=colour, markersize=5,
                        elinewidth=1, color=colour)
    
    ax.set_xlabel(SURF_LIST.COMPLEX.columns[x_axis], fontsize='large')
    ax.set_ylabel(r'$\langle R_g \rangle _{complex}$ [$\sigma$]',
                  fontsize='large')
    ax.tick_params(direction='in', labelsize='large')
    return

def _get_size(SURF, PMF, bound=5, mol='pot'):
    xi = _get_min(PMF, bound=bound) 
    data = SURF.surf
    data = data[data['xi'] > xi[0]]
    data = data[data['xi'] < xi[1]+1]
    mean = data[mol].mean()
    std = data[mol].std()
    return [mean, std]

    
class Surface():
    def __init__(self, fname, radius=2.2, bjerrum=2, T=1.2, dna=2):
        self.data = dr(fname).read(kind='positions-long')
        
        dna = _get_dna(self.data, 2)
        self.atoms = self.data[['x','y','z','q']].values
        self.surface = Molecule(dna, N=125, sigma=radius).points
        self.potential = _get_potential(self, bjerrum, T)

        self.data_nosalt = self.data[self.data['mol']!=3]
        self.atoms_nosalt = self.data_nosalt[['x','y','z','q']].values
        self.potential_nosalt = _get_potential(self, bjerrum, T, salt=False)

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.potential['x'].values,
                self.potential['y'].values,
                self.potential['z'].values,
                c=self.potential['pot'].values)
        ax.set_xlim(-25,25)
        ax.set_ylim(-25,25)
        ax.set_zlim(-25,25)                
        plt.show()

class SURF():
    def __init__(self, PMF, fname='surf', runs=10, root=None):
        if root != None:
            self.fname = '{}/{}.csv'.format(root, fname)
        else:
            self.fname = '{}.csv'.format(fname)
        self.surf = _get_surf()
        self.complex = {'mean':_get_size(self, PMF, bound=bound)[0],
                        'std':_get_size(self, PMF, bound=bound)[1]}
        self.complex_nosalt = {'mean':_get_size(self, PMF, bound=bound, mol='pot_nosalt')[0],
                               'std':_get_size(self, PMF, bound=bound, mol='pot_nosalt')[1]}

    def write(self):
        self.surf.to_csv(self.fname, index=False)
        return None

class SURF_LIST():
    def __init__(self, surf_list, variables=[''], parameters=[[]],
                 units=None, fname='surfs'):

        self.N = len(surf_list)
        self.SURF = _collate(surf_list)
        self.COMPLEX = _collate_complex(surf_list, variables, parameters)
        self.labels = _label_generator(variables, parameters)

    def plot_complex(self, fout='surf.pdf', legend_cols=2,
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