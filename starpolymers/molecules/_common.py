import pandas as pd
from starpolymers.tools import geometry

class _MoleculeRegistry():
    def __init__(self):
        self.lam_list = ['star', 'DNA']
        self.salt_list = ['salt']
        self.molecule_list = self.lam_list + self.salt_list
        self.columns = {
            'atoms' : ['type', 'q', 'x', 'y', 'z'],
            'bonds' : ['type', 'atom_1', 'atom_2'],
            'angles' : ['type', 'atom_1', 'atom_2', 'atom_3']
        }
        self.spacing = 1.
        self.directions = geometry.direction
        self.start = geometry.translation

    @property
    def settings(self):
        return {
            'spacing' : self.spacing,
            'direction' : self.directions
        }

registry = _MoleculeRegistry()

class AbstractMolecule:
    def __init__(self, item, settings=registry.settings, mol=1, **kwargs):
        
        self._item = item
        self._atoms = pd.DataFrame(
            columns = registry.columns['atoms']
        )
        self._bonds = pd.DataFrame(
            columns = registry.columns['bonds']
        )
        self._angles = pd.DataFrame(
            columns = registry.columns['angles']
        )
        self.settings = settings
        self.mol = mol
        self._kwargs = kwargs

    @property
    def kind(self):
        return self.item['molecule']

    @property
    def charge(self):
        _charge = sum(self._atoms['q'].values)
        return _charge

    @property
    def atoms(self):
        # import self._atoms as data
        data = self._atoms
    
        # check format of dataframe
        assert data.columns == registry.columns['atoms']

        # iterate over all lines
        lines = []
        for i in range(len(data)):
            # get values
            atom_type = data.loc[i]['type']
            charge = data.loc[i]['q']
            x = data.loc[i]['x']
            y = data.loc[i]['y']
            z = data.loc[i]['z']
            
            # create line
            line = '{} {} {} {} {} {} {}\n'.format(
                self.kind,
                atom_type, 
                charge,
                x, 
                y, 
                z
            )
            lines.append(line)
        
        # export string
        string = ''.join(lines)
        return string

    @property
    def n(self):
        return {
            'atoms' : len(self._atoms),
            'bonds' : len(self._bonds),
            'angles' : len(self._angles)
        }

class ParticleArray(AbstractMolecule):
    def __init__(self, item, **kwargs):
        AbstractMolecule.__init__(self, item, **kwargs)
        self.types = {'atom' : item.get('atom_type', 1)}

class Molecule(AbstractMolecule):
    def __init__(self, item, **kwargs):
        AbstractMolecule.__init__(self, item, **kwargs)
        self.types = {
            'atom' : item.get('atom_type', 1),
            'bond' : item.get('bond_type', 1),
            'angle' : item.get('angle_type', 1)
        }
