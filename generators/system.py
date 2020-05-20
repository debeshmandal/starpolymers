import numpy as np

class System():
    def __init__(self, box, molecules=[], atom_masses=[1.0], bond_types=1, angle_types=1):
        self._molecules = molecules
        self._box = box
        self._masses = atom_masses
        self.types = {
            'atom' : len(atom_masses),
            'bond' : bond_types,
            'angles' : angle_types
        }
        self.assert_neutral()

    def assert_neutral(self):
        assert True
    
    @property
    def box(self):
        return {
            'xlo' : -self._box / 2,
            'xhi' : self._box / 2,
            'ylo' : -self._box / 2,
            'yhi' : self._box / 2,
            'zlo' : -self._box / 2,
            'zhi' : self._box / 2
        }

    @property
    def atom_types(self):
        return self.types['atom']

    @property
    def bond_types(self):
        return self.types['bond']

    @property
    def angle_types(self):
        return self.types['angle']

    @property
    def n(self):
        return {
            'atoms' : len(self.atoms),
            'bonds' : len(self.bonds),
            'angles' : len(self.angles)
        }

    def add_molecule(self, molecule, mol=1):
        return