import numpy as np
import pandas as pd

from starpolymers.molecules._common import registry, AbstractMolecule

class System():
    def __init__(self, box, molecules=[], atom_masses=[1.0], bond_types=1, angle_types=1):

        self._molecules = int(0)
        self._box = float(box)
        self._masses = atom_masses
        self.types = {
            'atom' : len(atom_masses),
            'bond' : bond_types,
            'angle' : angle_types
        }
        self._atoms = pd.DataFrame(
            columns = ['mol'] + registry.columns['atoms']
        )
        self._bonds = pd.DataFrame(
            columns = registry.columns['atoms']
        )
        self._angles = pd.DataFrame(
            columns = registry.columns['angles']
        )
        self.assert_neutral()

        for molecule in molecules:
            self.add_molecule(molecule)

    def assert_neutral(self):
        assert True

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles
    
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
    def n(self):
        return {
            'atoms' : len(self.atoms),
            'bonds' : len(self.bonds),
            'angles' : len(self.angles),
            'molecules' : self._molecules
        }

    def add_molecule(self, molecule):
        #if not isinstance(molecule, AbstractMolecule):
        #    raise TypeError(
        #        'A molecule that is not an '
        #        'instance of abstract molecule was passed to system')
        def _atoms():
            _temp = molecule._atoms.copy()
            _temp['mol'] = self._molecules + 1
            self._atoms = pd.concat([self._atoms, _temp], sort=False)

        def _bonds():
            self._bonds = pd.concat([self._bonds, molecule._bonds], sort=False)
        def _angles():
            self._angles = pd.concat([self._angles, molecule._angles], sort=False)
        
        _atoms()
        _bonds()
        _angles()
        self._molecules += 1