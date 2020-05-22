import numpy as np
import pandas as pd

from starpolymers.molecules import Salt
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
        assert self.charge == 0

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

    @property
    def charge(self):
        _charge = sum(self._atoms['q'].values)
        return _charge

    def add_molecule(self, molecule):
        #if not isinstance(molecule, AbstractMolecule):
        #    raise TypeError(
        #        'A molecule that is not an '
        #        'instance of abstract molecule was passed to system')
        def _atoms():
            _temp = molecule._atoms.copy()
            _temp['mol'] = self._molecules + 1
            self._atoms = pd.concat([self._atoms, _temp], sort=False).reset_index(drop=True)

        def _bonds():
            self._bonds = pd.concat([self._bonds, molecule._bonds], sort=False).reset_index(drop=True)
        def _angles():
            self._angles = pd.concat([self._angles, molecule._angles], sort=False).reset_index(drop=True)
        
        _atoms()
        _bonds()
        _angles()
        self._molecules += 1

    def neutralise(self, salt):

        def difference(_delta, n_anions, anion, n_cations, cation):
            diff = (n_anions * -anion + n_cations * cation) + _delta
            return diff

        delta = int(self.charge)

        if delta == 0:
            return

        if delta < 0:
            n_cations = abs( delta / abs(int(salt.cation)) )
            n_anions = 0
            
        if delta > 0:
            n_anions = abs( delta / abs(int(salt.anion)))
            n_cations = 0

        diff = difference(delta, n_anions, salt.anion, n_cations, salt.cation)
        i = 0
        while diff != 0:
            if diff > 0:
                n_anions += 1
            elif diff < 0:
                n_cations += 1
            else:
                raise ValueError('diff={} which is not allowed'.format(diff))
            diff = difference(delta, n_anions, salt.anion, n_cations, salt.cation)
            i+=1
            if i > 20:
                raise StopIteration
            
        n = int(n_cations + n_anions)
        self._atoms = pd.concat([self._atoms,
            pd.DataFrame({
                    'mol' : n * [self._molecules],
                    'type' : n * [salt._atoms['type'][0]],
                    'x' : np.random.uniform(size=(n, 1)).T[0],
                    'y' : np.random.uniform(size=(n, 1)).T[0],
                    'z' : np.random.uniform(size=(n, 1)).T[0],
                    'q' : [salt.cation] * int(n_cations) + [-salt.anion] * int(n_anions)
                })
            ], sort=True).reset_index(drop=True)
        
        self.assert_neutral()
        return