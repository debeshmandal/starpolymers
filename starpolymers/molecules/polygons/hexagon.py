#!/usr/bin/env python
"""hexagon.py - auto-generated by softnanotools"""
from softnanotools.logger import Logger
logger = Logger(__name__)
import random

import numpy as np
import pandas as pd

from .._common import Molecule, registry

class Hexagon(Molecule):
    """
    Linear Polyelectrolyte with given length
    """
    def __init__(self, item, **kwargs):
        Molecule.__init__(self, item, **kwargs)
        self._atoms = self.generate_atoms()
        self._bonds = self.generate_bonds()
        self._angles = self.generate_angles()

    def generate_atoms(self):
        """
        Returns dataframe with the following columns
        """
        r3o2 = np.sqrt(3.0) / 2
        positions = np.array([
            [1.0, 0, 0],
            [0.5, r3o2, 0],
            [-0.5, r3o2, 0],
            [-1.0, 0, 0],
            [-0.5, -r3o2, 0],
            [0.5, -r3o2, 0]
        ])

        data = pd.DataFrame({
            'mol' : len(positions) * [self.mol],
            'type' : len(positions) * [self.types['atom']],
            'x' : positions[:, 0],
            'y' : positions[:, 1],
            'z' : positions[:, 2],
            'q' : self.generate_charges(positions)
        })

        return data

    def generate_bonds(self):

        data = pd.DataFrame({
            'type' : (self.n['atoms']) * [self.types['bond']],
            'atom_1' : list(range(1, self.n['atoms'] + 1)),
            'atom_2' : list(range(2, self.n['atoms'] + 1)) + [1],
        })

        return data

    def generate_angles(self):
        data = pd.DataFrame({
            'type' : (self.n['atoms']-2) * [self.types['angle']],
            'atom_1' : list(range(1, self.n['atoms']-1)),
            'atom_2' : list(range(2, self.n['atoms'])),
            'atom_3' : list(range(3, self.n['atoms']+1)),
        })

        return data

    def generate_charges(self, positions):

        # store length in a variable
        length = len(positions)

        # define the generation functions

        def _all():
            return [charge] * length

        def _random():
            _n = len(positions)
            _ratio = self._item['charge']['ratio']
            if _ratio > 1.0:
                raise ValueError('Charge ratio ({}) is greater than 1.'.format(_ratio))
            _n_charges = int((1-_ratio) * _n)
            indexes = sorted(random.sample(list(range(_n)), _n_charges))
            charge_list = [charge] * length
            for i in indexes:
                charge_list[i] = 0.0
            return charge_list

        def _diblock():
            return [charge] * length

        # store the functions in a hash table
        functions = {
                'all' : _all,
                'random' : _random,
                'diblock-regular' : _diblock
            }

        # manage charges if they exist
        charge = 0
        _c = self._item.get('charge', False)
        if _c:
            charge += _c['max']
            return functions[self._item['charge'].get('style', 'all')]()

        # default should be charge is zero of all the particle
        else:
            return functions['all']()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
