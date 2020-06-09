import random

import numpy as np
import pandas as pd

from ._common import Molecule
from starpolymers.tools.geometry import translation

class LinearPolyelectrolyte(Molecule):
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
        starting_point = translation[0]
        direction = self.settings['direction'][0]
        length = self._item['lam']
        positions = np.array([
            starting_point + direction * self.settings['spacing'] * i for \
                i in range(length)
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
            'type' : (self.n['atoms']-1) * [self.types['bond']],
            'atom_1' : list(range(1, self.n['atoms'])),
            'atom_2' : list(range(2, self.n['atoms']+1)),
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

        charge = float(self._item['charge']['max'])
        length = len(positions)

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

        functions = {
            'all' : _all,
            'random' : _random,
            'diblock-regular' : _diblock
        }

        return functions[self._item['charge'].get('style', 'all')]()

class StarPolyelectrolyte(Molecule):
    def __init__(self, item, **kwargs):
        Molecule.__init__(self, item, **kwargs)

        self.kap = self._item['kap']
        self.lam = self._item['lam']

        self._atoms = self.generate_atoms()
        self._bonds = self.generate_bonds()
        self._angles = self.generate_angles()

    def generate_atoms(self):
        spacing = self.settings['spacing']
        direction = self.settings['direction']
        mol_length = self.lam * spacing
        x = []
        y = []
        z = []
        for i in range(self.kap):
            for j in range(self.lam):
                x.append((mol_length-(j*spacing))*direction[i][0])
                y.append((mol_length-(j*spacing))*direction[i][1])
                z.append((mol_length-(j*spacing))*direction[i][2])

        x.append(0.0)
        y.append(0.0)
        z.append(0.0)

        data = pd.DataFrame({
                'mol' : len(x) * [self.mol],
                'type' : len(x) * [self.types['atom']],
                'x' : x,
                'y' : y,
                'z' : z,
                'q' : self.generate_charges(x)
            })

        return data

    def generate_bonds(self):

        m_bonds = self.kap * self.lam
        data = {
            'type' : [],
            'atom_1' : [],
            'atom_2' : []
        }
            
        for i in range(m_bonds):
            atom_1 = i+1
            
            if (i+1) % self.lam == 0:      
                atom_2 = self.n['atoms']                       
            else:
                atom_2 = i+2

            data['type'].append(self._item.get('bond_type', 1))
            data['atom_1'].append(atom_1)
            data['atom_2'].append(atom_2)

        return pd.DataFrame(data)

    def generate_angles(self):
        lam = self.lam
        kap = self.kap
        data = {
            'type' : [],
            'atom_1' : [],
            'atom_2' : [],
            'atom_3' : []
        }

        for i in range(kap):
            for j in range(lam-2):
                data['type'].append(self._item.get('angle_type', 1))
                data['atom_1'].append(lam*i+1+j)
                data['atom_2'].append(lam*i+2+j)
                data['atom_3'].append(lam*i+3+j)

        for i in range(kap):
            data['type'].append(self._item.get('angle_type', 1))
            data['atom_1'].append((i+1)*lam-1)
            data['atom_2'].append((i+1)*lam)
            data['atom_3'].append(self.n['atoms'])

        for j in reversed(list(range(kap+1))):
            if j > 2:
                for i in range((kap-j),kap-1):
                    data['type'].append(self._item.get('angle_type', 1))
                    data['atom_1'].append(lam*(kap+1-j))
                    data['atom_2'].append(self.n['atoms'])
                    data['atom_3'].append(lam*(i+2))
            elif j == 2:
                data['type'].append(self._item.get('angle_type', 1))
                data['atom_1'].append(lam*(kap-1))
                data['atom_2'].append(self.n['atoms'])
                data['atom_3'].append(lam*kap)

        return pd.DataFrame(data)

    def generate_charges(self, positions):
        
        charge = float(self._item['charge']['max'])
        length = len(positions)

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
            params = self._item['charge']
            q = params['max']
            ratio = params['ratio']
            position = params['position'] # should be core or tail
            if position not in ['core', 'tail']:
                raise ValueError(
                    "item['charge']['position']"
                    " should be 'core' or tail but is {}".format(position)
                )
            n_charges_per_arm = int(ratio * self.lam)
            arm_list = np.array([charge] * self.lam)

            # the PEs tail is where the 0th index is
            if position == 'core':
                arm_list[:n_charges_per_arm] = 0.0
            elif position == 'tail':
                arm_list[n_charges_per_arm:] = 0.0
        
            arm_list = np.concatenate([arm_list] * self.kap)
            arm_list = np.concatenate([arm_list, np.array([params.get('central', q)])])
            assert len(arm_list) == len(positions)
            return arm_list

        functions = {
            'all' : _all,
            'random' : _random,
            'diblock-regular' : _diblock
        }

        return functions[self._item['charge'].get('style', 'all')]()
