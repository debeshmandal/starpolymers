from ._common import ParticleArray
import numpy as np
import pandas as pd

class Salt(ParticleArray):
    def __init__(self, item, **kwargs):
        ParticleArray.__init__(self, item, **kwargs)
         # charge on each species
        self.anion = float(self._item['anion'])
        self.cation = float(self._item['cation'])
        self._atoms = self.generate_atoms()
    
    def generate_atoms(self):
        """
        Returns dataframe with the following columns
        """

        q = self.generate_charges()

        positions \
            = np.random.uniform(size=(q.shape[0], 3)) \
            * self._kwargs['box']

        data = pd.DataFrame({
            'mol' : len(positions) * [self.mol],
            'type' : len(positions) * [self.types['atom']],
            'x' : positions[:, 0],
            'y' : positions[:, 1],
            'z' : positions[:, 2],
            'q' : q
        })

        return data
    
    def generate_charges(self):

        # lowest common multiple
        lcm = np.lcm(abs(int(self.cation)), abs(int(self.anion)))

        # number of each species
        n_anions = int(self._item['concentration'] * lcm/abs(self.anion))
        n_cations = int(self._item['concentration']  * lcm/abs(self.cation))

        charge_array = np.concatenate(
                [
                    np.ones((n_anions, 1)) * -self.anion,
                    np.ones((n_cations, 1)) * self.cation
                ]
            ).T[0]

        return charge_array