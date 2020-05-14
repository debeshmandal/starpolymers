from polyelectrolyte import StarPolyelectrolyte, LinearPolyelectrolyte
from salt import Salt

class MoleculeFactory():
    def __init__(self, item):
        options = {
            'star' : StarPolyelectrolyte,
            'dna' : LinearPolyelectrolyte,
            'salt' : Salt
        }
        self._function = options[item['molecule']]
    
    @property
    def molecule(self):
        return self._function(self.item)