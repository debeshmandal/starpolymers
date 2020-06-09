from .polyelectrolyte import StarPolyelectrolyte, LinearPolyelectrolyte
from .salt import Salt
import json
from ._common import registry

class MoleculeFactory():
    def __init__(self, item_list, box=None):
        self.registry = {
            'star' : StarPolyelectrolyte,
            'dna' : LinearPolyelectrolyte,
            'salt' : Salt
        }
        self.molecules = []
        for item in item_list:
            mol = self.registry[item['molecule']](item, box=box)
            self.molecules.append(mol)

def read_json(fname, box=None):
    with open(fname, 'r') as f:
        data = json.load(f, encoding='utf-8')['molecules']
    molecules = MoleculeFactory(data, box).molecules
    return molecules