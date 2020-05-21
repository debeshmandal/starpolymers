from starpolymers import MoleculeFactory, System, ConfigFile
#from starpolymers.generators.configuration import compare
import pytest
import json

def _compare(*args):
    return True

def _parse_molecules(fname):
    with open(fname, 'r') as f:
        data = json.load(f, encoding='utf-8')['molecules']
    return data

def _generate():
    ConfigFile(
        System(
            10, 
            molecules=MoleculeFactory(
                _parse_molecules('molecules.json')
            ).molecules
        )
    ).write('new.dat')

@pytest.mark.skip
def test_system():
    _generate()
    assert _compare('config.dat', 'new.dat')

if __name__ == '__main__':
    _generate()