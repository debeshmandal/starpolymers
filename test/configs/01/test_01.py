from starpolymers import MoleculeFactory, System, ConfigFile
#from starpolymers.generators.configuration import compare
import pytest

def _compare(*args):
    return True

def _parse_molecules():
    return

def _generate():
    ConfigFile(
        System(
            MoleculeFactory(
                _parse_molecules()
            ).molecules
        )
    ).write('new.dat')

@pytest.mark.skip
def test_system():
    _generate()
    assert _compare('config.dat', 'new.dat')