import starpolymers.molecules.base as base
import starpolymers.molecules.brush as brush
import starpolymers.molecules._common as _common
import starpolymers.molecules.polyelectrolyte as polyelectrolyte

from starpolymers.molecules import MoleculeFactory

def test_base():
    return

def test_brush():
    return

def test__common():
    _common.Molecule({'molecule' : 'DNA'})
    _common.ParticleArray({'molecule' : 'star'})
    _common.registry
    return

def test_polyelectrolyte():
    polyelectrolyte.LinearPolyelectrolyte({
        'molecule': 'DNA',
        'kap': 1,
        'lam': 21,
        'charge_style': 'all',
        'charge_max': -1,
        'counterions': False,
        'angle_type': 1
    })
    #polyelectrolye.StarPolyelectrolyte({
    #    'molecule': 'star',
    #    'kap': 10,
    #    'lam': 3,
    #    'charge_style': 'all',
    #    'charge_max': 1,
    #    'central' : 'all',
    #    'counterions':False,
    #    'angle_type': 1
    #})
    
def test_MoleculeFactory():
    return