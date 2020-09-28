from starpolymers.generators.input_file import templates, Variable
from starpolymers.generators.configuration import ConfigFile, System
from starpolymers.generators.colvars import Colvars
from starpolymers.molecules import read_json, MoleculeFactory
from starpolymers.tools import Logger

logger = Logger(__name__)

def generate_input(row, fout='simulation.in'):
    input_file = templates['01_master']
    for item in ['k1', 'k2', 'K', 'theta0']:
        input_file.update(Variable(item, row[item]))
    input_file.write(fout)
    return input_file

def generate_config(row, fout='config.dat'):
    box = 100
    molecules = read_json('molecules.json', box=box)
    molecules.insert(0, MoleculeFactory([{
            'molecule' : 'star',
            'kap' : row['kap'],
            'lam' : row['lam'],
            'charge' : {
                'style' : 'diblock-regular',
                'max' : row['charge_max'],
                'ratio' : row['charge_ratio'],
                'position' : row['charge_position']
            }   
        }]).molecules[0]
    )
    system = System(box, molecules=molecules, atom_masses=[1.0, 1.0, 1.0], bond_types=2, angle_types=2, threshold=0.001)
    ConfigFile(system).write(fout)
    Colvars(row['kap']*row['lam']+1, 21, upper=70.0, k=2.0, stop=70, steps=1000000).write(fname='col.vars')

    return system
