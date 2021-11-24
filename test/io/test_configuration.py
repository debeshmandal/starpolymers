from pathlib import Path

from starpolymers.io.configuration import *

from starpolymers import ConfigFile
from starpolymers import System

from starpolymers.molecules import LinearPolyelectrolyte

FOLDER = Path(__file__).parent

def test_configuration(delete: bool = True):
    item = {
        'molecule': 'dna',
        'lam': 10
    }
    molecule = LinearPolyelectrolyte(item)
    system = System(10., molecules=[molecule])

    # try all styles except dihedral
    for style in STYLES:
        if style == 'dihedral': continue
        config_file = ConfigFile(system, style=style)

        OUT = FOLDER / f'lammps.{style}.conf'
        OUT.unlink(missing_ok=True)
        config_file.write(OUT)

        if delete:
            OUT.unlink()

    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--delete', action='store_true')
    test_configuration(parser.parse_args().delete)