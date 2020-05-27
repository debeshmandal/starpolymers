import sys

# Check correct Python version is being used
version = int(sys.version.split('|')[0].strip().split('.')[0])
if version != 2:
    raise SystemError('Using Python {} when you must use Python 2'.format(version))

from molecules import MoleculeFactory
from generators import System, ConfigFile