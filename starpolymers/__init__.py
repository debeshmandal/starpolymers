from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .molecules import MoleculeFactory
from .generators import System, ConfigFile
