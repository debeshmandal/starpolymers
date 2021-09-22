#!/usr/bin/env python
"""Storing template configurations
"""
import json

from softnanotools.logger import Logger
logger = Logger("CONFIG")

from starpolymers.generators import System, ConfigFile
from starpolymers.molecules import MoleculeFactory

class ConfigurationContainer:
    star = {
        'molecule': 'star',
        'kap': 10,
        'lam': 3,
        'charge' : {
            'max': 1,
            'style': 'all'
        },
        'central' : 'all',
        'counterions': False,
        'angle_type': 1
    }

    salt = {
        'molecule': 'salt',
        'anion': 1,
        'cation': 1,
        'concentration': 10,
        'box': 20.0
    }

    __molecules__ = [star, salt]

    @property
    def summary(self):
        string = ""
        for mol in self.__molecules__:
            string += json.dumps(mol, indent=2)
        return string

def main():
    configs = ConfigurationContainer()
    logger.info(configs.summary)
    box = 20
    molecules = MoleculeFactory(configs.__molecules__, box=box)
    system = System(
        box,
        molecules=molecules.molecules
    )
    config_file = ConfigFile(system, comment='Test')
    config_file.write('lammps.test.conf')
    return

if __name__ == '__main__':
    main()