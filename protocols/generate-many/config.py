#!/usr/bin/env python
"""Storing template configurations
"""
import json
from typing import List
import random

from softnanotools.logger import Logger
from starpolymers import molecules
logger = Logger("CONFIG")

from starpolymers.generators import System
from starpolymers.io.configuration import ConfigFile
from starpolymers.molecules import MoleculeFactory

class ConfigurationContainer:
    star = {
        'molecule': 'star',
        'kap': 10,
        'lam': 3,
        'charge' : {
            'max': 0,
            'style': 'all'
        },
        'central' : 'all',
        'counterions': False,
        'atom_type': 1,
    }

    dna = {
        'molecule': 'dna',
        'lam': 2,
        'charge' : {
            'max': 0,
            'style': 'all'
        },
        'central' : 'all',
        'counterions': False,
        'atom_type': 2,
    }

    salt = {
        'molecule': 'salt',
        'anion': 1,
        'cation': 1,
        'concentration': 10,
        'box': 20.0
    }

    def __init__(self, system: System):
        self.system = system
        self.items = []

    def add_stars(
        self,
        N: int,
        lam: int,
        kap: int,
    ):

        for i in range(N):
            item = ConfigurationContainer.star
            item['kap'] = kap
            item['lam'] = lam
            self.items.append(item)
        return

    def add_dna(
        self,
        N: int
    ):
        for i in range(N):
            self.items.append(ConfigurationContainer.dna)

    def initialise(self):
        random.shuffle(self.items)
        factory = MoleculeFactory(item_list=self.items, box=self.system._box)
        self.system.add_molecules(factory.molecules, distribute=True)
        return


def main(
    stars: int = 5,
    dna: int = 5,
    box: float = 20.0,
    lam: int = 5,
    kap: int = 4,
):
    system = System(
        box,
        atom_masses=[1, 1, 1, 1, 1],
        bond_types=2,
        angle_types=0
    )
    container = ConfigurationContainer(system)
    container.add_stars(stars, lam, kap)
    container.add_dna(dna)
    container.initialise()
    config_file = ConfigFile(container.system, comment='Test', include_angles=False)
    config_file.write('lammps.test.conf')
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('stars', type=int)
    parser.add_argument('dna', type=int)
    parser.add_argument('-b', '--box', default=20.0, type=float)
    parser.add_argument('-l', '--lam', default=5, type=int)
    parser.add_argument('-k', '--kap', default=4, type=int)
    main(**vars(parser.parse_args()))