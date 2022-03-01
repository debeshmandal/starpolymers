#!/usr/bin/env python
"""Example script for generating a LAMMPS system with many of the same
topology with no overlaps.
"""
from starpolymers import generators
import subprocess

from starpolymers.io.configuration import ConfigFile
from starpolymers.generators.system import System
from starpolymers.molecules import MoleculeFactory
from softnanotools.logger import Logger
logger = Logger('MAIN')

import config

def generate(n: int, box: float, fname: str = 'lammps.test.conf'):
    molecule_list = [config.ConfigurationContainer.star] * n
    molecule_list.append(
        config.ConfigurationContainer.salt
    )
    molecules = MoleculeFactory(molecule_list, box=box).molecules
    system = System(box, molecules=molecules)
    config_file = ConfigFile(
        system,
        comment=f"Example system with {n} star polymers"
    )
    config_file.write(fname=fname, angles=False)
    return

def run(lammps: str, input_file: str = 'lammps.test.in'):
    cmd = lammps.split() + ['-i', input_file]
    subprocess.check_output(cmd)
    return

def main(**kwargs):
    generate(5, 50.0)
    run('mpirun -np 2 /home/debesh/programs/lammps/build/lmp -sf gpu -pk gpu 1')
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    main(**vars(parser.parse_args()))
