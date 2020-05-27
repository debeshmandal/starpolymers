#!/usr/bin/env python
"""
Example file generation for a linear polyelectrolyte in a neutral system
"""
from starpolymers import MoleculeFactory, System, ConfigFile
import pytest
import json
import os.path as path

ROOT = '/'.join(path.abspath(__file__).split('/')[:-1])

def _compare(*args):
    return True

def _parse_molecules(fname):
    with open(fname, 'r') as f:
        data = json.load(f, encoding='utf-8')['molecules']
    return data

def _generate(fout):
    f_mol = '{}/molecules.json'.format(ROOT)
    box = 50
    lpe, salt = MoleculeFactory(_parse_molecules(f_mol), box=box).molecules
    system = System(50, molecules=[lpe, salt])
    system.neutralise(salt)
    config = ConfigFile(system, comment='Linear Polyelectrolyte in Neutral system')
    config.write(fout)

def test_02():
    fout = '{}/test_config.dat'.format(ROOT)
    _generate(fout)
    assert ConfigFile.validate(fout)
    assert ConfigFile.compare('{}/config.dat'.format(ROOT), fout)

if __name__ == '__main__':
    test_02()