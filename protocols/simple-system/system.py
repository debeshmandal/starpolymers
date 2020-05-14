#!/usr/bin/env python
"""
Script to generate a simple configuration
"""
from starpolymers.generators.configuration import FileGenerator

def main():
    # first initiate the FileGenerator
    fg = FileGenerator(30, atom_masses=[1.0])

    # store molecule as dictionary
    star = {
        'molecule': 'star',
        'kap': 5,
        'lam': 10,
        'charge_style': 'all',
        'charge_max': 0,
        'central' : 'all',
        'counterions':False,
        'angle_type': 1
    }
        
    # pass list of molecules (in this case one) to FileGenerator
    system = [star]
    fg.write_system_to_file(system)

if __name__=='__main__':
    main()
