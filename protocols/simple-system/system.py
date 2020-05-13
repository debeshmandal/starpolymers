#!/usr/bin/env python
"""
Script to generate a simple configuration
"""
from starpolymers.generators.configuration import FileGenerator

def main():
    fg = FileGenerator(30, atom_masses=[1.0])

    star = {'molecule': 'star',
            'kap': 5,
            'lam': 10,
            'charge_style': 'all',
            'charge_max': 0,
            'central' : 'all',
            'counterions':False,
            'angle_type': 1}

    system = [star]
    fg.write_system_to_file(system)

if __name__=='__main__':
    main()
