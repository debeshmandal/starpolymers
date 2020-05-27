#!/usr/bin/env python
"""
"""
import os
import pandas as pd
from argparse import ArgumentParser
from generate_doe_table import generate_table
from config import generate_config

def generate(index):
    generate_table()
    generate_config()
    return

def simulate(command):
    os.system(command)
    return

if __name__ == '__main__':
    parser = ArgumentParser(description='Design of Experiments example script')
    parser.add_argument(
            "--simulate", 
            default='',
            help='use this to run a command to run the simulation as well'
        )
    parser.add_argument(
            "--all", 
            action = 'store_true',
            help='use this to run a simulation as well'
        )

    args = parser.parse_args()
    generate()
    if args.simulate:
        simulate(args.simulate)
