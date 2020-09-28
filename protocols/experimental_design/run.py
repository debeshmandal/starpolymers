#!/usr/bin/env python
"""
Script to reproducibly generate a design of experiments table
and select a config/input file from the table and generate them
"""
import os
import pandas as pd
from argparse import ArgumentParser
from generate_doe_table import generate_table
from config import generate_config, generate_input
from starpolymers.tools import Logger

logger = Logger(__name__)

def generate(n, index, seed=1994):
    table = generate_table(n, seed)
    row = table.iloc[index, :]
    generate_config(row, fout='config.dat')
    generate_input(row, fout='simulation.in')
    return

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Design of Experiments example script'
    )
    parser.add_argument('--index', default=1, type=int)
    parser.add_argument('--seed', default=1994)
    parser.add_argument('-n', '--number', default=10, type=int)
    args = parser.parse_args()
    logger.info(
        f'Running experimental design protocol with:\n\tindex\t{args.index}'
        f'\n\tseed\t{args.seed}\n\tnumber\t{args.number}'
    )
    generate(args.number, args.index-1, args.seed)
