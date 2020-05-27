#!/usr/bin/env python
"""
"""
import os
import pandas as pd
from argparse import ArgumentParser
from generate_doe_table import generate_table
from config import generate_config, generate_input

def generate(index, seed=1994):
    table = generate_table(10, seed)
    row = table.iloc[index, :]
    generate_config(row, fout='config.dat')
    generate_input(row, fout='simulation.in')
    return

if __name__ == '__main__':
    parser = ArgumentParser(description='Design of Experiments example script')
    parser.add_argument('--index', default=0)
    parser.add_argument('--seed', default=1994)
    args = parser.parse_args()
    generate(args.index, args.seed)
