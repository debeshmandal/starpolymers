#!/usr/bin/env python
"""Example script for generating a LAMMPS system with many of the same
topology with no overlaps.
"""
from softnanotools.logger import Logger
logger = Logger('MAIN')

def main(*args, **kwargs):
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    main(**vars(parser.parse_args()))
