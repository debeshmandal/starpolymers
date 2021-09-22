#!/usr/bin/env python
"""Example script for generating a single topology within a geometry
with a given shape and size e.g. a sphere with radius=5, a cube with
side=5 etc.
"""
from softnanotools.logger import Logger
logger = Logger('MAIN')

def main(*args, **kwargs):
    return

if __name__ == '__main__':
    import argparse
    parser = ArgumentParser()
    main(**vars(parser.parse_args()))
