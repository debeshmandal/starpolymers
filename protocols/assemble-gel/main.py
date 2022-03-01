#!/usr/bin/env python
"""Combines the generate-many and squeeze protocols to create a system
that contains a high density of two different and opposite-ly charged
molecule types, and then simulates them to equilibrium
"""
from softnanotools.logger import Logger
logger = Logger('MAIN')

def main(*args, **kwargs):
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    main(**vars(parser.parse_args()))
