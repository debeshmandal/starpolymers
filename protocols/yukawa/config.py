#!/usr/bin/env python
"""config.py - auto-generated by softnanotools"""
from softnanotools.logger import Logger
logger = Logger('CONFIG')

def main(**kwargs):
    logger.info('Running config...')
    # insert code here
    logger.info('Done!')
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='config.py - auto-generated by softnanotools')
    main(**vars(parser.parse_args()))
