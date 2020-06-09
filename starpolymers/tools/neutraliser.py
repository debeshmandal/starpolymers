#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys


for i in range(1, len(sys.argv)):
    fname = sys.argv[i]

    # read number of atoms
    with open(fname, 'r') as f:
        for i in range(4):
            line = f.readline()
            #print line
            if i == 3:
                n_atoms = [int(s) for s in line.split() if s.isdigit()]

    start = 0
    # read from 2 lines after atoms
    with open(fname, 'r') as f:
        for number, line in enumerate(f, 1):
            #print line
            if 'Atoms' in line:
                start += number
                break  

    # until this line + number of atoms          

    data = pd.read_csv(fname, delim_whitespace=True, 
                    skiprows=start, nrows=n_atoms[0],
                    header=None)
    

    # get charge column
    # sum up

    overall = data.values[:,3].sum()

    # print

    if overall == 0:
        print(("\nSUCCESS: The overall charge of {} is {}".format(fname, overall)))
    else:
        print(("\nERROR: The overall charge of {} is {}".format(fname, overall)))
        print((data.tail(1)))    