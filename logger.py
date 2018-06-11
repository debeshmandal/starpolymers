# logger.py - for reading LAMMPS logfiles

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

def line_begins_with(words, fname):

    linelist = []
    with open(fname, 'r') as f:
        for line in f:
            linelist.append(line)

    numberlist = []
    for word in words:
        idx = [i for i, item in enumerate(linelist) if re.search(word, item)]
        numberlist.append(idx)
        
    return numberlist

class LogReader():
    
    def __init__(self, ID):
        self.ID = str(ID)
        self.fname = 'results/{0}/log.{0}.txt'.format(ID)

    def read(self, kind):

        idxlist = line_begins_with(['Step', 'Loop'], self.fname)
        start = max(idxlist[0])
        end = max(idxlist[1])
        length = (end-start)-1

        # read file from LBW 2-3 and write to temp file

        thermo = pd.read_csv(self.fname, delim_whitespace=True,
                             header=0, skiprows=start,  nrows=length) # read tempfile

        # delete temp file

        if kind == 'thermo':
            data = thermo

        elif kind == 'energy':
            energy = thermo[['Step', 'TotEng', 'PotEng']]
            data = energy

        elif kind == 'temp':
            temp = thermo[['Step', 'Temp']]
            data = temp

        return data


    def fast_plot(self, kind, savefig=False):

        data = self.read(kind)

        if kind == 'energy':
            plt.plot(data['Step'], data['PotEng'], label='Potential Energy')
            plt.plot(data['Step'],data['TotEng'], label='Total Energy')
            plt.ylabel('Energy $[ k_B T ]$')

        if kind == 'temp':
            plt.plot(data['Step'], data['Temp'], label='Temperature')
            plt.ylabel('Temperature $[T]$')
                
        else:
            raise ValueError('No such fast plot for kind: {}'.format(kind))
            return

        plt.xlabel(r'Timestep $[ \tau ]$')
        plt.title('{}'.format(self.ID))
        plt.legend()
        if savefig==True:
            plt.savefig('{}-{}'.format(self.ID, kind))
        plt.show()
