# logger.py - for reading LAMMPS logfiles

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

def line_begins_with(words, fname):

    linelist = []
    with open(fname, 'r') as f:
        for line in f:
            linelist.append(line)#.split(None, 1)[0])

    numberlist = []
    for word in words:
        idx = [i for i, item in enumerate(linelist) if re.search(word, item)]
        numberlist.append(idx)
        
    return numberlist

class LogReader():
    
    def __init__(self, ID):
        self.ID = str(ID)
        self.fname = 'results/{0}/log.{0}.txt'.format(ID)

    def change_path(self, path):
        self.fname = path

    def read_thermo(self):

        idxlist = line_begins_with(['Step', 'Loop'], self.fname)
        start = max(idxlist[0])
        end = max(idxlist[1])
        length = (end-start)-1

        # read file from LBW 2-3 and write to temp file

        thermo = pd.read_csv(self.fname, delim_whitespace=True,
                             header=0, skiprows=start,  nrows=length) # read tempfile

        # delete temp file

        return thermo

    def read_energy(self):

        # read_thermo

        data = self.read_thermo()
        data = data[['Step', 'TotEng', 'PotEng']]

        return data

    def fast_plot(self, kind, savefig=False):

        if kind == 'energy':
            data = self.read_energy()
            plt.plot(data['Step'], data['PotEng'], label='Potential Energy')
            plt.plot(data['Step'],data['TotEng'], label='Total Energy')
            plt.xlabel(r'Timestep $[ \tau ]$')
            plt.ylabel('Energy $[ k_B T ]$')
            plt.title('{}'.format(self.ID))
            plt.legend()
            if savefig==True:
                plt.savefig('{}-energy'.format(self.ID))
            plt.show()
                
        else:
            ValueError('No such fast plot for kind: {}'.format(self.ID))
