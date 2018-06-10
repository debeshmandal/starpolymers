# logger.py - for reading LAMMPS logfiles

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def line_begins_with(word):

    line_number = 0

    return line_number

class LogReader():
    
    def __init__(self, ID):
        self.ID = str(ID)
        self.fname = 'log.{}.txt'.format(ID)

    def read_thermo(self):

        line_begins_with('fix') # LBW 1

        line_begins_with('Step') # LBW 2
        
        line_begins_with('test') # LBW 3 - line doesn't begin with a number

        # read file from LBW 2-3 and write to temp file

        thermo = pd.read_csv(self.fname) # read tempfile

        # delete temp file

        return thermo

    def read_energy(self):

        # read_thermo

        data = self.read_thermo()
        data = pd.DataFrame([
            data['Step'].values,
            data['PotEng'].values,
            data['TotEng'].values
            ])

        return data

    def fast_plot(self, kind, savefig=False):

        if kind == 'energy':
            data = self.read_energy()
            fig, ax = plt.subplot(1)
            ax.plot(data['Step'],data['PotEng'], label='Potential Energy')
            ax.plot(data['Step'],data['TotEng'], label='Kinetic Energy')
            ax.legend()
            if savefig==True:
                plt.savefig('{}-energy'.format(self.ID))
            plt.show()
                
        else:
            ValueError('No such fast plot for kind: {}'.format(self.ID))
