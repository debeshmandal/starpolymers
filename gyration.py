# module to read LAMMPS gyr.out data
# generated through the fix ave/time command

import pandas as pd
import numpy as np

class Gyration():
    def __init__(self, fname):
        self.fname = fname

    def set_root(self, root):
        self.root = root

    def set_fout(self, fout):
        self.fout = fout

    def get(self):
        data = pd.read_csv(self.fname, comment='#',
                           delim_whitespace=True, header=None)
        items = int(data.loc[0,1])
        steps = data[0].values[0::items+1]
        master = pd.DataFrame()
        master['ts'] = steps

        for i in range(items):
            master[str(i+1)] = data[data[0]==i+1][1].values

        self.master=master
        return self.master

    def write(self, fout='out.csv'):
        data = self.get()
        try:
            data.to_csv(self.fout, index_col=False)
        except:
            data.to_csv(fout, index_col=False)

    
        
        
