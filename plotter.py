# plotter.py
import numpy as np
import matplotlib.pyplot as plt

class PLOT():
    def __init__(self, fname):
        self.fname = fname
    
    def write(self):
        self.fig.savefig(self.fname)


class TWO_AXES_SHAREX(PLOT):
    def __init__(self, fname):
        PLOT.__init__(self, fname)
        self.fig, self.ax = plt.subplots(2, sharex=True)
        

class ADVANCED_1(PLOT):
    # aim of the init is to create parallel 2x1, 1x3 plots
    def __init__(self):
        return