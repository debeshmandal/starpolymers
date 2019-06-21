# plotter.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

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
    """

    Uses GridSpec to create 5 sets of axes, one wide column with two
    rows for plotting data and next to it on the right, a column with
    half of the width and 3 rows for images.

    """
    def __init__(self, fname):
        PLOT.__init__(self, fname)
        self.fig = plt.figure(constrained_layout=True)
        self.gs = gridspec.GridSpec(ncols=3, nrows=6, figure=self.fig)
        self.bot_ax = self.fig.add_subplot(self.gs[3:, :-1])
        self.top_ax = self.fig.add_subplot(self.gs[0:3, :-1], 
                                           sharex=self.top_ax)
        self.top_ax.set_xticklabels([])
        self.bot_ax = self.fig.add_subplot(self.gs[3:, :-1], 
                                           sharex=self.top_ax)
        self.im = {'top':self.fig.add_subplot(self.gs[0:2, -1:]),
                   'mid':self.fig.add_subplot(self.gs[2:4, -1:]),
                   'bot':self.fig.add_subplot(self.gs[4:, -1:])}

        for pos in self.im:
            ax = self.im[pos]
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_yticks([])
            ax.set_xticks([])

    
    def add_image(self, fname, target):
        img = mpimg.imread(fname)
        ax = self.im[target]
        ax.imshow(img)


        