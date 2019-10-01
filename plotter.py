# plotter.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

_letters = ['a','b','c','d','e','g','h','i','j','k','l']

class PLOT():
    def __init__(self, fname):
        self.fname = fname
    
    def write(self):
        self.fig.savefig(self.fname)

    def add_image(self, fname, target):
        img = mpimg.imread(fname)
        ax = self.im[target]
        ax.imshow(img)


    def label_image(self, target_ax, xy_list, xytext_list):
        for i in range(3):
            ax = self.im[i]
            ax.text(20, 20, _letters[i],
                    horizontalalignment='left',
                    verticalalignment='top')
            target_ax.annotate(_letters[i], xy=xy_list[i],
                               xytext=xytext_list[i],
                               arrowprops=dict(arrowstyle="->",
                                               connectionstyle="arc3"))

class TWO_AXES_SHAREX(PLOT):
    def __init__(self, fname):
        PLOT.__init__(self, fname)
        self.fig, self.ax = plt.subplots(2, sharex=True)

class BACK_TO_BACK(PLOT):
    def __init__(self, fname):
        PLOT.__init__(self, fname)
        self.fig, ((self.ax_L), (self.ax_R)) = plt.subplot(1, 2, sharey=True, sharex=True)
        self.ax_R.invert_xaxis()
        

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
                                           sharex=self.bot_ax)
        #self.top_ax.set_xticklabels([])
        self.im = [self.fig.add_subplot(self.gs[0:2, -1:]),
                   self.fig.add_subplot(self.gs[2:4, -1:]),
                   self.fig.add_subplot(self.gs[4:, -1:])]
        
        for i in range(3):
            ax = self.im[i]
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_yticks([])
            ax.set_xticks([])

class ADVANCED_2(PLOT):
    """
    
    Uses GridSpec to create a 3xn array attached to a 2x1 array
    The 3xn array has three axes in a column and the 2x1 array 
    has two axes in a column.

    """
    def __init__(self, fname):
        PLOT.__init__(self, fname)
        self.fig = plt.figure(constrained_layout=True)
        self.gs = gridspec.GridSpec(ncols=3, nrows=6, figure=self.fig)
        #self.top_ax.set_xticklabels([])
        self.sub_ax = [self.fig.add_subplot(self.gs[0:2, -1:]),
                       self.fig.add_subplot(self.gs[2:4, -1:]),
                       self.fig.add_subplot(self.gs[4:, -1:])]
        #self.top_ax.set_xticklabels([])
        self.im = [self.fig.add_subplot(self.gs[0:2, -1:]),
                   self.fig.add_subplot(self.gs[2:4, -1:]),
                   self.fig.add_subplot(self.gs[4:, -1:])]
        
        #self.top_ax.set_xticklabels([])
        self.ax = [self.fig.add_subplot(self.gs[0:2, -1:]),
                       self.fig.add_subplot(self.gs[2:4, -1:]),
                       self.fig.add_subplot(self.gs[4:, -1:])]
        
 
