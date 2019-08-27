import numpy as np
import pandas as pd

markers = ['o','d','*','P','h','x','+','p','D','v','^','>','<','s','X']
base = ['blue', 'red', 'green', 'orange']
cdict = {'blue': ((0.024, 0.451, 0.573),
                  (0.043, 0.365, 0.455),
                  (0.031, 0.294, 0.369),
                  (0.008, 0.204, 0.259),
                  (0.004, 0.118, 0.153)),

         'red' : ((0.831, 0.004, 0.039),
                  (0.659, 0.039, 0.208),
                  (0.533, 0.027, 0.165),
                  (0.373, 0.000, 0.102),
                  (0.220, 0.000, 0.059)),

         'green': ((0.341, 0.812, 0.004),
                   (0.290, 0.643, 0.039),
                   (0.231, 0.522, 0.027),
                   (0.153, 0.365, 0.000),
                   (0.009, 0.216, 0.000)),

         'orange': ((0.922, 0.467, 0.008),
                    (0.729, 0.388, 0.043),
                    (0.592, 0.314, 0.031),
                    (0.412, 0.208, 0.000),
                    (0.243, 0.122, 0.000))}

def label_generator(variables, parameters, units=None):
    dictionary = dict()
    for i in range(len(parameters)):
        label = ''
        params = parameters[i]
        for j in range(len(variables)):
            label += variables[j]
            label += '='
            label += str(params[j])
            if units != None:
                    label += units[j]
            if j+1 != len(variables):
                label += ', '
        dictionary[i+1] = label
    
    return dictionary

def colour(dataframe, order_by=0):
    colours = []
    params = dataframe.values[:, order_by]
    group, counts = np.unique(params, return_counts=True)
    if len(group)/len(params) == 1:
        return False
    for i in range(len(group)):
        c = base[i]
        for j in range(counts[i]):
            colours.append(cdict[c][j])
    return colours
