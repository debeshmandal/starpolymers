import numpy as np
import math
import pandas as pd

def base_gen(dims, z=0, spac=2, mol=1, shift=0, out='output'):
    lines = []

    x = dims[0]
    y = dims[1]

    xx = np.arange(-x, x, spac)
    yy = np.arange(-y, y, spac*math.sqrt(3))

    mesh = np.meshgrid(xx,yy)

    data=pd.DataFrame()
    for i in range(len(mesh[0])):

        # create first mesh in column form
        temp = pd.DataFrame()
        temp['x'] = mesh[0][i]
        temp['y'] = mesh[1][i]
        data = data.append(temp, ignore_index=True)

        # Do the same but accounting for a shift in the {1000} direction
        temp = pd.DataFrame()
        temp['x'] = mesh[0][i] + 0.5*spac
        temp['y'] = mesh[1][i] + 0.5*spac*math.sqrt(3)
        data = data.append(temp, ignore_index=True)


    data = data[abs(data['x'])<dims[0]]
    data = data[abs(data['y'])<dims[1]]
    data['z'] = np.zeros(len(data))+z

    if out=='df':
        return data

    for i in range(len(data)):
        atom = i+1
        atom_type = 1
        charge = 0
        x = data.loc[i]['x']
        y = data.loc[i]['y']
        z = data.loc[i]['z']
        line = '{} {} {} {} {} {} {}\n'.format(atom+shift, mol, 
                                            atom_type, charge, 
                                            x, y, z)
        lines.append(line)

    output = str()
    for line in lines:
        output += line
    return [output, len(lines)]
        

#print base_gen([10, 10], z=10, spac=2)