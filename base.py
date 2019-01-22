import numpy as np
import math
import pandas as pd

def base_gen(dims, z=0, spac=0.5, mol=1, shift=0):
    lines = []

    x = dims[0]
    y = dims[1]

    xlo = -x/2.0
    xhi = x/2.0

    ylo = -y/2.0
    yhi = y/2.0

    c = spac * math.sin(math.pi/3.0) # height of triangles


    ax = np.arange(xlo,xhi,0.5)        # arg 1 is xlo, arg 2 is xhi, arg 3 is distance between atoms [width]
    bx = np.arange(xlo+spac/2.0,xhi,0.5)     # same but for second row
    ay = np.arange(ylo,yhi,2*c)        # arg 1 is ylo, arg 2 is yhi, arg 3 is distance between atoms [height]
    by = np.arange(c,yhi,2*c)        # same but for second row

    xa, ya = np.meshgrid(ax, ay)
    xb, yb = np.meshgrid(bx, by)
    xx = np.hstack(np.concatenate((ya, yb)))
    yy = np.hstack(np.concatenate((xa, xb)))
    zz = np.zeros((len(xx))) + z    # z is z-value of plane

    #randomfunc = np.random.choice(np.arange(0,len(xx),1), 50)

    df = pd.DataFrame()
    df['x']=xx
    df['y']=yy
    df['z']=zz

    for i in range(len(df)):
        atom = i+1
        atom_type = 1
        charge = 0
        x = df.loc[i]['x']
        y = df.loc[i]['y']
        z = df.loc[i]['z']
        line = '{} {} {} {} {} {} {}\n'.format(atom+shift, mol, 
                                            atom_type, charge, 
                                            x, y, z)
        lines.append(line)

    output = str()
    for line in lines:
        output += line
    return [output, len(lines)]
        

#print base_gen([10, 10], z=10, spac=2)