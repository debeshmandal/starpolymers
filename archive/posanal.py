## position based analysis of LAMMPS output data

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import math

# Function to read a LAMMPS output file and allocate positions to DataFrame

def get_results(k, l, step='100'):
    
    """
    get_results reads LAMMPS output files and returns a dataframe of
    all of the atoms by ID and their positions in cartesian space

    Use this for LAMMPS.timestep --> DataFrame conversions
    """
    
    folder = 'results'
    results = 'star' + str(k) + '_' + str(l)
    destination = os.path.join(folder, results)
    fname_dump = 'dump.' + results + '.' + step
    floc = os.path.join(destination, fname_dump)
    with open(floc, 'r') as f:
        contents = f.read()
    fname_new = results + '-' + step + '.csv'
    fdest_new = os.path.join(destination, fname_new)
    with open(fdest_new, 'w+') as f:
        f.write(contents)
    with open(fdest_new, 'r') as fin:
        contents_new = fin.read().splitlines(True)
    with open(fdest_new, 'w') as fout:
        fout.writelines(contents_new[9:])
    data = pd.read_csv(fdest_new, delimiter = ' ', header=None,encoding="utf-8-sig")
    data = data.drop([5],1)
    data = data.sort_values([0])
    return data

# Function to calculate distance between to rows in the DataFrame

def dist_calc(Ri, Rj):
    
    """
    dist_calc returns the distance between two sets of cartesian co-ordinates

    Use this for calculating distances between vectors in array form
    """

    distance = math.sqrt(sum(np.subtract(Ri, Rj)**2))
    return distance

# Create DataFrame of all distances between atoms

def dist_tabgen(kap, lam, step='100'):
    
    """
    dist_tabgen returns a dataframe of the distance between atoms for a
    star of length lam and number of arms kap, at a given timestep (default str(100))

    Use this to convert raw data from get_results to useful pair position data
    """

    raw_data = get_results(kap, lam, step)
        
    DistanceBetweenAtoms = pd.DataFrame()

    count = 0
    
    # column 1 = Atom1; column 2 = Atom2; column 3 = r; column 4 = r2

    for i in range(len(raw_data)):
        x_i = raw_data.iloc[i][2]
        y_i = raw_data.iloc[i][3]
        z_i = raw_data.iloc[i][4]
        vector_i = np.array([x_i, y_i, z_i])
        for j in range(len(raw_data)):
            temp = np.zeros([1,4])
            temp[0][0] = int(raw_data.iloc[i][0])
            temp[0][1] = int(raw_data.iloc[j][0])
            x_j = raw_data.iloc[j][2]
            y_j = raw_data.iloc[j][3]
            z_j = raw_data.iloc[j][4]
            vector_j = np.array([x_j, y_j, z_j])            
            temp[0][2] = dist_calc(vector_i, vector_j)
            temp[0][3] = (temp[0][2]) ** 2
            tempdf = pd.DataFrame(temp)
            DistanceBetweenAtoms = DistanceBetweenAtoms.append(tempdf, ignore_index=True)
    
    return DistanceBetweenAtoms

# Calculate Rg2

def Rg2_calc(kap, lam, step='100'):
    
    """
    Rg2_calc returns the radius of gyration for a star of length lam
    and number of arms kap, at a given timestep (default str(100))

    Use this for calculating the radius of gyration for a datafile
    """
    distance_sum = 0
    count = 0
    d = dist_tabgen(kap, lam, step)
    N = d[0].max()
    for i in d:
        count = count + 1
        distance_sum = distance_sum + d[i][4]
    radius_of_gyration = (distance_sum / count)/(N**2)
    
    return radius_of_gyration  

# Calculate g(r)

def raddist_calc(kap, lam, step='100', resolution=10):
    
    """
    raddist_cal returns a dataframe of the radial distribution for a group
    of cartesian co-ordinates

    Use this to return a DataFrame that is equal to g(r) for a polymer
    """

    # Choose bin resolution (set to 100), initialise bin DataFrame using this resolution

    bin_res = resolution
    size = 1.0

    d = dist_tabgen(kap, lam, step)
    DistancesFromCentre = pd.DataFrame()
    N = d[0].max()
    DistancesFromCentre = d.loc[d[0] == N]
    DistancesFromCentre = DistancesFromCentre.reset_index(drop=True)
            
    # find largest distance between two atoms in d DataFrame

    MAX_length = DistancesFromCentre[2].max()

    print 'Maximum length from the central atom =', MAX_length

    # Use this distance as a normaliser... and scale to use integers corresponding to bin index

    scaled_distances = DistancesFromCentre[2] * (1/MAX_length) * bin_res
    
    # for all i in scaled_d, for all j in bins, if

    bins = np.ndarray([bin_res+1,1])
    bins.fill(0)

    for i in range(len(scaled_distances)):
        for j in range(bin_res+1):
            try:
                if int(scaled_distances[i]+size) == j:
                    bins[j] = bins[j]+1
                if int(scaled_distances[i]-size) == j:
                    bins[j] = bins[j]+1
            except:
                print 'error'

    g_r = pd.DataFrame(bins)
    g_r = g_r.reset_index(drop=False)
    g_r['index'] = g_r['index'] * (MAX_length/100)
    g_r = g_r.rename(columns={'index':'r',0:'g(r)'})

    return g_r





