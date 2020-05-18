import random

import numpy as np
import pandas as pd

from _common import Molecule
from _common import PE_SPACING as spacing

class LinearPolyelectrolyte(Molecule):
    """
    Linear Polyelectrolyte with given length
    """
    def __init__(self, item):
        super(Molecule, self).__init__(item)

        self._atoms = self.generate_atoms()
        self._bonds = self.generate_bonds()
        self._angles = self.generate_angles()

    def generate_atoms(self):
        """
        Returns dataframe with the following columns
        """
        starting_point = self.settings['spacing']
        direction = self.settings['direction'][0]
        length = self._item['lam']
        positions = np.array([
            starting_point + direction * spacing * i for \
                i in range(length)
        ])

        data = pd.DataFrame({
            'mol' : len(positions) * [self.mol],
            'type' : len(positions) * [self.atom_type],
            'x' : positions[:, 0],
            'y' : positions[:, 1],
            'z' : positions[:, 2],
            'q' : self.generate_charges()
        })

        return data

    def generate_bonds(self):
        return

    def generate_angles(self):
        return

class StarPolyelectrolyte(Molecule):
    def __init__(self, item):
        super(Molecule, self).__init__(item)
        self._atoms = self.generate_atoms()
        self._bonds = self.generate_bonds()
        self._angles = self.generate_angles()

    def generate_atoms(self):
        return

    def generate_bonds(self):
        return

    def generate_angles(self):
        return

def central_centre_gen(n_atoms, kap, lam, angle_shift, atom_shift):

    """

    Returns string of angle topology for a LAMMPS config data file.
    
    """

    angle_list_central = str()
    n_start = lam+1
    k=0
    for j in reversed(range(kap+1)):
        if j > 2:
            for i in range((kap-j),kap-1):
                angle_ID = kap+1 + k + angle_shift
                angle_type = 1
                atom1 = lam*(kap+1-j) + atom_shift
                atom2 = n_atoms
                atom3 = lam*(i+2) + atom_shift
                next_line = str()
                next_line += str("{} ".format(angle_ID))
                next_line += str("{} ".format(angle_type))
                next_line += str("{} ".format(atom1))
                next_line += str("{} ".format(atom2))
                next_line += str("{}".format(atom3))
                next_line += "\n"
                angle_list_central += str(next_line)
                k+=1


    for j in reversed(range(kap+1)):
        if j == 2:
            angle_ID = kap+1+k + angle_shift
            angle_type = 1
            atom1 = lam*(kap-1) + atom_shift
            atom2 = n_atoms
            atom3 = lam*kap + atom_shift
            next_line = str()
            next_line += str("{} ".format(angle_ID))
            next_line += str("{} ".format(angle_type))
            next_line += str("{} ".format(atom1))
            next_line += str("{} ".format(atom2))
            next_line += str("{}".format(atom3))
            next_line += "\n"
            angle_list_central += str(next_line)

    return angle_list_central

def item_charge(item, system):
    """

    Returns list of atom numbers that should have a charge.
    The length of the list can be fed into the max calculator.
    
    
    """

    # get number of atoms

    n_atoms = MaxCalculator(item).atoms(system)

    # create list from range

    atom_list = []

    if item['charge_style'] == 'all':
        atom_list = range(1, n_atoms+1)

    if item['charge_style'] == 'random':

    # calculate number of charges
    # and sample the list

        atom_list = range(1, n_atoms+1)
        n_charges = int(n_atoms*item['charge_params']['ratio'])
        atom_list = random.sample(atom_list, n_charges)

    if item['charge_style'] == 'diblock-regular':
        
        arm_charges = item['lam'] * item['charge_params']['ratio']
        arm_charges = int(arm_charges)
        for i in range(item['kap']):
            arm_list = np.array(range(1, item['lam']+1)) + (i*item['lam'])
            if item['charge_params']['block_position'] == 'centre':
                atom_list.extend(list(arm_list[:arm_charges]))
            elif item['charge_params']['block_position'] == 'end':
                arm_charges_neg = -1 * arm_charges
                atom_list.extend(list(arm_list[arm_charges_neg:]))
        if item['charge_params']['centre']:
            atom_list.append(n_atoms)
            print "atom_list: {}, {}".format(atom_list, len(atom_list))
    return atom_list

def charge_gen(item, atom_number, charge_list):

    """

    Returns a float that represents the charge on a single bead.
    
    """
    if atom_number in charge_list:
        return float(item['charge_max'])
    else:
        return 0.0

    # use dictionary for complicated things
    #
    #
    # : {'style' : 'random'|'block',
    #    'homo' : True | False,
    #    'ratio' : 0.0 < x < 1.0,
    #    'arms' : pattern list form that recurs e.g. [0, 1, 0, 1]
    #    'blocks' : pattern but dict w/ block size e.g. {3: 1, 2: 0, 1: -1}
    #    'het' : {arm_index <int> : group <a-z> # check if groups exist
    #    'a' : 'blocks'|'arms' : {}|[] ... }
    #
    #    return float(charge)

def neutraliser(system):
    """

    Returns the charge imbalance of a system by summing the charges of all
    none salt items
    
    """

    sys_charge = int()
    for item in system:
        if item['molecule'] != 'salt':
            # find out what the total charge is
            if item['charge_style'] == 'all':
                q = item['charge_max']
                n_atoms = MaxCalculator(item).atoms(system)
                sys_charge += n_atoms * -q
            elif item['charge_style'] == 'random':
                q = item['charge_max']
                n_atoms = MaxCalculator(item).atoms(system)
                sys_charge += n_atoms * -q * item['charge_params']['ratio']
            elif item['charge_style'] == 'diblock-regular':
                q = item['charge_max']
                n_atoms = MaxCalculator(item).atoms(system)
                sys_charge += item['kap'] * -q * int(item['charge_params']['ratio'] * item['lam']) -1
    print "sys_charge: {}".format(sys_charge)
    return int(sys_charge)

def salt(item, system, neutralise=True):
    "Returns number of anions and cations for any given combinations of salt valencies"

    conc = item['concentration']            # e.g. 100
    cation = item['cation']                 # e.g. +2 ... +2
    anion = item['anion']                   # e.g. -1 ... -3

    lcm = np.lcm(abs(cation), abs(anion))     # e.g.  2 ...  6

    # figure out n_anions and n_cations

    n_anions = conc * lcm/abs(anion)         # e.g. 200...300
    n_cations = conc * lcm/abs(cation)         # e.g. 100...200

    if neutralise:
        MAX_charge = neutraliser(system)

        if MAX_charge > 0:
            charge = cation
                
        elif MAX_charge < 0:
            charge = anion

        # set n_neut
        n_neut = abs(MAX_charge)/abs(charge)

        if MAX_charge > 0:
            n_cations += n_neut
            extra = n_neut % cation
                
        elif MAX_charge < 0:
            n_anions += n_neut
            extra = n_neut % anion              
    return [n_anions, n_cations, extra]

class MaxCalculator():

    def __init__(self, item):
        self.item = item
    
    def atoms(self, system):
        
        if self.item['molecule'] == 'star':
            kap = self.item['kap']
            lam = self.item['lam']
            if self.item['counterions'] == True:
                
                max_atoms = 2*(kap*lam+1)
            elif self.item['counterions'] == False:
                max_atoms = kap*lam + 1
        elif self.item['molecule'] == 'dummy':
            max_atoms = 0
        elif self.item['molecule'] == 'DNA':
            lam = self.item['lam']
            if self.item['counterions'] == True:
                max_atoms = 2*lam
            elif self.item['counterions'] == False:
                max_atoms = lam
        if self.item['molecule'] == 'salt':

            max_atoms = salt(self.item, system)[0] + salt(self.item, system)[1]

            if salt(self.item, system)[2] != 0:
                max_atoms+=1
        return max_atoms

    def bonds(self):
        
        if self.item['molecule'] == 'star':
            kap = self.item['kap']
            lam = self.item['lam']
            max_bonds = kap*lam
        elif self.item['molecule'] == 'dummy':
            kap = self.item['kap']
            lam = self.item['lam']
            max_bonds = 0
        elif self.item['molecule'] == 'DNA':
            kap = self.item['kap']
            lam = self.item['lam']
            max_bonds = lam-1
        elif self.item['molecule'] == 'salt':
            max_bonds = 0
        elif self.item['molecule'] == 'base':
            max_bonds=0
        return max_bonds

    def angles(self):
        
        if self.item['molecule'] == 'star':
            kap = self.item['kap']
            lam = self.item['lam']
            max_angles = kap*(kap-3+2*lam)/2
        elif self.item['molecule'] == 'dummy':
            max_angles = 0
        elif self.item['molecule'] == 'DNA':
            kap = self.item['kap']
            lam = self.item['lam']
            max_angles = lam-2
        elif self.item['molecule'] == 'salt':
            max_angles = 0
        elif self.item['molecule'] == 'brush':
            max_angles=0
        elif self.item['molecule'] == 'base':
            max_angles=0
        return max_angles