"""
Script for generating LAMMPS config files for a set of molecules.
"""
import math
import os
import random
import datetime

import numpy as np
import pandas as pd

from starpolymers.molecules import base
from starpolymers.molecules.brush import Brush

from starpolymers.molecules import registry, MoleculeFactory

from starpolymers.tools.geometry import direction, translation

spacing = 1.0

def central_centre_gen(n_atoms, kap, lam, angle_shift, atom_shift):

    """

    Returns string of angle topology for a LAMMPS config data file.
    
    """

    angle_list_central = str()
    n_start = lam+1
    k=0
    for j in reversed(list(range(kap+1))):
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


    for j in reversed(list(range(kap+1))):
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
        atom_list = list(range(1, n_atoms+1))

    if item['charge_style'] == 'random':

    # calculate number of charges
    # and sample the list

        atom_list = list(range(1, n_atoms+1))
        n_charges = int(n_atoms*item['charge_params']['ratio'])
        atom_list = random.sample(atom_list, n_charges)

    if item['charge_style'] == 'diblock-regular':
        
        arm_charges = item['lam'] * item['charge_params']['ratio']
        arm_charges = int(arm_charges)
        for i in range(item['kap']):
            arm_list = np.array(list(range(1, item['lam']+1))) + (i*item['lam'])
            if item['charge_params']['block_position'] == 'centre':
                atom_list.extend(list(arm_list[:arm_charges]))
            elif item['charge_params']['block_position'] == 'end':
                arm_charges_neg = -1 * arm_charges
                atom_list.extend(list(arm_list[arm_charges_neg:]))
        if item['charge_params']['centre']:
            atom_list.append(n_atoms)
            print(("atom_list: {}, {}".format(atom_list, len(atom_list))))
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
    print(("sys_charge: {}".format(sys_charge)))
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

def make_brush(item, mol=1):
    brush = Brush(item['trunk']['lam'], mol=mol,
                  starting_position=item['start'],
                  direction=item['direction'], base_id=item.get('base_id'),
                  graft_type=item.get('graft_type', 1))
    for branch in item['branches']:
        brush.create_branch(branch['site'],branch['lam'])
    return brush

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

        if self.item['molecule']=='brush':
            brush=make_brush(self.item)
            max_atoms = brush.write_atoms()[1]
        if self.item['molecule']=='base':
            max_atoms = base.base_gen(self.item['dims'],
                                      spac=self.item['spacing'])[1]
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
        elif self.item['molecule'] == 'brush':
            brush = make_brush(self.item)
            max_bonds = brush.write_bonds()[1]
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

    #def charges(self, system):
    #    if self.item['molecule'] != 'salt':
    #        charge_mag = self.item['charge_max']
    #        n_atoms = self.atoms(system)
    #        max_charge = n_atoms * charge_mag
    #    else:
    #        max_charge = 0
    #    return max_charge
                
class FileGenerator():
    """
    Input:

    system: list of dictionaries, with each item generally either a star polymer or DNA molecule:

    star = {'molecule': 'star',
            'kap': kap,
            'lam': lam,
            'charge_style': 'random' or 'alternating' or 'all'
            'charge_max': integer}

    DNA = {'molecule': 'DNA',
                'kap': kap
                'lam': int(number of base pairs),
                'charge_style': 'all'
                'charge': -1}
                
    """

    def __init__(
        self, 
        box, 
        atom_masses=[1.0], 
        bond_types=1, 
        angle_types=1,
        suppress_printing = False
        ):

        self._box = box
        self._atom_masses = atom_masses
        self._bond_types = bond_types
        self._angle_types = angle_types
        self._suppress_printing = suppress_printing

    def write_comments(self, system):

        """

        Returns string that is formatted as the first few lines of a LAMMPS config data file

        """
        comments = ''
        first_line = 'LAMMPS config file [{}]'.format(
            str(datetime.datetime.today()).split('.')[0]
        )
        second_line = 'Molecules: {}'.format(
            ', '.join([i['molecule'] for i in system[:3]])
        )
        comments += '# {}\n'.format(first_line)
        comments += '# {}\n'.format(second_line)
        comments += '\n'

        return comments
    
    def write_header(self, system):

        """

        Returns string that is formatted as the header of a LAMMPS config data file

        """
        
        MAX_length = float()
        MAX_atoms = int()
        MAX_bonds = int()
        MAX_angles = int()

        atom_type_list = list(range(1, len(self._atom_masses)+1))
        bond_type_list = list(np.array(list(range(self._bond_types)))+1)
        angle_type_list = list(np.array(list(range(self._angle_types)))+1)
            

        for item in system:
            HeadGen = MaxCalculator(item)
            #MAX_length += int(item['lam']) * spac
            MAX_atoms += HeadGen.atoms(system)
            MAX_bonds += HeadGen.bonds()
            MAX_angles += HeadGen.angles()

            try:
                atom_type_list.append(item['atom_type'])
            except:
                None
                
            try:
                bond_type_list.append(item['bond_type'])
            except:
                None
                
            try:
                angle_type_list.append(item['angle_type'])
            except:
                None
        
        n_atoms = MAX_atoms
        m_bonds = MAX_bonds
        l_angles = MAX_angles

        a_atom_types = max(atom_type_list)
        b_bond_types = max(bond_type_list)
        c_angle_types = max(angle_type_list)

        xlo = -self._box
        xhi = self._box
        ylo = -self._box
        yhi = self._box
        zlo = -self._box
        zhi = self._box
        
        header = str()

        next_line = str()
        next_line += str("{} atoms\n".format(n_atoms))
        next_line += str("{} bonds\n".format(m_bonds))
        next_line += str("{} angles\n".format(l_angles))
        next_line += "\n"
        next_line += str("{} atom types\n".format(a_atom_types))
        next_line += str("{} bond types\n".format(b_bond_types))
        next_line += str("{} angle types\n".format(c_angle_types))
        next_line += "\n"
        next_line += str("{} {} xlo xhi\n".format(xlo, xhi))
        next_line += str("{} {} ylo yhi\n".format(ylo, yhi))
        next_line += str("{} {} zlo zhi\n".format(zlo, zhi))
        header += next_line

        return header

    def write_masses(self):

        masses = str()

        for i in range(len(self._atom_masses)):
            next_line = str()
            next_line += str("{} {}\n".format(i+1, self._atom_masses[i]))
            masses += next_line

        return masses
        
    def write_atoms(self, system, system_index):

        """

        Returns string that is formatted as the atoms section of a LAMMPS input data file

        The format is:

        atom-ID mol-ID atom-type q x y z

        """

        atom_list = str()
        item = system[system_index]
        box = self._box
        charge_list = item_charge(item, system)

        # write function to check that item obeys rules

        CUMU_atoms = int()        
        for i in range(system_index):
            CUMU_atoms += MaxCalculator(system[i]).atoms(system)
        
        atom_ID_shift = CUMU_atoms
        atom_pos_shift = translation
        try:
            lam = item['lam']
        except:
            None
        molecule_id = system_index + 1
        if item['molecule'] in ['star', 'DNA']:
            counterions = item['counterions']
        else:
            counterions = False

        if item['molecule'] == 'star':
            # atom_list = Molecule(item).atoms

            kap = item['kap']           
            atom_type = 1
            mol_length = lam * spacing
            n_atoms = kap * lam + 1 + atom_ID_shift

            for i in range(kap):
                for j in range(lam):
                    atom_id = j+1 + (i*lam) + atom_ID_shift
                    
                    # check charge_ratio of item (if present), if not assume
                    # = 1
                    # if ratio < item[charge_ratio], do

                    charge = charge_gen(item, atom_id-atom_ID_shift, charge_list)
                    
                    # else continue
                    x_pos = (mol_length-(j*spacing))*direction[i][0] + atom_pos_shift[system_index][0]
                    y_pos = (mol_length-(j*spacing))*direction[i][1] + atom_pos_shift[system_index][1]
                    z_pos = (mol_length-(j*spacing))*direction[i][2] + atom_pos_shift[system_index][2]
                    next_line = str()
                    next_line += str("{} ".format(atom_id))
                    next_line += str("{} ".format(molecule_id))
                    next_line += str("{} ".format(atom_type))
                    next_line += str("{} ".format(charge))
                    next_line += str("{} ".format(x_pos))
                    next_line += str("{} ".format(y_pos))
                    next_line += str("{}".format(z_pos))
                    next_line += "\n"
                    atom_list += next_line
            atom_list += str("{} {} {} {} {} {} {} \n".format(n_atoms, molecule_id,
                                                              atom_type,
                                                              charge_gen(item, n_atoms, charge_list),
                                                              atom_pos_shift[system_index][0],
                                                              atom_pos_shift[system_index][1],
                                                              atom_pos_shift[system_index][2]))
        
        if item['molecule'] == 'DNA':
            
            # create DNA atoms
            atom_type = 1
            mol_length = lam*spacing
            n_atoms = lam
            for i in range(lam):
                atom_id = i+1 + atom_ID_shift
                charge = charge_gen(item, atom_id-atom_ID_shift, charge_list)
                x_pos = (mol_length-(i*spacing))*direction[0][0] + atom_pos_shift[system_index][0]
                y_pos = (mol_length-(i*spacing))*direction[0][1] + atom_pos_shift[system_index][1]
                z_pos = (mol_length-(i*spacing))*direction[0][2] + atom_pos_shift[system_index][2]
                next_line = str()
                next_line += str("{} ".format(atom_id))
                next_line += str("{} ".format(molecule_id))
                next_line += str("{} ".format(atom_type))
                next_line += str("{} ".format(charge))
                next_line += str("{} ".format(x_pos))
                next_line += str("{} ".format(y_pos))
                next_line += str("{}".format(z_pos))
                next_line += "\n"
                atom_list += next_line
            
            # create counterions
        if counterions == True:                
            for i in range(lam):
                atom_id = i+1 + atom_ID_shift + lam
                charge = -1 * charge_gen(item, atom_id-atom_ID_shift, charge_list)
                x_pos = (mol_length-(i*spacing))*direction[0][0] + atom_pos_shift[system_index][0] + spacing
                y_pos = (mol_length-(i*spacing))*direction[0][1] + atom_pos_shift[system_index][1] + spacing
                z_pos = (mol_length-(i*spacing))*direction[0][2] + atom_pos_shift[system_index][2] + spacing
                next_line = str()
                next_line += str("{} ".format(atom_id))
                next_line += str("{} ".format(molecule_id))
                next_line += str("{} ".format(atom_type))
                next_line += str("{} ".format(charge))
                next_line += str("{} ".format(x_pos))
                next_line += str("{} ".format(y_pos))
                next_line += str("{}".format(z_pos))
                next_line += "\n"
                atom_list += next_line

        if item['molecule'] == 'salt':
             
        # generates salt ions for a given concentration
            n_anions = salt(item, system)[0] 
            n_cations = salt(item, system)[1]

            # create cations

            for i in range(n_cations):
                atom_id = i+1 + atom_ID_shift
                charge = 1.0 * abs(item['cation'])
                x_pos = (random.random()-0.5)*box*2
                y_pos = (random.random()-0.5)*box*2
                z_pos = (random.random()-0.5)*box*2
                atom_type = 1
                next_line = str()
                next_line += str("{} ".format(atom_id))
                next_line += str("{} ".format(molecule_id))
                next_line += str("{} ".format(atom_type))
                next_line += str("{} ".format(charge))
                next_line += str("{} ".format(x_pos))
                next_line += str("{} ".format(y_pos))
                next_line += str("{}".format(z_pos))
                next_line += "\n"
                atom_list += next_line

            atom_ID_shift += n_cations

            # create anions

            for i in range(n_anions):
                atom_id = i+1 + atom_ID_shift
                charge = -1.0 * abs(item['anion'])
                x_pos = (random.random()-0.5)*box*2
                y_pos = (random.random()-0.5)*box*2
                z_pos = (random.random()-0.5)*box*2
                atom_type = 1
                next_line = str()
                next_line += str("{} ".format(atom_id))
                next_line += str("{} ".format(molecule_id))
                next_line += str("{} ".format(atom_type))
                next_line += str("{} ".format(charge))
                next_line += str("{} ".format(x_pos))
                next_line += str("{} ".format(y_pos))
                next_line += str("{}".format(z_pos))
                next_line += "\n"
                atom_list += next_line   

            atom_ID_shift += n_anions     

            # even out any charge discrepancies
            # 

            if salt(item, system)[2] != 0:
                atom_id = atom_ID_shift + 1
                x_pos = (random.random()-0.5)*box*2
                y_pos = (random.random()-0.5)*box*2
                z_pos = (random.random()-0.5)*box*2
                next_line = str()
                next_line += str("{} ".format(atom_id))
                next_line += str("{} ".format(molecule_id))
                next_line += str("{} ".format(atom_type))
                next_line += str("{} ".format(salt(item, system)[2]))
                next_line += str("{} ".format(x_pos))
                next_line += str("{} ".format(y_pos))
                next_line += str("{}".format(z_pos))
                next_line += "\n"
                atom_list += next_line

        if item['molecule'] == 'brush':
            brush = make_brush(item, mol=molecule_id)
            atom_list = brush.write_atoms(shift=atom_ID_shift)[0]

        if item['molecule'] == 'base':
            atom_list = base.base_gen(item['dims'], z=item['plane'],
                                      spac=item['spacing'], 
                                      mol=molecule_id, 
                                      shift=atom_ID_shift)[0]      
                
        return atom_list
        
    def write_bonds(self, system, system_index):

        """

        Returns string that is formatted as the bonds section of a LAMMPS input data file

        The format is:

        bond-ID bond-type atom-ID_1 atom-ID_2

        """

        bond_list = str()
        item = system[system_index]
        bond_type = item.get('bond_type', 1)

        CUMU_atoms = int()
        CUMU_bonds = int()
        for i in range(system_index):
            CUMU_atoms += MaxCalculator(system[i]).atoms(system)
            CUMU_bonds += MaxCalculator(system[i]).bonds()
            
        atom_ID_shift = CUMU_atoms
        bond_ID_shift = CUMU_bonds
        try:
            lam = item['lam'] 
        except:
            None

        if item['molecule'] == 'star':
            # bond_list = Molecule(item).bonds

            kap = item['kap']           
            n_atoms = kap * lam + 1 + atom_ID_shift
            m_bonds = kap * lam
            
            for i in range(m_bonds):
                
                if (i+1) % lam == 0:
                    bond_ID = i+1 + bond_ID_shift
                    
                    atom1 = i+1 + atom_ID_shift
                    atom2 = n_atoms                   
                    
                else:
                    bond_ID = i+1 + bond_ID_shift
                    
                    atom1 = i+1 + atom_ID_shift
                    atom2 = i+2 + atom_ID_shift
                    
                next_line = str()
                next_line += str("{} ".format(bond_ID))
                next_line += str("{} ".format(bond_type))
                next_line += str("{} ".format(atom1))
                next_line += str("{}".format(atom2))
                next_line += "\n"
                bond_list += next_line

        if item['molecule'] == 'DNA':
            
            for i in range(lam-1):
                bond_ID = i+1 + bond_ID_shift
                
                atom1 = i+1 + atom_ID_shift
                atom2 = i+2 + atom_ID_shift
                next_line = str()
                next_line += str("{} ".format(bond_ID))
                next_line += str("{} ".format(bond_type))
                next_line += str("{} ".format(atom1))
                next_line += str("{}".format(atom2))
                next_line += "\n"
                bond_list += next_line

        if item['molecule'] == 'brush':
            brush = make_brush(item)
            bond_list = brush.write_bonds(bond_shift=bond_ID_shift,
                                          atom_shift=atom_ID_shift)[0]

        return bond_list

    def write_angles(self, system, system_index):

        """

        Returns string that is formatted as the angles section of a LAMMPS input data file

        The format is:

        angle-ID angle-type atom1 atom2 atom3

        """
        angle_list = str()
        item = system[system_index]

        CUMU_atoms = int()
        CUMU_angles = int()
        for i in range(system_index):
            CUMU_atoms += MaxCalculator(system[i]).atoms(system)
            CUMU_angles += MaxCalculator(system[i]).angles()
            
        atom_ID_shift = CUMU_atoms
        angle_ID_shift = CUMU_angles
        lam = item['lam']
        
        
        if item['molecule'] == 'star':
            # angle_list = Molecule(item).angles

            kap = item['kap']
            n_atoms = kap * lam + 1 + atom_ID_shift

            # Central Atom is k in ijk angular potential (x-x-o)
            
            if item['central'] != 'centre' or 'none':
                for i in range(kap):
                    angle_ID = i+1 + angle_ID_shift
                    angle_type = item['angle_type']
                    atom1 = (i+1)*lam-1 + atom_ID_shift
                    atom2 = (i+1)*lam + atom_ID_shift
                    atom3 = n_atoms
                    next_line = str()
                    next_line += str("{} ".format(angle_ID))
                    next_line += str("{} ".format(angle_type))
                    next_line += str("{} ".format(atom1))
                    next_line += str("{} ".format(atom2))
                    next_line += str("{}".format(atom3))
                    next_line += "\n"
                    angle_list += next_line

            # Central Atom is j in ijk angular potential (x-o-x) 
            
            if item['central'] != 'end' or 'none':
                angle_list += central_centre_gen(n_atoms, kap, lam, angle_ID_shift, atom_ID_shift)

            # Angular topology for arms
            
            for i in range(kap):
                for j in range(lam-2):
                    angle_ID = kap*(kap+1)/2 + j+1 + i*(lam-2) + angle_ID_shift
                    angle_type = item['angle_type']
                    atom1 = lam*i+1+j + atom_ID_shift
                    atom2 = lam*i+2+j + atom_ID_shift
                    atom3 = lam*i+3+j + atom_ID_shift
                    next_line = str()
                    next_line += str("{} ".format(angle_ID))
                    next_line += str("{} ".format(angle_type))
                    next_line += str("{} ".format(atom1))
                    next_line += str("{} ".format(atom2))
                    next_line += str("{}".format(atom3))
                    next_line += "\n"
                    angle_list += next_line
                    
        if item['molecule'] == 'DNA':
            for i in range(lam-2):
                angle_ID = i+1 + angle_ID_shift
                angle_type = item['angle_type']
                atom1 = i+1 + atom_ID_shift
                atom2 = i+2 + atom_ID_shift
                atom3 = i+3 + atom_ID_shift
                next_line = str()
                next_line += str("{} ".format(angle_ID))
                next_line += str("{} ".format(angle_type))
                next_line += str("{} ".format(atom1))
                next_line += str("{} ".format(atom2))
                next_line += str("{}".format(atom3))
                next_line += "\n"
                angle_list += next_line

        return angle_list

    def write_system_to_file(self, system, fname='config.dat'):

        """

        High-level function:

        Writes configdatafile for a system. Takes list of dictionaries as the input.

        star = {'molecule': 'star',
                'kap': kap,
                'lam': lam,
                'charge_style': 'random' or 'alternating' or 'all'
                'charge_max': integer}

        DNA = {'molecule': 'DNA',
                    'lam': int(number of base pairs),
                    'charge_style': 'all'
                    'charge': -1}

        """
        #print "Running system checks ... \n"
        #for i in range(len(system)):
        #    print "Checking item index {} ...".format(i)
        #    if check_item(system[i]):
        #        print "Item is valid ... \n"
        #    else:
        #        print "Error in system, see above messages"
        #        return
        #print "System is valid, proceeding to write ..."
            
        with open(fname, 'w') as f:
            f.write(self.write_comments(system))
            f.write(self.write_header(system))

            # write masses
            f.write('\nMasses\n\n')
            f.write(self.write_masses())

            # write atoms
            f.write('\nAtoms\n\n')
            for i in range(len(system)):
                if i == 0:
                    f.write(self.write_atoms(system, i))
                else:
                    f.write(self.write_atoms(system, i))   
            
            # write bonds
            if self._bond_types > 0: 
                f.write('\nBonds\n\n')
                for i in range(len(system)):
                    if i == 0:
                        f.write(self.write_bonds(system, i))
                    else:
                        f.write(self.write_bonds(system, i))

            # write angles
            if self._angle_types > 0:
                f.write('\nAngles\n\n')
                for i in range(len(system)):
                    if i == 0:
                        f.write(self.write_angles(system, i))
                    else:
                        f.write(self.write_angles(system, i))

        if not self._suppress_printing:
            print(("Writing complete for {}".format(fname)))
