# Star Generator Module

import numpy as np
import pandas as pd
import math
import os
import random

direction = np.array([[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1],
                      [-1, 0, 0],
                      [0, -1, 0],
                      [0, 0, -1],
                      [1/math.sqrt(2), 1/math.sqrt(2), 0],
                      [0, 1/math.sqrt(2), 1/math.sqrt(2)],
                      [1/math.sqrt(2), 0, 1/math.sqrt(2)],
                      [-1/math.sqrt(2), 1/math.sqrt(2), 0],
                      [0, -1/math.sqrt(2), 1/math.sqrt(2)],
                      [-1/math.sqrt(2), 0, 1/math.sqrt(2)],
                      [1/math.sqrt(2), -1/math.sqrt(2), 0],
                      [0, 1/math.sqrt(2), -1/math.sqrt(2)],
                      [1/math.sqrt(2), 0, -1/math.sqrt(2)],
                      [-1/math.sqrt(2), -1/math.sqrt(2), 0],
                      [0, -1/math.sqrt(2), -1/math.sqrt(2)],
                      [-1/math.sqrt(2), 0, -1/math.sqrt(2)]])

spacing = 2.0

translation = np.array([[0, 0, 0],
                        [0.2, 0.2, 0],
                        [2, 1.5, -0.5],
                        [1.2, -0.5, 0.5],
                        [-0.5, -0.5, -0.5]])

# ---- #

# for star: use dictionary to read properties? #

test_star = {'molecule': 'star',
                 'kap': 3,
                 'lam': 10,
                 'charge_style': 'all',
                 'charge_max': 1}

test_star2 = {'molecule': 'star',
                 'kap': 5,
                 'lam': 8,
                 'charge_style': 'all',
                 'charge_max': 1}


# for DNA: use dictionary to read properties? #

test_DNA = {'molecule': 'DNA',
                'kap': 1,
                'lam': 21,
                'charge_style': 'all',
                'charge_max': -1,
            'counterions': True}

# overall system should be list of item dictionaries #

test_system = [test_star, test_DNA]

dummy_item = {'molecule': 'dummy',
              'kap': 0,
              'lam': 0,
              'charge_style': 'all',
              'charge_max': 0}

# ---- #

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

def charge_gen(item):

    """

    Returns an float that represents the charge on a single bead.
    
    """
    
    if item['charge_style'] == 'random':
        return
    if item['charge_style'] == 'alternating':
        return
    if item['charge_style'] == 'all':
        return float(item['charge_max'])

def neutraliser(system):
    """

    Returns the charge imbalance of a system by summing the charges of all
    none salt items
    
    """

    sys_charge = int()
    for item in system:
        if item['molecule'] != 'salt':
            q = item['charge_max']
            n_atoms = MaxCalculator(item).atoms(system)
            sys_charge += n_atoms * -q
    return sys_charge
        
    

class MaxCalculator():

    def __init__(self, item):
        self.item = item
    
    def atoms(self, system):
        kap = self.item['kap']
        lam = self.item['lam']
        if self.item['molecule'] == 'star':
            if self.item['counterions'] == True:
                max_atoms = 2*(kap*lam+1)
            elif self.item['counterions'] == False:
                max_atoms = kap*lam + 1
        elif self.item['molecule'] == 'dummy':
            max_atoms = 0
        elif self.item['molecule'] == 'DNA':
            if self.item['counterions'] == True:
                max_atoms = 2*lam
            elif self.item['counterions'] == False:
                max_atoms = lam
        if self.item['molecule'] == 'salt':
            max_atoms = 2*self.item['concentration']
            if self.item['neutralise'] == True:
                max_atoms += abs(neutraliser(system))
        return max_atoms

    def bonds(self):
        kap = self.item['kap']
        lam = self.item['lam']
        if self.item['molecule'] == 'star':
            max_bonds = kap*lam
        elif self.item['molecule'] == 'dummy':
            max_bonds = 0
        elif self.item['molecule'] == 'DNA':
            max_bonds = lam-1
        elif self.item['molecule'] == 'salt':
            max_bonds = 0
        return max_bonds

    def angles(self):
        kap = self.item['kap']
        lam = self.item['lam']
        if self.item['molecule'] == 'star':
            max_angles = kap*(kap-3+2*lam)/2
        elif self.item['molecule'] == 'dummy':
            max_angles = 0
        elif self.item['molecule'] == 'DNA':
            max_angles = lam-2
        elif self.item['molecule'] == 'salt':
            max_angles = 0
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

    def __init__(self, box, fstyle):
        self.box = box
        self.fstyle = fstyle

    def write_comments(self, system):

        """

        Returns string that is formatted as the first few lines of a LAMMPS config data file

        """
        comments = str()
        kap = system[0]['kap']
        lam = system[0]['lam']
        first_line = str('Star Polymer with {} arms which are {} beads in length'.format(kap, lam))
        second_line = str('secondline')
        comments += str('# {}\n'.format(first_line))
        comments += str('# {}\n'.format(second_line))
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

        spac = spacing

        for item in system:
            HeadGen = MaxCalculator(item)
            MAX_length += item['lam'] * spac
            MAX_atoms += HeadGen.atoms(system)
            MAX_bonds += HeadGen.bonds()
            MAX_angles += HeadGen.angles()        
        
        n_atoms = MAX_atoms
        m_bonds = MAX_bonds
        l_angles = MAX_angles

        a_atom_types = 1
        b_bond_types = 1
        c_angle_types = 1

        xlo = -self.box
        xhi = self.box
        ylo = -self.box
        yhi = self.box
        zlo = -self.box
        zhi = self.box
        
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

        atom_mass = 1.0

        next_line = str()
        next_line += str("1 {}\n".format(atom_mass))
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
        box = self.box

        CUMU_atoms = int()        
        for i in range(system_index):
            CUMU_atoms += MaxCalculator(system[i]).atoms(system)
        
        atom_ID_shift = CUMU_atoms
        spac = spacing
        shift_length = item['lam'] * spac
        atom_pos_shift = translation
        lam = item['lam']
        molecule_id = system_index + 1
        if item['molecule'] != 'salt':
            counterions = item['counterions']
        else:
            counterions = False

        if item['molecule'] == 'star':

            kap = item['kap']           
            atom_type = 1
            mol_length = lam*spac
            n_atoms = kap * lam + 1 + atom_ID_shift

            for i in range(kap):
                for j in range(lam):
                    atom_id = j+1 + (i*lam) + atom_ID_shift
                    charge = charge_gen(item)
                    x_pos = (mol_length-(j*spac))*direction[i][0] + atom_pos_shift[system_index][0]
                    y_pos = (mol_length-(j*spac))*direction[i][1] + atom_pos_shift[system_index][1]
                    z_pos = (mol_length-(j*spac))*direction[i][2] + atom_pos_shift[system_index][2]
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
                                                                 atom_type, charge_gen(item),
                                                              atom_pos_shift[system_index][0],
                                                              atom_pos_shift[system_index][1],
                                                              atom_pos_shift[system_index][2]))
        
        if item['molecule'] == 'DNA':
            
            # create DNA atoms
            atom_type = 1
            mol_length = lam*spac
            n_atoms = lam
            for i in range(lam):
                atom_id = i+1 + atom_ID_shift
                charge = charge_gen(item)
                x_pos = (mol_length-(i*spac))*direction[0][0] + atom_pos_shift[system_index][0]
                y_pos = (mol_length-(i*spac))*direction[0][1] + atom_pos_shift[system_index][1]
                z_pos = (mol_length-(i*spac))*direction[0][2] + atom_pos_shift[system_index][2]
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
                charge = -1 * charge_gen(item)
                x_pos = (mol_length-(i*spac))*direction[0][0] + atom_pos_shift[system_index][0] + spac
                y_pos = (mol_length-(i*spac))*direction[0][1] + atom_pos_shift[system_index][1] + spac
                z_pos = (mol_length-(i*spac))*direction[0][2] + atom_pos_shift[system_index][2] + spac
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
        
            conc = item['concentration']
            for i in range(conc):
                for j in range(2):
                    atom_id = 2*i+1 + j + atom_ID_shift
                    if j == 0:
                        charge = 1
                    elif j == 1:
                        charge = -1
                    x_pos = random.random()*box
                    y_pos = random.random()*box
                    z_pos = random.random()*box
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

            if item['neutralise'] == True:
            
                MAX_charge = neutraliser(system)
                    
                # set n_neut
                n_neut = int(abs(MAX_charge)/item['charge_max'])

                # set charge_sign

                if MAX_charge > 0:
                    charge_sign = 1.0
                else:
                    charge_sign = -1.0
                atom_id_start = 2*conc + atom_ID_shift
                atom_type = 1
                for i in range(n_neut):
                    atom_id = atom_id_start + i+1
                    charge = charge_sign * item['charge_max']
                    x_pos = random.random()*box
                    y_pos = random.random()*box
                    z_pos = random.random()*box
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
                
        return atom_list
        
    def write_bonds(self, system, system_index):

        """

        Returns string that is formatted as the bonds section of a LAMMPS input data file

        The format is:

        bond-ID bond-type atom-ID_1 atom-ID_2

        """

        bond_list = str()
        item = system[system_index]

        CUMU_atoms = int()
        CUMU_bonds = int()
        for i in range(system_index):
            CUMU_atoms += MaxCalculator(system[i]).atoms(system)
            CUMU_bonds += MaxCalculator(system[i]).bonds()
            
        atom_ID_shift = CUMU_atoms
        bond_ID_shift = CUMU_bonds
        lam = item['lam'] 


        if item['molecule'] == 'star':

            kap = item['kap']           
            n_atoms = kap * lam + 1 + atom_ID_shift
            m_bonds = kap * lam
            
            for i in range(m_bonds):
                
                if (i+1) % lam == 0:
                    bond_ID = i+1 + bond_ID_shift
                    bond_type = 1
                    atom1 = i+1 + atom_ID_shift
                    atom2 = n_atoms                   
                    
                else:
                    bond_ID = i+1 + bond_ID_shift
                    bond_type = 1
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
                bond_type = 1
                atom1 = i+1 + atom_ID_shift
                atom2 = i+2 + atom_ID_shift
                next_line = str()
                next_line += str("{} ".format(bond_ID))
                next_line += str("{} ".format(bond_type))
                next_line += str("{} ".format(atom1))
                next_line += str("{}".format(atom2))
                next_line += "\n"
                bond_list += next_line

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

            kap = item['kap']
            n_atoms = kap * lam + 1 + atom_ID_shift

            # Central Atom is k in ijk angular potential (x-x-o)
            
            if item['central'] != 'centre' or 'none':
                for i in range(kap):
                    angle_ID = i+1 + angle_ID_shift
                    angle_type = 1
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
                    angle_type = 1
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
                angle_type = 1
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

    def create_filename(self, system):

        """

        Returns string that is the filename for the system

        """
        if self.fstyle == 'al':
            star = system[0]
            salt = system[1]
            filename = str('al_'+star['kap']+'_'+star['lam']+'_'+salt['conc'])
        if self.fstyle == 'ssr':
            f_kap = str(system[0]['kap'])
            f_lam = str(system[0]['lam'])
            f_conc = str(system[2]['concentration'])
            filename = 'ssr_'+f_kap+'_'+f_lam+'_'+f_conc+'.dat'
        if self.fstyle == 'svl':
            f_kap = str(system[0]['kap'])
            f_lam = str(system[0]['lam'])
            filename = 'svl_{}_{}.dat'.format(f_kap, f_lam)
        if self.fstyle == 'ca':
            f_kap = str(system[0]['kap'])
            f_lam = str(system[0]['lam'])
            f_ang = str(system[0]['central'])
            filename = 'ca_{}_{}_{}.dat'.format(f_kap, f_lam, f_ang)
        elif self.fstyle == None:
            if len(system) == 1:
                item = system[0]
                filename = item['molecule']+str(item['kap'])+'_'+str(item['lam'])+'.dat'
            elif len(system) == 3:
                f_kap = str(system[2]['kap'])
                f_lam = str(system[2]['lam'])
                f_conc = str(system[0]['concentration'])
                filename = 'pd_'+f_kap+'_'+f_lam+'_'+f_conc+'.dat'
            elif len(system) == 2 and system[0]['molecule'] == 'star':
                f_kap = str(system[0]['kap'])
                f_lam = str(system[0]['lam'])
                f_conc = str(system[1]['concentration'])
                filename = 'sl_'+f_kap+'_'+f_lam+'_'+f_conc+'.dat'
            else:
                filename = str('exp.dat')

        return filename

    def write_system_to_file(self, system):

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
        
        with open(self.create_filename(system), 'w') as f:
            f.write(self.write_comments(system))
            f.write(self.write_header(system))
            f.write('\nMasses\n\n')
            f.write(self.write_masses())
            f.write('\nAtoms\n\n')
            for i in range(len(system)):
                if i == 0:
                    f.write(self.write_atoms(system, i))
                else:
                    f.write(self.write_atoms(system, i))   
            f.write('\nBonds\n\n')
            for i in range(len(system)):
                if i == 0:
                    f.write(self.write_bonds(system, i))
                else:
                    f.write(self.write_bonds(system, i))
            f.write('\nAngles\n\n')
            for i in range(len(system)):
                if i == 0:
                    f.write(self.write_angles(system, i))
                else:
                    f.write(self.write_angles(system, i))                
        return
