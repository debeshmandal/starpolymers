import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def last_branch(Brush, branch_index):
    n = Brush.trunk['lam']
    if branch_index == 0:
        return n
    branches = Brush.branches[0:branch_index]
    for branch in branches:
        n += branch['lam']
    return n

class Brush():
    
    def __init__(self, trunk_length, mol=1, starting_position=np.array([0,0,0]),
                 direction='up', base_id=None, graft_type=1):
        self.trunk = {'lam': trunk_length,
                      'atoms': list(range(1, trunk_length+1))}
        self.branches = []
        self.n_atoms = trunk_length
        self.generate_atom_list()
        self.mol = mol
        self.q_branches = len(self.branches)
        self.start = starting_position
        self.direction=direction
        self.base_id = base_id
        self.graft_type = graft_type
        
    def generate_atom_list(self):
        self.atom_list = list(range(1, self.n_atoms+1))
        
    def change_starting_position(self, starting_position=np.array([0,0,0])):
        self.start = starting_position
        
    def create_branch(self, branch_site, branch_length):
        self.branches.append({'lam': branch_length,
                              'atoms': list(range(self.n_atoms+1, self.n_atoms+branch_length+1)),
                              'site': branch_site})
        self.n_atoms += branch_length
        self.q_branches += 1
        self.generate_atom_list()
        
    def write_atoms(self, shift=0):
        
        lines = []
        graft = True
        
        charge=0
        existing_branches=0
        direction=self.direction
        for atom in self.atom_list:
            if graft==True:
                atom_type=self.graft_type
            else:
                atom_type=1
            x = self.start[0]
            y = self.start[1]
            z = self.start[2]
            if atom in self.trunk['atoms']:
                x += 0
                y += 0
                if direction == 'up':
                    z += atom-1
                elif direction == 'down':
                    z -= atom-1
            for i in range(self.q_branches):
                existing_branches = 0
                if i != 0:
                    for j in self.branches[:-1]:
                        existing_branches += j['lam']
                if atom in self.branches[i]['atoms']:
                    x += atom-self.trunk['lam']-existing_branches
                    y += 0
                    if direction =='up':
                        z += self.branches[i]['site']-1
                    elif direction =='down':        
                        z -= self.branches[i]['site']-1
            line = "{} {} {} {} {} {} {}\n".format(atom+shift, self.mol, atom_type, charge, x, y, z)
            lines.append(line)
            graft = False
        output = str()
        for line in lines:
            output += line
        return [output, len(lines)]
        
    def write_bonds(self, bond_shift=0, atom_shift=0):
        
        lines = []
        
        bonds=0
        bond_type=1
        branches_bonds=0
        for branch in self.branches:
            branches_bonds+=branch['lam']-1
            
        sites_bonds = len(self.branches)

        bonds += self.trunk['lam']-1 # bonds in the trunk
        bonds += branches_bonds # bonds within the branches
        bonds += sites_bonds # bonds between branches 
        
        branches_bonded = 0
        site_bonds_made = 0
        
        
        
        for bond in range(1, bonds+1): #i.e. for i in [1,2,3,4,5 ...]
            
            # build trunk bonds
            
            if bond < self.trunk['lam']:
                atom_1 = list(range(1,self.trunk['lam']+1))[bond-1]
                atom_2 = list(range(1,self.trunk['lam']+1))[bond]
                
            # build bonds in each branch
            
            elif bond <= bonds-sites_bonds:
                
                atom_1 = last_branch(self, branches_bonded) + bond - self.trunk['lam'] + 1
                atom_2 = last_branch(self, branches_bonded) + bond - self.trunk['lam'] + 2
        
            # build connections between branches and bonds
            
            else:
                
                if self.branches[site_bonds_made].get('site')!= None:
                    atom_1 = self.branches[site_bonds_made]['site']
                    atom_2 = last_branch(self, site_bonds_made) + 1
                    site_bonds_made+=1
        
            line = "{} {} {} {}\n".format(bond + bond_shift, bond_type,
                                          atom_1 + atom_shift, 
                                          atom_2 + atom_shift)
            lines.append(line)
        if self.base_id != None:
            base_bond = "{} {} {} {}\n".format(bonds+1 + bond_shift, 
                                               bond_type,
                                               1 + atom_shift, 
                                               self.base_id)
            lines.append(base_bond)
        output = str()
        for line in lines:
            output+=line
        return [output, len(lines)]