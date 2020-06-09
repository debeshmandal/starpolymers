from .system import System
import datetime

import pandas as pd

class ConfigFile():
    def __init__(self, system, comment=''):
        self._system = system
        self.comment = comment

    def write(self, fname='config.dat'):
        with open(fname, 'w') as f:
            f.write(self.comments)
            f.write(self.header)
            f.write(self.masses)
            f.write(self.atoms)
            f.write(self.bonds)
            f.write(self.angles)

    @property
    def comments(self):
        first_line = 'LAMMPS config file [{}]'.format(
            str(datetime.datetime.today()).split('.')[0]
        )
        second_line = self.comment
        comments = '# {}\n# {}\n\n'.format(first_line, second_line)
        return comments

    @property
    def header(self):
        box = self._system.box
        header = ''
        header += str("{} atoms\n".format(self._system.n['atoms']))
        header += str("{} bonds\n".format(self._system.n['bonds']))
        header += str("{} angles\n".format(self._system.n['angles']))
        header += "\n"
        header += str("{} atom types\n".format(self._system.types['atom']))
        header += str("{} bond types\n".format(self._system.types['bond']))
        header += str("{} angle types\n".format(self._system.types['angle']))
        header += "\n"
        header += str("{} {} xlo xhi\n".format(box['xlo'], box['xhi']))
        header += str("{} {} ylo yhi\n".format(box['ylo'], box['yhi']))
        header += str("{} {} zlo zhi\n".format(box['zlo'], box['zhi']))
        return header
    
    @property
    def masses(self):
        masses = '\nMasses\n\n'
        for i, mass in enumerate(self._system._masses):
            masses += "{} {}\n".format(i + 1, mass)
        return masses

    @property
    def atoms(self):
        atoms = '\nAtoms\n\n' 
        table = self._system.atoms
        for i in range(len(table)):
            atoms += "{} {} {} {} {} {} {}\n".format(
                i + 1,
                table['mol'][i],
                table['type'][i],
                table['q'][i],
                table['x'][i],
                table['y'][i],
                table['z'][i]
            )
        return atoms

    @property
    def bonds(self):
        bonds = '\nBonds\n\n' 
        table = self._system.bonds
        for i in range(len(table)):
            bonds += "{} {} {} {}\n".format(
                i + 1,
                table['type'][i],
                int(table['atom_1'][i]),
                int(table['atom_2'][i])
            )
        return bonds

    @property
    def angles(self):
        angles = '\nAngles\n\n' 
        table = self._system.angles
        for i in range(len(table)):
            angles += "{} {} {} {} {}\n".format(
                i + 1,
                table['type'][i],
                table['atom_1'][i],
                table['atom_2'][i],
                table['atom_3'][i]
            )
        return angles

    @staticmethod
    def stats(fname):
        stats = {}

        with open(fname, 'r') as f:
            [f.readline() for i in range(3)]
            stats['n_atoms'] = int(f.readline().split()[0])
            stats['n_bonds'] = int(f.readline().split()[0])
            stats['n_angles'] = int(f.readline().split()[0])
            f.readline()
            stats['atom_types'] = int(f.readline().split()[0])
            stats['bond_types'] = int(f.readline().split()[0])
            stats['angle_types'] = int(f.readline().split()[0])

        skip_atoms = 21+stats['atom_types']
        skip_bonds = 24+stats['atom_types']+stats['n_atoms']
        skip_angles = 26+stats['atom_types']+stats['n_atoms']+stats['n_bonds']

        stats['atoms'] = pd.read_csv(
                fname, 
                delim_whitespace=True, 
                header=None, 
                skiprows=skip_atoms,
                nrows=stats['n_atoms']
            )

        stats['bonds'] = pd.read_csv(
                fname, 
                delim_whitespace=True, 
                header=None, 
                skiprows=skip_bonds,
                nrows=stats['n_bonds']
            )

        stats['angles'] = pd.read_csv(
                fname, 
                delim_whitespace=True, 
                header=None, 
                skiprows=skip_angles,
                nrows=stats['n_angles']
            )

        return stats

    @staticmethod
    def validate(fname):
        stats = ConfigFile.stats(fname)
        for name in ['atom', 'bond', 'angle']:
            header_n = stats['n_{}s'.format(name)]
            length_n = len(stats['{}s'.format(name)])
            if header_n != length_n:
                return False
        return True

    @staticmethod
    def compare(fname_1, fname_2):
        stats_1 = ConfigFile.stats(fname_1)
        stats_2 = ConfigFile.stats(fname_2)
        for key, value in list(stats_1.items()):
            try:
                if stats_2[key] != value:
                    return False
            except ValueError:
                if isinstance(value, pd.DataFrame):
                    continue
                else:
                    raise ValueError
        return True