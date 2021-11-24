from ..generators.system import System
import datetime
import time

import pandas as pd

STYLES = [
    'atomic',
    'bond',
    'angle',
    'dihedral',
    'full',
    'molecular',
    'charge',
]

class ConfigFile():
    def __init__(
        self,
        system: System,
        comment: str = '',
        style: str = 'full',
        include_bonds: bool = True,
        include_angles: bool = True,
        include_dihedrals: bool = False,
        include_charges: bool = False,
    ):
        # important attributes
        self._system = system
        self.comment = comment
        self.style = style

        # list of things to include
        # note that style will dynamically EXCLUDE some of these
        self.include_bonds = include_bonds
        self.include_angles = include_angles
        self.include_dihedrals = include_dihedrals
        self.include_charges = include_charges

        self.implement_style()

    def implement_style(self, style: str = None):
        """Changes all obj.include_* boolean storage units to ensure
        only the proper parts of the system are included"""
        if style:
            self.style = style

        def atomic():
            self.include_bonds = False
            self.include_angles = False
            self.include_dihedrals = False
            self.include_charges = False
            return

        def bond():
            # skip bonds - defined in class
            self.include_angles = False
            self.include_dihedrals = False
            self.include_charges = False
            return

        def angle():
            # skip bonds - defined in class
            # skip angles - defined in class
            self.include_dihedrals = False
            self.include_charges = False
            return

        def molecular():
            # skip bonds - defined in class
            # skip angles - defined in class
            # skip dihedrals - defined in class
            self.include_charges = False
            return

        def dihedral():
            # skip bonds - defined in class
            # skip angles - defined in class
            # skip dihedrals - defined in class
            self.include_charges = False
            raise NotImplementedError
            return

        def full():
            # skip bonds - defined in class
            # skip angles - defined in class
            # skip dihedrals - defined in class
            self.include_charges = True
            return

        def charge():
            self.include_bonds = False
            self.include_angles = False
            self.include_dihedrals = False
            self.include_charges = True
            return

        # create hash-table of style functions
        style_functions = {
            'atomic': atomic,
            'bond': bond,
            'angle': angle,
            'dihedral': dihedral,
            'full': full,
            'molecular': molecular,
            'charge': charge
        }

        # choose function and call
        style_functions[self.style]()

        return

    def write(
        self,
        fname='lammps.main.conf',
        bonds: bool = None,
        angles: bool = None,
        dihedrals: bool = None,
        style: str = None,
    ):

        if style:
            self.implement_style(style=style)

        with open(fname, 'w') as f:
            f.write(self.comments)
            f.write(self.header)
            f.write(self.masses)
            f.write(self.atoms)
            if bonds or self.include_bonds:
                f.write(self.bonds)
            if angles or self.include_angles:
                f.write(self.angles)
            if dihedrals or self.include_dihedrals:
                raise NotImplementedError


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
        if self.include_angles:
            header += str("{} angles\n".format(self._system.n['angles']))
        header += "\n"
        header += str("{} atom types\n".format(self._system.types['atom']))
        header += str("{} bond types\n".format(self._system.types['bond']))
        if self.include_angles:
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
            if self.include_charges:
                atoms += "{} {} {} {} {} {} {}\n".format(
                    i + 1,
                    table['mol'][i],
                    table['type'][i],
                    table['q'][i],
                    table['x'][i],
                    table['y'][i],
                    table['z'][i]
                )
            else:
                atoms += "{} {} {} {} {} {}\n".format(
                    i + 1,
                    table['mol'][i],
                    table['type'][i],
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