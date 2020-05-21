from system import System
import datetime

class ConfigFile():
    def __init__(self, system):
        self._system = system

    def write(self, fname='config.dat'):
        with open(fname, 'w') as f:
            f.write(self.comments)
            f.write(self.headers)
            f.write(self.atoms)
            f.write(self.bonds)
            f.write(self.angles)

    @property
    def comments(self):
        first_line = 'LAMMPS config file [{}]'.format(
            str(datetime.datetime.today()).split('.')[0]
        )
        second_line = ''
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
            masses += "{} {}\n".format(i, mass)
        return masses

    @property
    def atoms(self):
        return '' 

    @property
    def bonds(self):
        return '' 

    @property
    def angles(self):
        return '' 