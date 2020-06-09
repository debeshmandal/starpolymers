import os.path as path
import re

"""
Objects used for scripting the generation of LAMMPS input files
"""
ROOT = '/'.join(path.abspath(__file__).split('/')[:-1])

class Entry(object):
    """
    Stores a setting used in a LAMMPS input file
    """
    def __init__(self, name, params):
        """
        >>> Entry('variable', 100, name='kap')
        Entry(variable[kap]: 100)
        """
        self.name = name
        self.params = params
        
    def __repr__(self):
        if self.name == None:
            name = ''
        else:
            name = self.name
        return "{}([{}]: {})".format(
            type(self).__name__,
            name,
            self.params
        )

    def is_similar(self, other):
        """
        Returns true if two entries are substitutable

        >>> Entry('run', 1000).is_similar(Entry('run', 2000))
        True
        """
        if self.kind == other.kind:
            if other.name != None:
                if self.name == other.name:
                    return True
                else:
                    return False
            else:
                return True
                
class Variable(Entry):
    def __init__(self, name, value):
        Entry.__init__(self, name, value)
        try:
            self.params = float(self.params)
            self.equals = 'equal'
        except ValueError:
            self.equals = 'string'

    @property
    def string(self):
        return "{} {} {} {}".format(
            type(self).__name__.lower(),
            self.name,
            self.equals,
            self.params
        )

    @property
    def regex(self):
        pattern = r"(variable {} {}) [A-z|0-9|\.]+(\n)".format(
            self.name,
            self.equals
        )
        replacement = r"\1 {}\2".format(self.params)
        return pattern, replacement

class Fix(Entry):
    def __init__(self):
        raise NotImplementedError

class Compute(Entry):
    def __init__(self):
        raise NotImplementedError

class Dump(Entry):
    def __init__(self):
        raise NotImplementedError

class Template():
    """
    Read a template input file, return an object that has parsed all of the 
    variables which are the only changeable things in a new file that can
    be written.
    """
    def __init__(self, f_template):
        self._string = ''
        self._variables = []
        self._fixes = []
        self._dumps = []
        self._computes = []

        self.read(f_template)
        
    def read(self, fname):
        with open(fname, 'r') as f:
            self._string += f.read()
            for i, line in enumerate(self._string.split('\n')):
                line_list = line.split()
                if len(line_list) == 0: continue

                if line_list[0] == 'variable':
                    self._variables.append(
                        Variable(
                            line_list[1],
                            line_list[3]
                        ))
        return

    def write(self, fout):
        with open(fout, 'w') as f:
            f.write(self._string)

    def update(self, entry_obj):
        if not isinstance(entry_obj, Entry):
            raise TypeError('entry_obj should be an Entry object')
        if isinstance(entry_obj, Variable):
            self._set_variable(entry_obj.name, entry_obj.params)
            #print entry_obj.regex[0], entry_obj.regex[1]
            #print re.search(entry_obj.regex[0], self._string)
            self._string = re.sub(
                    entry_obj.regex[0],
                    entry_obj.regex[1],
                    self._string
                )

    @property
    def variables(self):
        return self._variables
        
    def _set_variable(self, name, value):
        for variable in self._variables:
            if variable.name == name:
                variable.params = value
                print(('{} has been updated to {}'.format(variable.name, variable.params)))


    @variables.getter
    def variables(self):
        return self._variables

    @property
    def fixes(self):
        raise NotImplementedError

    @fixes.setter
    def fixes(self):
        raise NotImplementedError

    @fixes.getter
    def fixes(self):
        raise NotImplementedError

    @property
    def computes(self):
        raise NotImplementedError

    @computes.setter
    def computes(self):
        raise NotImplementedError

    @computes.getter
    def computes(self):
        raise NotImplementedError

    @property
    def dumps(self):
        raise NotImplementedError

    @dumps.setter
    def dumps(self):
        raise NotImplementedError

    @dumps.getter
    def dumps(self):
        raise NotImplementedError


class InputFile():
    """
    Takes a template LAMMPS input file and parses the lines
    to return entries that can be edited
    """
    def __init__(self, template, comment='# LAMMPS input file'):
        """
        File handler for LAMMPS input files.
        """
        self._entries = []
        self.settings = self._parse_template(template)
        

    def _parse_template(self, fname):
        """
        Read a template file and create Entry objects
        """
        # open file
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                # skip comments
                if line[0] == '#': continue
                line = line.split()
                kind = line[0]
                if kind in ['variable', 'fix', 'compute' 'dump']:
                    name = line[1]
                    params = line[2:]
                else:
                    name = None
                    params = line[1:]
                entry = Entry(kind, params, name=name)
                self.entries.append(entry)

    @property
    def entries(self):
        return self._entries

    @entries.setter
    def entries(self, entry, index=-1, skip=-1):
        """
        Reads through all entries and finds a similar entry to
        the one given.
        """
        count = 0
        for i, old_entry in enumerate(self._entries):
            if entry.is_similar(old_entry):
                if count != skip:
                    self._entries[i] = entry
                    return
        
        self._entries.insert(index, entry)

    @entries.getter
    def entries(self, index=0):
        return
    
    def write(self, fout):
        """
        Write entries to file
        """
        return

templates = {
    '00_test' : Template('{}/input_files/_00_test.in'.format(ROOT)),
    '01_master' : Template('{}/input_files/_01_master.in'.format(ROOT))
}