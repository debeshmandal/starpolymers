"""
Objects used for scripting the generation of LAMMPS input files
"""
class Template():
    """
    Read a template input file, return an object that has parsed all of the 
    variables which are the only changeable things in a new file that can
    be written.
    """
    def __init__(self, f_template):
        pass

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

class Entry(object):
    """
    Stores a setting used in a LAMMPS input file
    """
    def __init__(self, kind, params, name=None):
        """
        >>> Entry('variable', 100, name='kap')
        Entry(variable[kap]: 100)
        """
        self.kind = kind
        self.params = params
        self.name = name

    def __repr__(self):
        if self.name == None:
            name = ''
        else:
            name = self.name
        return "Entry({}[{}]: {})".format(
            self.kind,
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
                
