"""
Objects used for scripting the generation of LAMMPS input files
"""
class InputFile():
    """
    Takes a template LAMMPS input file and parses the lines
    to return entries that can be edited
    """
    def __init__(template, comment='# LAMMPS input file'):
        """
        File handler for LAMMPS input files.
        """
        self._entries = []
        self.settings = self._parse_template(template)
        

    def _parse_template(fname):
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
    
    def write(fout):
        """
        Write entries to file
        """
        return

class Entry():
    """
    Stores a setting used in a LAMMPS input file
    """
    def __init__(kind, params, name=None):
        """
        """
        self.kind = kind
        self.params = params
        self.name = name

    def __repr__(self):
        return "Entry({}[{}]: {})".format(
            self.kind,
            self.name,
            ' '.join([i for i in self.params])
        )

    def is_similar(self, other):
        """
        Returns true if two entries are substitutable

        >>> Entry('run', 1000).is_similar(Entry('run', 2000))
        True
        """
        if type(other) != Entry:
            raise TypeError

        if self.kind == other.kind:
            if other.name != None:
                if self.name == other.name:
                    return True
                else:
                    return False
            else:
                return True
                
