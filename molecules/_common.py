import pandas as pd

PE_SPACING = 1.0

class _Registry():
    def __init__(self):
        self.lam_list = ['star', 'DNA']
        self.salt_list = ['salt']
        self.molecule_list = self.lam_list + self.salt_list

registry = _Registry()

class AbstractMolecule:
    def __init__(self, item, molecule=1, start=1):
        self.molecule = item['molecule']
        self._atoms = pd.DataFrame()
        self._molecule = molecule

    @property
    def atoms(self):
        # import self._atoms as data
        data = self._atoms
    
        # check format of dataframe
        assert data.columns == ['type', 'q', 'x', 'y', 'z']

        # iterate over all lines
        lines = []
        for i in range(len(data)):
            # get values
            atom_type = data.loc[i]['type']
            charge = data.loc[i]['q']
            x = data.loc[i]['x']
            y = data.loc[i]['y']
            z = data.loc[i]['z']
            
            # create line
            line = '{} {} {} {} {} {} {}\n'.format(
                molecule,
                atom_type, 
                charge,
                x, 
                y, 
                z
            )
            lines.append(line)
        
        # export string
        string = ''.join(lines)
        return string

class ParticleArray(AbstractMolecule):
    def __init__(self, item):
        super().__init__(item)

class Molecule(AbstractMolecule):
    def __init__(self, item):
        super().__init__(item)

    @property
    def bonds(self):
        # import self._atoms as data
        data = self._atoms
    
        # check format of dataframe
        assert data.columns == ['id', 'type', 'atom_1', 'atom_2']
        
        # iterate over all lines
        lines = []
        for i in range(len(data)):
            # get values
            bond = data.loc[i]['id']
            bond_type = data.loc[i]['type']
            atom_1 = data.loc[i]['atom_2']
            atom_2 = data.loc[i]['atom_1']
            
            # create line
            line = '{} {} {} {}\n'.format(
                bond, 
                bond_type,
                atom_1, 
                atom_2
            )
            lines.append(line)
        
        # export string
        string = ''.join(lines)
        return string

    @property
    def angles(self):
        # import self._atoms as data
        data = self._atoms
    
        # check format of dataframe
        assert data.columns == ['id', 'type', 'atom_1', 'atom_2', 'atom_3']
        
        # iterate over all lines
        lines = []
        for i in range(len(data)):
            # get values
            angle = data.loc[i]['id']
            angle_type = data.loc[i]['type']
            atom_1 = data.loc[i]['atom_1']
            atom_2 = data.loc[i]['atom_2']
            atom_3 = data.loc[i]['atom_3']
            
            # create line
            line = '{} {} {} {} {}\n'.format(
                angle, 
                angle_type,
                atom_1, 
                atom_2,
                atom_3
            )
            lines.append(line)
        
        # export string
        string = ''.join(lines)
        return string

    @property
    def charge(self):
        return