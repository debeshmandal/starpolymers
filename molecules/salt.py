from _common import ParticleArray

class Salt(ParticleArray):
    def __init__(self, item):
        super(Salt, self).__init__(item)
        self._atoms = self.generate_atoms()
    
    def generate_atoms(self):
        return