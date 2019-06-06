def _create(star, dna, lower, upper, k, steps, start, stop):
    group1_end = star
    group2_start = star+1
    group2_end = star+dna
    colvar = "colvar {{\n  name dist\n  distance {{\n    group1 {{atomNumbersRange 1-{}}}\n    group2 {{atomNumbersRange {}-{}}}\n}}\n  lowerBoundary {}\n  upperBoundary {}\n}}".format(group1_end, group2_start, group2_end, lower, upper)
    harmonic = "harmonic {{\n  colvars dist\n  forceConstant {}\n  centers {}\n  targetCenters{}\n  targetNumSteps {}\n  writeTIPMF on\n}}".format(k, start, stop, steps)
    return "{}\n\n{}".format(colvar, harmonic)

def _write(s, fout):
    with open(fout, 'w') as f:
        f.write(s)
    return

class Colvars():
    def __init__(self, star, dna, lower=0.0, upper=70.0, k=5.0, steps=1000000, start=0, stop=70):
        self.file=_create(star, dna, lower, upper, k, steps, start, stop)
    
    def write(self, fname='col.vars'):
        _write(self.file, fname)