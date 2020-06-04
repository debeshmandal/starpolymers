"""
Script to generate Colvars files
"""
def _create(star, dna, lower, upper, k, steps, start, stop):
    """
    Create a string to be written to a Colvars file.

    Parameters
    ----------
    star : int 
        number of atoms in star polymer
    dna : int
        number of atoms in dna molecule
    lower : float
        starting point for the collective variables
    upper : float
        endpoint for collective variable
    k : float
        spring constant
    steps : int
        number of simulation steps
    start : float
        starting point for spring position
    stop : float
        end point for spring position

    Returns
    -------
    output_string : str
        string that can be written to file
    """
    group1_end = star
    group2_start = star+1
    group2_end = star+dna
    colvar = "colvar {{\n  name dist\n  distance {{\n    group1 {{atomNumbersRange 1-{}}}\n    group2 {{atomNumbersRange {}-{}}}\n}}\n  lowerBoundary {}\n  upperBoundary {}\n}}".format(group1_end, group2_start, group2_end, lower, upper)
    harmonic = "harmonic {{\n  colvars dist\n  forceConstant {}\n  centers {}\n  targetCenters {}\n  targetNumSteps {}\n  writeTIPMF on\n}}".format(k, start, stop, steps)
    return "{}\n\n{}".format(colvar, harmonic)

class Colvars():
    """
    Wrapper that creates a Colvars file string and can write it to file.
    
    Note that this assumes that the star argument is the first molecule and
    the DNA molecule is the second.

    Parameters
    ----------
    star : int 
        number of atoms in star polymer
    dna : int
        number of atoms in dna molecule
    lower : float
        starting point for the collective variables
    upper : float
        endpoint for collective variable
    k : float
        spring constant
    steps : int
        number of simulation steps
    start : float
        starting point for spring position
    stop : float
        end point for spring position
    """
    def __init__(
        self, 
        star, 
        dna, 
        lower=0.0, 
        upper=70.0, 
        k=5.0, 
        steps=1000000, 
        start=0, 
        stop=70):
        self.string = _create(star, dna, lower, upper, k, steps, start, stop)
    
    def write(self, fname='col.vars'):
        with open(fname, 'w') as f:
            f.write(self.string)
