# `starpolymers`

`starpolymers` is a Python Library used for generating LAMMPS input files and analysing LAMMPS output files, in dump formats and other formats. Additionally, there are scripts used for interacting with Colvars input/output files. 

The package is used to generate molecules and LAMMPS configuration files that follow the `atom_style full` convention. User-defined molecules can be added if the package is installed manually

## Features

- Create LAMMPS system containing 'molecules'
- Molecules include salt, linear polyelectrolytes, star polymers
- Use a template LAMMPS input file and change the variables
- Generate a Colvars input file for free energy calculations
- Import LAMMPS dump files

### Upcoming

- Analysis on dump files
- Plotting tools
- Custom molecules without manual installation

## Installation

### pip

Star Polymers is available on PyPi and can be installed using:

    pip install starpolymers


### conda

From v1.1.3 it will be possible to install `starpolymers` using anaconda, via conda-forge:

    conda install starpolymers -c conda-forge

### manual installation

To install the package manually run the following commands:

    git clone https://github.com/debeshmandal/starpolymers.git
    cd starpolymers
    pip install . -vv

Alternatively, since the package is written entirely in Python, one can clone the package and append the path to the `$PYTHONPATH` environment variable.

## Usage

A variety of examples can be found in the `./test` folder which contains tests for generating files. Notably the `config` and `input_files` subfolders contain examples of generating files needed to create LAMMPS systems and run simulations.

### Example

    import starpolymers
    lammps_system = starpolymers.generators.System(50)
    molecules = starpolymers.generators.MoleculeFactory([
            {
                'molecule' : 'salt', 
                'concentration' : 100, 
                'anion' : 1, 
                'cation' : 1
            }
        ])
    lammps_system.add_molecules(molecules)
    config_file = starpolymers.generators.ConfigFile(system)
    config_file.write('config.dat')
