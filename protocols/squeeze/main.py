#!/usr/bin/env python
"""Example script for generating a single topology within a geometry
with a given shape and size e.g. a sphere with radius=5, a cube with
side=5 etc.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from numpy.core.fromnumeric import squeeze
from starpolymers.molecules import (
    MoleculeFactory,
    StarPolyelectrolyte,
    LinearPolyelectrolyte,
    AbstractMolecule
)

from itertools import product, combinations

from softnanotools.logger import Logger
logger = Logger('MAIN')

class SqueezeGeometry:
    __valid_geometries__ = {
        'cube': {'L'},
        'sphere': {'R'}
    }

    def __init__(
        self,
        kind: str = 'cube',
        molecule: AbstractMolecule = None,
        spacing: float = 1.0,
        **kwargs
    ):
        assert kind in SqueezeGeometry.__valid_geometries__
        #assert set(kwargs.keys()) == SqueezeGeometry.__valid_geometries__.values()
        self.kind = kind
        self.spacing = spacing

        self._molecule = molecule
        self._bound = 0.0

    @property
    def bound(self):
        return self._bound

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, obj: AbstractMolecule):
        assert isinstance(obj, AbstractMolecule)
        self._molecule = obj

    def initialise(self, obj: AbstractMolecule):
        self.molecule = obj
        return

    def build(
        self,
        strategy: str = 'zig-zag'
    ):
        """Takes a Molecule's parameters and rebuilds it to ensure that it is
        centered (by its COM) in the centre of the geometry, and that
        all of the atoms are within the bounds of the geometry"""

        if strategy == 'zig-zag':
            spacing = self.spacing
            # we need to calculate the min and max angle

            # the min angle is the minimum angle needed to avoid the
            # next-but-one particle in a sequence
            maximum_angle = np.deg2rad(60.0)

            # the max angle is a maximum angle needed to avoid the edge
            # of the box

            # the min angle is given by the value of ((0.5 * lambda) / box)
            # times the spacing, the ratio gives the
            minimum_angle = float()

            # vertical spacing = hypotenuse
            # the other two sides = spacing

            # therefore the angle = cosine of (spacing / (0.5 * vertical))
            vertical = 2.0 * min([self.bound / self.molecule.lam, spacing])
            print(vertical)
            minimum_angle = np.arccos(spacing / 0.5 * vertical)

            print(np.degrees(maximum_angle), np.degrees(minimum_angle))

            assert maximum_angle >= minimum_angle, \
                f"{maximum_angle} not greater than {minimum_angle}"

            deviation = minimum_angle

        else:
            raise NotImplementedError

        print(f"Deviation Angle is {deviation}")

        return deviation

class SqueezeCube(SqueezeGeometry):
    def __init__(self, length: float):
        super().__init__('cube', L=length)
        self._bound = length

    def add_to_ax(self, ax: plt.Axes):
        # draw cube
        r = np.array([-1, 1]) * self._bound
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s-e)) == r[1]-r[0]:
                ax.plot3D(*zip(s, e), color="b", alpha=0.2)

class SqueezeSphere(SqueezeGeometry):
    def __init__(self, radius: float):
        super().__init__('sphere', R=radius)
        self._bound = radius

    def add_to_ax(self, ax: plt.Axes):
        # draw sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = self._bound * np.cos(u)*np.sin(v)
        y = self._bound * np.sin(u)*np.sin(v)
        z = self._bound * np.cos(v)
        ax.plot_wireframe(x, y, z, color="r", alpha=0.2)

class SqueezedStarPolyelectrolyte(StarPolyelectrolyte):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def generate_atoms(self, squeeze: SqueezeGeometry = None):
        if squeeze:
            squeeze.initialise(self)
            squeeze.build()

        spacing = self.settings['spacing']
        direction = self.settings['direction']
        mol_length = self.lam * spacing
        x = []
        y = []
        z = []
        for i in range(self.kap):
            for j in range(self.lam):
                x.append((mol_length-(j*spacing))*direction[i][0])
                y.append((mol_length-(j*spacing))*direction[i][1])
                z.append((mol_length-(j*spacing))*direction[i][2])

        x.append(0.0)
        y.append(0.0)
        z.append(0.0)

        data = pd.DataFrame({
                'mol' : len(x) * [self.mol],
                'type' : len(x) * [self.types['atom']],
                'x' : x,
                'y' : y,
                'z' : z,
                'q' : self.generate_charges(x)
            })

        return data

star = {
        'molecule': 'star',
        'kap': 4,
        'lam': 10,
        'charge' : {
            'max': 1,
            'style': 'all'
        },
        'central' : 'all',
        'counterions': False,
        'angle_type': 1
    }

def visualise(molecule: AbstractMolecule, squeeze: SqueezeGeometry):

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect("auto")
    squeeze.add_to_ax(ax)
    positions = molecule.generate_atoms(squeeze=squeeze)[['x', 'y', 'z']]
    positions = positions.to_numpy()
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], 'o')
    fig.savefig('squeeze.png')
    pass


def main(*args, **kwargs):
    box = 20
    molecule = SqueezedStarPolyelectrolyte(star)
    squeeze = SqueezeSphere(8.0)
    #print(molecule._atoms)
    visualise(molecule, squeeze)
    #print(molecule.generate_atoms())
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    main(**vars(parser.parse_args()))
