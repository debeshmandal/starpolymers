from .bjerrum import bjerrum
from math import pi

avo = 6.02e23

def I(salts):
    strength = 0
    for salt in salts:
        strength += salt[0] ** 2 * salt[1]
    return strength

def debye(salts, epsr=74.5, T=310, reduced=False):
    length = 8 * pi * bjerrum(epsr, T, reduced=False) * I(salts) * avo * 1e3
    length = length ** -0.5
    print(("The Debye length is {}".format(length)))
    return length
