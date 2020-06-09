from math import pi
e = 1.602e-19
eps0 = 8.85e-12
kb = 1.38e-23

def bjerrum(epsr, T, reduced=False):
    if reduced == True:
        k = 1
        eps = epsr
        top = 1
    else:
        k = kb
        eps = epsr * eps0
        top = e ** 2
    bottom = eps * k * T * 4 * pi
    length = top/bottom
    print(("The Bjerrum length is {}".format(length)))
    return length
