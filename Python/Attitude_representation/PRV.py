import numpy as np
from math import cos, sin, acos, pi


def explode_prv(prv):
    return prv["phi"], prv["e"]


def condense_prv(phi, e):
    prv = {"phi": phi,
           "e": e}
    return prv


def get_dcm(prv, rad=False):
    phi, e = explode_prv(prv)
    if not (e.shape == (3, 1)):
        print("Error. Principal axis vector must be of shape (3, 1)")
        return -1
    if not rad:
        phi = phi / (180 / pi)
    s = 1 - cos(phi)
    c = np.array([[e[0] ** 2 * s + cos(phi), e[0] * e[1] * s + e[2] * sin(phi), e[0] * e[2] * s - e * sin(phi)],
                  [e[1] * e[0] * s - e[2] * sin(phi), e[1] ** 2 * s + cos(phi), e[1] * e[2] * s + e[0] * sin(phi)],
                  [e[2] * e[0] * s + e[1] * sin(phi), e[2] * e[1] * s - e[0] * sin(phi), e[2] ** 2 * s + cos(phi)]])
    return c


def get_prv(c, rad=False):
    if not (c.shape == (3, 3)):
        print("Error. DCM matrix must be of shape (3, 3)")
        return -1
    phi = acos(0.5 * (c[0][0] + c[1][1] + c[2][2] - 1))
    e = (1 / (2 * sin(phi))) * np.array([[c[1][2] - c[2][1]], [c[2][0] - c[0][2]], [c[0][1] - c[1][0]]])
    if not rad:
        phi = phi * (180 / pi)
    prv = condense_prv(phi, e)
    return prv


def add_prv(prv1, prv2, rad=False):
    phi1, e1 = explode_prv(prv1)
    phi2, e2 = explode_prv(prv2)

    if not (e1.shape == e2.shape == (3, 1)):
        print("Error. Principal axis vector must be of shape (3, 1)")
        return -1
    if not rad:
        phi1 = phi1 / (180 / pi)
        phi2 = phi2 / (180 / pi)
    c1 = get_dcm(phi1, e1)
    c2 = get_dcm(phi2, e2)
    c = np.dot(c2, c1)
    prv = get_prv(c)
    return prv
