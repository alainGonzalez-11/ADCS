import numpy as np
from math import sqrt, pi, tan, sin, cos
import basic_operators as bo
import Quaternions


def get_dcm(sigmas):
    if not (sigmas.shape == (3, 1)):
        print("Error. Modified Rodrigues Parameters' vector must be of shape (3, 1)")
        return -1
    st = bo.tilde(sigmas)
    c = np.identity(3) + (8 * np.dot(st, st) - 4 * (1 - np.dot(sigmas.T, sigmas)) * st) / (
            (1 + np.dot(sigmas.T, sigmas)) ** 2)
    return c


def get_mrp(c, source="DCM", rad=False):
    if source == "DCM":
        if not (c.shape == (3, 3)):
            print("Error. DCM matrix must be of shape (3, 3)")
            return -1
        gi = sqrt(np.trace(c) + 1)
        s = 1 / (gi * (gi + 2)) * np.array([[c[1][2] - c[2][1],
                                            c[2][0] - c[0][2],
                                            c[0][1] - c[1][0]]]).T
    elif source == "Quaternions":
        if not (c.shape == (4, 1)):
            print("Error. Quaternions' vector must be of shape (4, 1)")
            return -1
        s = np.array([c[1] / (1 + c[0]),
                      c[2] / (1 + c[0]),
                      c[3] / (1 + c[0])])
    elif source == "PRV":
        phi = c["phi"]
        e = c["e"]
        if not (e.shape == (3, 1)):
            print("Error. Principal axis vector must be of shape (3, 1)")
            return -1
        if not rad:
            phi = phi / (180 / pi)
        s = tan(phi / 4) * e
    return s


def add(mrp1, mrp2):
    mrp1_n = np.linalg.norm(mrp1)
    mrp2_n = np.linalg.norm(mrp2)
    den = (1 + mrp1_n ** 2 * mrp2_n ** 2 - 2 * np.dot(mrp1, mrp2))
    if den < 0.01:
        mrp1 = get_shadow(mrp1)
    mrp = ((1 - mrp1_n ** 2) * mrp2 + (1 - mrp2_n ** 2) * mrp1 - 2 * np.cross(mrp2, mrp1)) / (
            1 + mrp1_n ** 2 * mrp2_n ** 2 - 2 * np.dot(mrp1, mrp2))
    return mrp


def get_shadow(mrp):
    s = - mrp / (np.linalg.norm(mrp) ** 2)
    return s


def get_rates(mrp, w):
    if not (mrp.shape == (3, 1)):
        print("Error. Modified Rodrigues Parameters' vector must be of shape (3, 1)")
        return -1
    if not (w.shape == (3, 1)):
        print("Error. Angular rates' vector must be of shape (3, 1)")
        print(w.shape)
        return -1
    ds = (1 / 4) * ((1 - np.linalg.norm(mrp)**2) * np.identity(3) + 2 * bo.tilde(mrp) + 2 * np.dot(mrp, mrp.T))
    ds = np.dot(ds, w)
    return ds


def integrator():
    t = 42
    N = [500000]
    for n in N:
        s = np.array([[0.4], [0.2], [-0.1]])
        dt = t / n
        for i in range(n + 1):
            t = dt * i
            w = np.array([[sin(0.1 * t)],
                          [0.01],
                          [cos(0.1 * t)]]) * 20 / (180 / pi)
            ds = get_rates(s, w)
            s = s + ds * dt
            if (np.linalg.norm(s) > 1.1):
                print("Invert", t)
                s = get_shadow(s)

        print(n,np.linalg.norm(s))
        print("#" * 45)


def rotation_Sum(AB, BC):
    """
    AC = [AB][BC]
    :param AB:
    :param BC:
    :return:
    """
    num = (1 - np.linalg.norm(BC) ** 2) * AB + (1 - np.linalg.norm(AB) ** 2) * BC - 2  * bo.cross(AB,BC)
    den = 1 + np.linalg.norm(AB) ** 2 * np.linalg.norm(BC) ** 2 - 2 * np.dot(BC, AB)
    AC = num / den
    return AC

def rotation_Subs(AC, BC):
    """
    AB = AC - BC
    AC = [AB][BC]
    [AB] = [AC][BC]-1
    :param AC:
    :param BC:
    :return:
    """
    AC = get_dcm(AC)
    BC = get_dcm(BC)
    AB = np.dot(AC, BC.T)
    AB = get_mrp(AB)

    return AB