import numpy as np
from math import cos, sin, acos, asin, atan, pi


def get_dcm(angles, group, rad=False):
    """
    Creation of a DCM from a vector of Euler angles according to a group.
    :param angles: Euler angles vector of shape (3, 1). Can be in radians or degrees according to the rad flag.
    :param group: String indicating the group used for the Euler angles. (3 digits).
    :param rad: (optional) If true the angles must be given in radians, otherwise in degrees, default is false.
    :return: DCM of shape (3, 3).
    """
    if not (angles.shape == (3,1)):
        print("Error. Euler angles vector must be of shape (3, 1)")
        return -1
    if not rad:
        angles = angles / (180 / pi)
    if group == "121":
        return c121(angles)
    elif group == "123":
        return c123(angles)
    elif group == "131":
        return c131(angles)
    elif group == "132":
        return c132(angles)
    elif group == "212":
        return c212(angles)
    elif group == "213":
        return c213(angles)
    elif group == "231":
        return c231(angles)
    elif group == "232":
        return c232(angles)
    elif group == "312":
        return c312(angles)
    elif group == "313":
        return c313(angles)
    elif group == "321":
        return c321(angles)
    elif group == "323":
        return c323(angles)
    else:
        print("Error. There is no", group, "group")
        return -1


def get_euler_angles(dcm, group, rad=False):
    """
    Creation of a Euler angles vector from a DDCM according to a group.
    :param dcm: DCM of shape (3, 3).
    :param group: String indicating the group used for the Euler angles. (3 digits).
    :param rad: (optional) If true the angles must be given in radians, otherwise in degrees, default is false.
    :return: Euler angles vector of shape (3, 1). Can be in radians or degrees according to the rad flag.
    """
    if not (dcm.shape == (3,3)):
        print("Error. DCM's must be of shape (3, 3)")
        return -1
    if group == "121":
        angles = angles121(dcm)
    elif group == "123":
        angles = angles123(dcm)
    elif group == "131":
        angles = angles131(dcm)
    elif group == "132":
        angles = angles132(dcm)
    elif group == "212":
        angles = angles212(dcm)
    elif group == "213":
        angles = angles213(dcm)
    elif group == "231":
        angles = angles231(dcm)
    elif group == "232":
        angles = angles232(dcm)
    elif group == "312":
        angles = angles312(dcm)
    elif group == "313":
        angles = angles313(dcm)
    elif group == "321":
        angles = angles321(dcm)
    elif group == "323":
        angles = angles323(dcm)
    else:
        print("Error. There is no", group, "group")
        return -1
    if not rad:
        angles = angles * 180 / pi
    return angles


def get_rates(angles, w, group, rad=False):
    if not (angles.shape == (3,1)):
        print("Error. Euler angles vector must be of shape (3, 1)")
        return -1
    if not (w.shape == (3,1)):
        print("Error. Angular rates' vector must be of shape (3, 1)")
        return -1
    if not rad:
        angles = angles / (180 / pi)
    if group == "121":
        return np.dot(rates_121(angles), w)
    elif group == "123":
        return np.dot(rates_123(angles), w)
    elif group == "131":
        return np.dot(rates_131(angles), w)
    elif group == "132":
        return np.dot(rates_132(angles), w)
    elif group == "212":
        return np.dot(rates_212(angles), w)
    elif group == "213":
        return np.dot(rates_213(angles), w)
    elif group == "231":
        return np.dot(rates_231(angles), w)
    elif group == "232":
        return np.dot(rates_232(angles), w)
    elif group == "312":
        return np.dot(rates_312(angles), w)
    elif group == "313":
        return np.dot(rates_313(angles), w)
    elif group == "321":
        return np.dot(rates_321(angles), w)
    elif group == "323":
        return np.dot(rates_323(angles), w)
    else:
        print("Error. There is no", group, "group")
        return -1


def add_rotations(rot1, rot2, group1, group2=None, rad=False):
    if not (rot1.shape == rot2.shape == (3,1)):
        print("Error. Euler angles vector must be of shape (3, 1)")
        return -1
    if not rad:
        rot1 = rot1 / (180 / pi)
        rot2 = rot2 / (180 / pi)
    g = list(group1)
    if (group2 is None or group1 == group2) and g[0] == g[2]:
        a2 = acos(cos(rot1[1]) * cos(rot2[1]) - sin(rot1[1]) * sin(rot2[1]) * cos(rot1[2] + rot2[0]))
        angles = np.array([[rot1[0] + arctan((sin(rot1[1]) * sin(rot2[1]) * sin(rot1[2] + rot2[0])),
                                             (cos(rot2[1] - cos(rot1[1]) * cos(a2))))],
                           [a2],
                           [rot2[2] + arctan((sin(rot1[1]) * sin(rot2[1]) * sin(rot1[2] + rot2[0])),
                                             (cos(rot1[1] - cos(rot2[1]) * cos(a2))))]])
        if not rad:
            angles = angles * (180 / pi)
    else:
        if group2 is None:
            group2 = group1
        dcm1 = get_dcm(rot1, group1)
        dcm2 = get_dcm(rot2, group2)
        dcm = np.dot(dcm2, dcm1)
        angles = get_euler_angles(dcm, group1, rad)
    return angles


def c121(angles):
    return np.array([
        [
            cos(angles[1]),
            sin(angles[1]) * sin(angles[0]),
            -sin(angles[1]) * cos(angles[0])
        ],
        [
            sin(angles[2]) * sin(angles[1]),
            -sin(angles[2]) * cos(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * cos(angles[1]) * cos(angles[0]) + cos(angles[2]) * sin(angles[0])
        ],
        [
            cos(angles[2]) * sin(angles[1]),
            -cos(angles[2]) * cos(angles[1]) * sin(angles[0]) - sin(angles[2]) * cos(angles[0]),
            cos(angles[2]) * cos(angles[1]) * cos(angles[0]) - sin(angles[2]) * sin(angles[0])
        ]])


def c123(angles):
    return np.array([
        [
            cos(angles[2]) * cos(angles[1]),
            cos(angles[2]) * sin(angles[1]) * sin(angles[0]) + sin(angles[2]) * cos(angles[0]),
            -cos(angles[2]) * sin(angles[1]) * cos(angles[0]) + sin(angles[2]) * sin(angles[0])
        ],
        [
            -sin(angles[2]) * cos(angles[1]),
            -sin(angles[2]) * sin(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * sin(angles[1]) * cos(angles[0]) + cos(angles[2]) * sin(angles[0])
        ],
        [
            sin(angles[1]),
            -cos(angles[1]) * sin(angles[0]),
            cos(angles[1]) * cos(angles[0])
        ]])


def c131(angles):
    return np.array([
        [
            cos(angles[1]),
            sin(angles[1]) * cos(angles[0]),
            sin(angles[1]) * sin(angles[0])
        ],
        [
            -cos(angles[2]) * sin(angles[1]),
            cos(angles[2]) * cos(angles[1]) * cos(angles[0]) - sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * cos(angles[1]) * sin(angles[0]) + sin(angles[2]) * cos(angles[0])
        ],
        [
            sin(angles[2]) * sin(angles[1]),
            -sin(angles[2]) * cos(angles[1]) * cos(angles[0]) - cos(angles[2]) * sin(angles[0]),
            -sin(angles[2]) * cos(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0])
        ]])


def c132(angles):
    return np.array([
        [
            cos(angles[2]) * cos(angles[1]),
            cos(angles[2]) * sin(angles[1]) * cos(angles[0]) + sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * sin(angles[1]) * sin(angles[0]) - sin(angles[2]) * cos(angles[0])
        ],
        [
            -sin(angles[1]),
            cos(angles[1]) * cos(angles[0]),
            cos(angles[1]) * sin(angles[0])
        ],
        [
            sin(angles[2]) * cos(angles[1]),
            sin(angles[2]) * sin(angles[1]) * cos(angles[0]) - cos(angles[2]) * sin(angles[0]),
            sin(angles[2]) * sin(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0])
        ]])


def c212(angles):
    return np.array([
        [
            -sin(angles[2]) * cos(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * sin(angles[1]),
            -sin(angles[2]) * cos(angles[1]) * cos(angles[0]) - cos(angles[2]) * sin(angles[0])
        ],
        [
            sin(angles[1]) * sin(angles[0]),
            cos(angles[1]),
            sin(angles[1]) * cos(angles[0])
        ],
        [
            cos(angles[2]) * cos(angles[1]) * sin(angles[0]) + sin(angles[2]) * cos(angles[0]),
            -cos(angles[2]) * sin(angles[1]),
            cos(angles[2]) * cos(angles[1]) * cos(angles[0]) - sin(angles[2]) * sin(angles[0])
        ]])


def c213(angles):
    return np.array([
        [
            sin(angles[2]) * sin(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * cos(angles[1]),
            sin(angles[2]) * sin(angles[1]) * cos(angles[0]) - cos(angles[2]) * sin(angles[0])
        ],
        [
            cos(angles[2]) * sin(angles[1]) * sin(angles[0]) - sin(angles[2]) * cos(angles[0]),
            cos(angles[2]) * cos(angles[1]),
            cos(angles[2]) * sin(angles[1]) * cos(angles[0]) + sin(angles[2]) * sin(angles[0])
        ],
        [
            cos(angles[1]) * cos(angles[0]),
            -sin(angles[1]),
            cos(angles[1]) * cos(angles[0])
        ]])


def c231(angles):
    return np.array([
        [
            cos(angles[1]) * cos(angles[0]),
            sin(angles[1]),
            -cos(angles[1]) * sin(angles[0])
        ],
        [
            -cos(angles[2]) * sin(angles[1]) * cos(angles[0]) + sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * cos(angles[1]),
            cos(angles[2]) * sin(angles[1]) * sin(angles[0]) + sin(angles[2]) * cos(angles[0])
        ],
        [
            sin(angles[2]) * sin(angles[1]) * cos(angles[0]) + cos(angles[2]) * sin(angles[0]),
            -sin(angles[2]) * cos(angles[1]),
            -sin(angles[2]) * sin(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0])
        ]])


def c232(angles):
    return np.array([
        [
            cos(angles[2]) * cos(angles[1]) * cos(angles[0]) - sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * sin(angles[1]),
            -cos(angles[2]) * cos(angles[1]) * sin(angles[0]) - sin(angles[2]) * cos(angles[0])
        ],
        [
            -sin(angles[1]) * cos(angles[0]),
            cos(angles[1]),
            sin(angles[1]) * sin(angles[0])
        ],
        [
            sin(angles[2]) * cos(angles[1]) * cos(angles[0]) + cos(angles[2]) * sin(angles[0]),
            sin(angles[2]) * sin(angles[1]),
            -sin(angles[2]) * cos(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0])
        ]])


def c312(angles):
    return np.array([
        [
            -sin(angles[2]) * sin(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * sin(angles[1]) * cos(angles[0]) + cos(angles[2]) * sin(angles[0]),
            -sin(angles[2]) * cos(angles[1])
        ],
        [
            -cos(angles[1]) * sin(angles[0]),
            cos(angles[1]) * cos(angles[0]),
            sin(angles[1])
        ],
        [
            cos(angles[2]) * sin(angles[1]) * sin(angles[0]) + sin(angles[2]) * cos(angles[0]),
            -cos(angles[2]) * sin(angles[1]) * cos(angles[0]) + sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * cos(angles[1])
        ]])


def c313(angles):
    return np.array([
        [
            -sin(angles[2]) * cos(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * cos(angles[1]) * cos(angles[0]) + cos(angles[2]) * sin(angles[0]),
            sin(angles[2]) * sin(angles[1])
        ],
        [
            -cos(angles[2]) * cos(angles[1]) * sin(angles[0]) - sin(angles[2]) * cos(angles[0]),
            cos(angles[2]) * cos(angles[1]) * cos(angles[0]) - sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * sin(angles[1])
        ],
        [
            sin(angles[1]) * sin(angles[0]),
            -sin(angles[1]) * cos(angles[0]),
            cos(angles[1])
        ]])


def c321(angles):
    return np.array([
        [
            cos(angles[1]) * cos(angles[0]),
            cos(angles[1]) * sin(angles[0]),
            -sin(angles[1])
        ],
        [
            sin(angles[2]) * sin(angles[1]) * cos(angles[0]) - cos(angles[2]) * sin(angles[0]),
            sin(angles[2]) * sin(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * cos(angles[1])
        ],
        [
            cos(angles[2]) * sin(angles[1]) * cos(angles[0]) + sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * sin(angles[1]) * sin(angles[0]) - sin(angles[2]) * cos(angles[0]),
            cos(angles[2]) * cos(angles[1])
        ]])


def c323(angles):
    return np.array([
        [
            cos(angles[2]) * cos(angles[1]) * cos(angles[0]) - sin(angles[2]) * sin(angles[0]),
            cos(angles[2]) * cos(angles[1]) * sin(angles[0]) + sin(angles[2]) * cos(angles[0]),
            -cos(angles[2]) * sin(angles[1])
        ],
        [
            -sin(angles[2]) * cos(angles[1]) * cos(angles[0]) - cos(angles[2]) * sin(angles[0]),
            -sin(angles[2]) * cos(angles[1]) * sin(angles[0]) + cos(angles[2]) * cos(angles[0]),
            sin(angles[2]) * sin(angles[1])
        ],
        [
            sin(angles[1]) * cos(angles[0]),
            sin(angles[1]) * sin(angles[0]),
            cos(angles[1])
        ]])


def arctan(num, den):
    r = atan(num / den)
    if num >= 0 and den >= 0:
        return r
    elif num >= 0 > den:
        return pi + r
    elif num < 0 and den < 0:
        return -pi + r
    else:
        return r


def angles121(dcm):
    angles = np.array([[arctan(dcm[0][1], -dcm[0][2])],
                       [acos(dcm[0][0])],
                       [arctan(dcm[1][0], dcm[2][0])]])
    return angles


def angles123(dcm):
    angles = np.array([[arctan(-dcm[2][1], dcm[2][2])],
                       [asin(dcm[2][0])],
                       [arctan(-dcm[1][0], dcm[0][0])]])
    return angles


def angles131(dcm):
    angles = np.array([[arctan(dcm[0][2], dcm[0][1])],
                       [acos(dcm[0][0])],
                       [arctan(dcm[2][0], -dcm[1][0])]])
    return angles


def angles132(dcm):
    angles = np.array([[arctan(dcm[1][2], dcm[1][1])],
                       [asin(-dcm[0][0])],
                       [arctan(dcm[2][0], dcm[0][0])]])
    return angles


def angles212(dcm):
    angles = np.array([[arctan(dcm[1][0], dcm[1][2])],
                       [acos(dcm[1][1])],
                       [arctan(dcm[0][1], -dcm[2][1])]])
    return angles


def angles213(dcm):
    angles = np.array([[arctan(dcm[2][0], dcm[2][2])],
                       [-asin(dcm[2][1])],
                       [arctan(dcm[0][1], dcm[1][1])]])
    return angles


def angles231(dcm):
    angles = np.array([[arctan(-dcm[0][2], dcm[0][0])],
                       [asin(dcm[0][1])],
                       [arctan(-dcm[2][1], dcm[1][1])]])
    return angles


def angles232(dcm):
    angles = np.array([[arctan(dcm[1][2], -dcm[1][0])],
                       [acos(dcm[1][1])],
                       [arctan(dcm[2][1], dcm[0][1])]])
    return angles


def angles312(dcm):
    angles = np.array([[arctan(dcm[1][1], -dcm[1][0])],
                       [asin(dcm[1][2])],
                       [arctan(-dcm[0][2], dcm[2][2])]])
    return angles


def angles313(dcm):
    angles = np.array([[arctan(dcm[2][0], -dcm[2][1])],
                       [acos(dcm[2][2])],
                       [arctan(dcm[0][2], dcm[1][2])]])
    return angles


def angles321(dcm):
    angles = np.array([[arctan(dcm[0][1], dcm[0][0])],
                       [-asin(dcm[0][2])],
                       [arctan(dcm[1][2], dcm[2][2])]])
    return angles


def angles323(dcm):
    angles = np.array([[arctan(dcm[2][1], dcm[2][0])],
                       [acos(dcm[2][2])],
                       [arctan(dcm[1][2], -dcm[0][2])]])
    return angles


def rates_121(angles):
    r = np.array([[0, sin(angles[2]), cos(angles[2])],
                  [0, sin(angles[1]) * cos(angles[2]), - sin(angles[1]) * sin(angles[2])],
                  [sin(angles[1]), -cos(angles[1]) * sin(angles[2]), -cos(angles[1]) * cos(angles[2])]]) / sin(
        angles[1])
    return r


def rates_123(angles):
    r = np.array([[cos(angles[2]), -sin(angles[2]), 0],
                  [cos(angles[1]) * sin(angles[2]), cos(angles[1]) * cos(angles[2]), 0],
                  [-sin(angles[1]) * cos(angles[2]), sin(angles[1]) * sin(angles[2]), cos(angles[1])]]) / cos(angles[1])
    return r


def rates_131(angles):
    r = np.array([[0, -cos(angles[2]), sin(angles[2])],
                  [0, sin(angles[1]) * sin(angles[2]), sin(angles[1]) * cos(angles[2])],
                  [sin(angles[1]), cos(angles[1]) * cos(angles[2]), - cos(angles[1]) * sin(angles[2])]]) / sin(
        angles[1])
    return r


def rates_132(angles):
    r = np.array([[cos(angles[2]), 0, sin(angles[2])],
                  [-cos(angles[1]) * sin(angles[2]), 0, cos(angles[1]) * cos(angles[2])],
                  [sin(angles[1]) * cos(angles[2]), cos(angles[1]), sin(angles[1]) * sin(angles[2])]]) / cos(angles[1])
    return r


def rates_212(angles):
    r = np.array([[sin(angles[2]), 0 - cos(angles[2])],
                  [sin(angles[1]) * cos(angles[2]), 0, sin(angles[1]) * sin(angles[2])],
                  [-cos(angles[1]) * sin(angles[2]), sin(angles[1]), cos(angles[1]) * cos(angles[2])]]) / sin(angles[1])
    return r


def rates_213(angles):
    r = np.array([[sin(angles[2]), cos(angles[2]), 0],
                  [cos(angles[1]) * cos(angles[2]), -cos(angles[1]) * sin(angles[2]), 0],
                  [sin(angles[1]) * sin(angles[2]), sin(angles[1]) * cos(angles[2]), cos(angles[1])]]) / cos(angles[1])
    return r


def rates_231(angles):
    r = np.array([[0, cos(angles[2]), -sin(angles[2])],
                  [0, cos(angles[1]) * sin(angles[2]), cos(angles[1]) * cos(angles[2])],
                  [cos(angles[1]), -sin(angles[1]) * cos(angles[2]), sin(angles[1]) * sin(angles[2])]]) / cos(angles[1])
    return r


def rates_232(angles):
    r = np.array([[cos(angles[2]), 0, sin(angles[2])],
                  [-sin(angles[1]) * sin(angles[2]), 0, sin(angles[1]) * cos(angles[2])],
                  [-cos(angles[1]) * cos(angles[2]), sin(angles[1]), -cos(angles[1]) * sin(angles[2])]]) / sin(
        angles[1])
    return r


def rates_312(angles):
    r = np.array([[-sin(angles[2]), 0, cos(angles[2])],
                  [cos(angles[1]) * cos(angles[2]), 0, cos(angles[1]) * sin(angles[2])],
                  [sin(angles[1]) * sin(angles[2]), cos(angles[1]), -sin(angles[1]) * cos(angles[2])]]) / cos(angles[1])
    return r


def rates_313(angles):
    r = np.array([[sin(angles[2]), cos(angles[2]), 0],
                  [sin(angles[1]) * cos(angles[2]), -sin(angles[1]) * sin(angles[2]), 0],
                  [-cos(angles[1]) * sin(angles[2]), -cos(angles[1]) * cos(angles[2]), sin(angles[1])]]) / sin(
        angles[1])
    return r


def rates_321(angles):
    r = np.array([[0.0, sin(angles[2]), cos(angles[2])],
                  [0.0, cos(angles[1]) * cos(angles[2]), -cos(angles[1]) * sin(angles[2])],
                  [cos(angles[1]), sin(angles[1]) * sin(angles[2]), sin(angles[1]) * cos(angles[2])]]) / cos(angles[1])
    print(cos(angles[1]))
    return r


def rates_323(angles):
    r = np.array([[-cos(angles[2]), sin(angles[2]), 0],
                  [sin(angles[1]) * sin(angles[2]), sin(angles[1]) * cos(angles[2]), 0],
                  [cos(angles[1]) * cos(angles[2]), -cos(angles[1]) * sin(angles[2]), sin(angles[1])]]) / sin(angles[1])
    return r
