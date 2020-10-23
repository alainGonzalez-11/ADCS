import numpy as np


def get_dcm(beta):
    if not (beta.shape == (4, 1)):
        print("Error. Quaternions' vector must be of shape (4, 1)")
        return -1
    c = np.array([[beta[0][0] ** 2 + beta[1][0] ** 2 - beta[2][0] ** 2 - beta[3][0] ** 2,
                   2 * (beta[1][0] * beta[2][0] + beta[0][0] * beta[3][0]), 
                   2 * (beta[1][0] * beta[3][0] - beta[0][0] * beta[2][0])],
                  [2 * (beta[1][0] * beta[2][0] - beta[0][0] * beta[3][0]),
                   beta[0][0] ** 2 - beta[1][0] ** 2 + beta[2][0] ** 2 - beta[3][0] ** 2,
                   2 * (beta[2][0] * beta[3][0] + beta[0][0] * beta[1][0])],
                  [2 * (beta[1][0] * beta[3][0] + beta[0][0] * beta[2][0]), 
                   2 * (beta[2][0] * beta[3][0] - beta[0][0] * beta[1][0]),
                   beta[0][0] ** 2 - beta[1][0] ** 2 - beta[2][0] ** 2 + beta[3][0] ** 2]])
    return c



def add_quaternions(b1, b2, short_rotation=True):
    if not (b1.shape == b2.shape == (4, 1)):
        print("Error. Quaternions' vector must be of shape (4, 1)")
        return -1
    m = np.array([[b2[0][0], -b2[1][0], -b2[2][0], -b2[3][0]],
                  [b2[1][0], b2[0][0], b2[3][0], -b2[2][0]],
                  [b2[2][0], -b2[3][0], b2[0][0], b2[1][0]],
                  [b2[3][0], b2[2][0], -b2[1][0], b2[0][0]]])
    b = np.dot(m, b1)
    if short_rotation and b[0][0] < 0:
        b = b * -1
    return b


def get_quaternions(c, source="DCM", short_rotation=True):
    if source == "DCM":
        if not (c.shape == (3, 3)):
            print("Error. DCM matrix must be of shape (3, 3)")
            return -1
        b = np.zeros((4,1))
        tc = np.trace(c)
        b0_2 = 0.25 * (1 + tc)
        b1_2 = 0.25 * (1 + 2 * c[0][0] - tc)
        b2_2 = 0.25 * (1 + 2 * c[1][1] - tc)
        b3_2 = 0.25 * (1 + 2 * c[2][2] - tc)
        if b0_2 >= b1_2 and b0_2 >= b2_2 and b0_2 >= b3_2:
            b[0][0] = b0_2 ** 0.5
            b[1][0] = (c[1][2] - c[2][1]) / (4 * b[0])
            b[2][0] = (c[2][0] - c[0][2]) / (4 * b[0])
            b[3][0] = (c[0][1] - c[1][0]) / (4 * b[0])
        elif b1_2 >= b0_2 and b1_2 >= b2_2 and b1_2 >= b3_2:
            b[1] [0]= b1_2 ** 0.5
            b[0][0] = (c[1][2] - c[2][1]) / (4 * b[1])
            b[2][0] = (c[0][1] + c[1][0]) / (4 * b[1])
            b[3][0] = (c[2][0] + c[0][2]) / (4 * b[1])
        elif b2_2 >= b0_2 and b2_2 >= b1_2 and b2_2 >= b3_2:
            b[2][0] = b2_2 ** 0.5
            b[0][0]= (c[2][0] - c[0][2]) / (4 * b[2])
            b[1][0] = (c[0][1] + c[1][0]) / (4 * b[2])
            b[3] [0]= (c[1][2] + c[2][1]) / (4 * b[2])
        else:
            b[3][0] = b3_2 ** 0.5
            b[0][0] = (c[0][1] - c[1][0]) / (4 * b[3])
            b[1][0] = (c[2][0] + c[0][2]) / (4 * b[3])
            b[2][0] = (c[1][2] + c[2][1]) / (4 * b[3])

        if short_rotation and b[0] < 0:
            b = b * -1
    elif source == "MRP":
        if not (c.shape == (3, 1)):
            print("Error. Modified Rodrigues Parameters' vector must be of shape (3, 1)")
            return -1
        b = np.array((4, 1))
        b[0][0] = (1 - np.dot(c, c)) / (1 + np.dot(c, c))
        for i in range(1, 4):
            b[i][0] = 2 * c[i] / (1 + np.dot(c, c))
    return b

def substract_quaternions(b, b1, n, short_rotation=True):
    if not (b.shape == b1.shape == (4, 1)):
        print("Error. Quaternions' vector must be of shape (4, 1)")
        return -1
    if n == 1:
        m = np.array([[b1[0][0], -b1[1][0], -b1[2][0], -b1[3][0]],
                      [b1[1][0], b1[0][0], b1[3][0], -b1[2][0]],
                      [b1[2][0], -b1[3][0], b1[0][0], b1[1][0]],
                      [b1[3][0], b1[2][0], -b1[1][0], b1[0][0]]])
    elif n == 2:
        m = np.array([[b1[0][0], -b1[1][0], -b1[2][0], -b1[3][0]],
                      [b1[1][0], b1[0][0], -b1[3][0], b1[2][0]],
                      [b1[2][0], b1[3][0], b1[0][0], -b1[1][0]],
                      [b1[3][0], -b1[2][0], b1[1][0], b1[0][0]]])
    else:
        return -1
    b_subs = np.dot(m.T, b)
    if short_rotation and b_subs[0] < 0:
        b_subs = b_subs * -1
    return b_subs


def get_rates(b, w):
    if not (b.shape == (4, 1)):
        print("Error. Quaternions' vector must be of shape (4, 1)")
        return -1
    if not (w.shape == (3, 1)):
        print("Error. Angular rates' vector must be of shape (3, 1)")
        return -1
    w = np.concatenate((np.array([[0]]), w))
    m = np.array([[b[0][0], -b[1][0], -b[2][0], -b[3][0]],
                  [b[1][0], b[0][0], -b[3][0], b[2][0]],
                  [b[2][0], b[3][0], b[0][0], -b[1][0]],
                  [b[3][0], -b[2][0], b[1][0], b[0][0]]])
    rate = 0.5 * np.dot(m, w)
    return rate
