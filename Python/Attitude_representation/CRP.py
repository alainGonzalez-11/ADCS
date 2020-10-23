import numpy as np
import Quaternions


def get_dcm(q):
    if not (q.shape == (3, 1)):
        print("Error. Classical Rodrigues Parameters' vector must be of shape (3, 1)")
        return -1
    c = np.array([[1 + q[0][0] ** 2 - q[1][0] ** 2 - q[2][0] ** 2, 2 * (q[0][0] * q[1][0] + q[2][0]), 2 * (q[0][0] * q[2][0] - q[1][0])],
                  [2 * (q[1][0] * q[0][0] - q[2][0]), 1 - q[0][0] ** 2 + q[1][0] ** 2 - q[2][0] ** 2, 2 * (q[1][0] * q[2][0] + q[0][0])],
                  [2 * (q[2][0] * q[0][0] + q[1][0]), 2 * (q[2][0] * q[1][0] - q[0][0]), 1 - q[0][0] ** 2 - q[1][0] ** 2 + q[2][0] ** 2]]) / (
                    1 + np.dot(q.T, q))
    return c
def get_rates(q, w):
    if not (q.shape == (3, 1)):
        print("Error. Classical Rodrigues Parameters' vector must be of shape (3, 1)")
        return -1
    if not (q.shape == (3, 1)):
        print("Error. Angular rates' vector must be of shape (3, 1)")
    q = q.reshape(3)
    m = np.array([[1 + q[0][0] ** 2, q[0][0] * q[1][0] - q[2][0], q[0][0] * q[2][0] + q[1][0]],
                  [q[1][0] * q[0][0] + q[2][0], 1 + q[1][0] ** 2, q[1][0] * q[2][0] - q[0][0]],
                  [q[2][0] * q[0][0] - q[1][0], q[2][0] * q[1][0] + q[0][0], 1 + q[2][0] ** 2]])
    dq = np.dot(m, w) / 2
    return dq

def get_crp(c):
    if not (c.shape == (3, 3)):
        print("Error. DCM matrix must be of shape (3, 3)")
        return -1
    b = Quaternions.get_quaternions(c)
    q = np.zeros((3,1))
    for i in range(q.size):
        q[i][0] = b[i + 1][0] / b[0][0]
    print(q)
    return q


def addition(q1, q2):
    if not (q1.shape == q2.shape == (3, 1)):
        print("Error. Classical Rodrigues Parameters' vector must be of shape (3, 1)")
        return -1
    q = (q2 + q1 - np.cross(q2, q1)) / (1 - np.dot(q2, q1))
    return q


def substraction(q, q1):
    if not (q.shape == q1.shape == (3, 1)):
        print("Error. Classical Rodrigues Parameters' vector must be of shape (3, 1)")
        return -1
    q2 = (q - q1 + np.cross(q, q1)) / (1 - np.dot(q, q1))
    return q2



