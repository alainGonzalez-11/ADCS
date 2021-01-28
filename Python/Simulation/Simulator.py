import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


from math import sqrt, pi, tan, sin, cos


def CrossMatrix(x):
    """Tilde operator for a vector defined by:
        |  0   -x3   x2  |
        |  x3   0   -x1  |
        | -x2   x1   0   |
    Parameters
    ----------
    x : ndarray
        Vector of size (3, 1).

    Returns
    -------
    ndarray
        Matrix of shape (3, 3).

    """
    x = np.array([[0, -x[2][0], x[1][0]],
                  [x[2][0], 0, -x[0][0]],
                  [-x[1][0], x[0][0], 0]])
    return x


def cross(x, y):
    """Cross product of vectors x and y (x X y).

    Parameters
    ----------
    x : ndarray
        First vector of size (3, 1).
    y : ndarray
        Second vector of size (3, 1).

    Returns
    -------
    ndarray
        Cross product vector of size (3,1).

    """
    x = x.reshape(3)
    y = y.reshape(3)
    z = np.cross(x, y)
    z = z.reshape((3, 1))
    return z


def get_dcm(sigmas):
    if not (sigmas.shape == (3, 1)):
        print("Error. Modified Rodrigues Parameters' vector must be of shape (3, 1)")
        return -1
    st = tilde(sigmas)
    c = np.identity(3) + (8 * np.dot(st, st) - 4 * (1 - np.dot(sigmas.T, sigmas)) * st) / (
        (1 + np.dot(sigmas.T, sigmas)) ** 2)
    return c


def get_mrp(c):
    if not (c.shape == (3, 3)):
        print("Error. DCM matrix must be of shape (3, 3)")
        return -1
    gi = sqrt(np.trace(c) + 1)
    s = 1 / (gi * (gi + 2)) * np.array([[c[1][2] - c[2][1],
                                         c[2][0] - c[0][2],
                                         c[0][1] - c[1][0]]]).T
    return s


def get_shadow(mrp):
    s = - mrp / (np.linalg.norm(mrp) ** 2)
    return s


# Constant definitions
J = np.array([[4 / 120, 0, 0],
              [0, 4 / 120, 0],
              [0, 0, 8 / 1200]])

Jw = np.array([[2.65258238E-7, 0, 0],
               [0, 2.65258238E-7, 0],
               [0, 0, 2.65258238E-7]])

altitud = 657
mu = 3.986005e14


a = 6371000 + altitud * 1000
T = 2 * pi * sqrt(a**3 / mu)
w0 = 2 * pi / T
w0
im = ((altitud + 2.3850e+04) / (57.3 * (900 - 500) / (99 - 97.4)))
W_RSW = np.array([[0], [-w0], [0]])

CrossMatrix(W_RSW)
A11 =  np.dot(-0.5 *CrossMatrix(W_RSW), (4 * CrossMatrix(W_RSW) + 16 * np.identity(3)))
A12 = .25 * np.identity(3)
A21 = np.array([[16 * w0 ** 2 * (J[2][2] - J[1][1]) / J[0][0], 0, 0],
                [0, 12 * w0 ** 2 * (J[2][2] - J[0][0]) /     J[1][1], 0],
                [0, 0, 4 * w0 ** 2 * (J[0][0] - J[1][1]) / J[2][2]]])
A22 = np.array([[0, 0, w0 * (J[0][0] - J[1][1] + J[2][2]) / J[0][0]],
                [0, 0, 0],
                [-w0 * (J[0][0] - J[1][1] + J[2][2]) / J[2][2], 0, 0]])
A23 = np.array([[0, 0, w0 * Jw[2][2] / J[0][0]],
                [0, 0, 0],
                [-w0 * Jw[0][0] / J[2][2], 0, 0]])
A = np.concatenate((np.concatenate((A11, A12, np.zeros((3, 3))), axis=1), np.concatenate((A21, A22, A23), axis=1), np.concatenate((np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))), axis=1)), axis=0)

def B_matrix(t):
    mf = 7.9e15
    b = (mf / a ** 3 )* np.array([[cos(w0 * t) * sin(im)],
                                [-cos(im)],
                                [2 * sin(w0 * t) * sin(im)]])

    B = np.concatenate((np.zeros((3, 6)),
                        np.concatenate(
                            (np.linalg.inv(J), -np.dot(np.linalg.inv(J), CrossMatrix(b))), axis=1),
                        np.concatenate((-np.linalg.inv(Jw), np.zeros((3, 3))), axis=1)), axis=0)
    return B, b

B_matrix(0)[0]

t = 3000

dt = .1

N = int(t / dt) + 1
x_arch = []
x_arch1 = []
x_arch2 = []
t_arch = []
x = np.array([[0.1],
              [0.1],
              [0.1],
              [0.1],
              [0.1],
              [0.1],
              [0],
              [0],
              [0]])

df = pd.read_excel('k1.xlsx')

K = df.to_numpy()
K
index =0
K[6*index:6*(index+1)][:]
index = 1
K[6*index:6*(index+1)][:]


w_RSW = np.array([[0],  [-w0],         [0]])
for i in range(N):
    t = i * dt
    B, b = B_matrix(t)
    index = int((t - t//T * T) // (T/100))
    u = -np.dot(K[6*index:6*(index+1)][:], x)


    s = x[0: 3]
    w = x[3:6]
    O = x[6:9]

    tw = u[0:3]
    m = u[3:6]


    s = x[0: 3]

    B = np.array(B)
    dx = np.dot(A, x) #+ np.dot(B, u)


    x = x + dt*dx
    s = x[0: 3]
    if np.linalg.norm(s) >1:
        s = get_shadow(s)
        x[0:3] = s

    x_arch = x_arch + [x[0][0]]
    x_arch1 = x_arch1 + [x[3][0]]
    x_arch2 = x_arch2 + [x[6][0]]

    t_arch = t_arch + [t]

    if i%20==0:
        print(x)
        plt.plot(t_arch, x_arch)
        plt.show()

        plt.plot(t_arch, x_arch1)
        plt.show()

        plt.plot(t_arch, x_arch2)
        plt.show()
