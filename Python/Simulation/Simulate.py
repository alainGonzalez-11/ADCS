import numpy as np
import matplotlib.pyplot as plt
import controlpy
from math import sqrt, pi, tan, sin, cos
from scipy import integrate
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

B, b = B_matrix(0)



t = 3000

dt = .1

N = int(t / dt) + 1
x_arch = []
x_arch1 = []
x_arch2 = []
t_arch = []
x = np.array([[0.05],
              [0.05],
              [0.05],
              [0.1],
              [0.1],
              [0.1],
              [0],
              [0],
              [0]])
q1 = .01
q2 = .01
q3 = 1

r1 = 1000
r2 = 1000
Q = np.concatenate((np.concatenate((q1 * np.identity(3), np.zeros((3, 6))), axis=1),  np.concatenate((np.zeros((3, 3)), q2 * np.identity(3), np.zeros((3, 3))), axis=1),   np.concatenate((np.zeros((3, 6)), q3 * np.identity(3)), axis=1)), axis=0)
R = np.concatenate((np.concatenate((r1 * np.identity(3), np.zeros((3, 3))), axis=1),  np.concatenate((np.zeros((3, 3)), r2 * np.identity(3)), axis=1)), axis=0)


Q
R

K = controlpy.synthesis.controller_lqr(A, B, Q, R)[0]
K
def get_shadow(mrp):
    s = - mrp / (np.linalg.norm(mrp) ** 2)
    return s

w_RSW = np.array([[0],  [-w0],         [0]])


A
B

for i in range(1, N):
    t = i * dt
    u = -np.dot(K, x)


    dx = np.dot(A, x) + np.dot(B, u)


    x = x + dt*dx
    s = x[0: 3]
    print(np.linalg.norm(s))
    if np.linalg.norm(s) >1:
        s = get_shadow(s)
        x[0:3] = s

    x_arch = x_arch + [x[0][0]]
    x_arch1 = x_arch1 + [x[3][0]]
    x_arch2 = x_arch2 + [x[6][0]]

    t_arch = t_arch + [t]

    if i%50==0:
        print(x)
        plt.plot(t_arch, x_arch)
        plt.show()

        plt.plot(t_arch, x_arch1)
        plt.show()

        plt.plot(t_arch, x_arch2)
        plt.show()
