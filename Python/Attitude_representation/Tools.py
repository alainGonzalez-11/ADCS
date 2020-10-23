import numpy as np

def tilde(x):
    """
    Tilde operator for a vector defined by:
        |  0   -x3   x2  |
        |  x3   0   -x1  |
        | -x2   x1   0   |
    :param x: Vector of size (3, 1)
    :return: Tilde matrix of shape (3, 3)
    """
    x = np.array([[0, -x[2][0], x[1][0]],
                  [x[2][0], 0, -x[0][0]],
                  [-x[1][0], x[0][0], 0]])
    return x

def normalize(x):
    """
    Normalizes a vector to a unit vector.
    :param x: Original vector
    :return: Normalized vector
    """
    x = x / np.linalg.norm(x)
    return x

def cross(x, y):
    """
    Cross product of vectors x and y (x X y).
    :param x: First vector of shape (3, 1)
    :param y: Second vector of shape (3, 1)
    :return: Cross product vector
    """
    x=x.reshape(3)
    y=y.reshape(3)
    z = np.cross(x,y)
    z =z.reshape((3,1))
    return z
