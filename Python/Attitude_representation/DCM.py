import numpy as np
import tools


def get_rates(dcm, w):
    """
    Obtains the DCM rates from the angular rates
    :param dcm: Direction cosines matrix. Numpy array of size (3, 3)
    :param w: Angular rates (Omega). Numpy array of size (3, 1)
    :return: DCM rates
    """
    w = tools.tilde(w)
    rates = -np.dot(w, dcm)
    return rates


def add(dcm1, dcm2, dcm3 = np.identity(3)):
    """
    DCMs addition
    :param dcm1: Numpy array of size (3, 3)
    :param dcm2: Numpy array of size (3, 3)
    :param dcm2: (optional) Numpy array of size (3, 3)
    :return: Final DCM
    """
    if not (dcm1.shape == dcm2.shape == dcm3.shape == (3,3)):
        raise ValueError("DCM's must be of shape (3, 3)")
    r = np.dot(dcm3, dcm2)
    r = np.dot(r, dcm1)
    return r
