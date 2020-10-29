import numpy as np

def static_attitude(a_i, v_i, s_i):
    """Short summary.

    Parameters
    ----------
    a_i : type
        Description of parameter `a_i`.
    v_i : type
        Description of parameter `v_i`.
    s_i : type
        Description of parameter `s_i`.

    Returns
    -------
    type
        Description of returned object.

    """

    B = np.dot(np.multiply(a_i, s_i), v_i.T)
    S = B + B.T
    sigma = np.trace(B)
    Z = np.array([[B[1][2]-B[2][1]],
              [B[2][0]-B[0][2]],
              [B[0][1]-B[1][0]]])
    K1 = np.append(np.array([[sigma]]),Z.T,axis=1)
    K2 = np.append(Z,(S - np.identity(3)*sigma), axis=1)
    K=np.append(K1, K2, axis=0)
    ## Check K
    A = np.sum(a_i)
    error = 1
    while error > 1e-10:
        A_old = A
        r = np.linalg.det(K - A * np.identity(4))
        r1 = np.linalg.det(K - A * np.identity(4)) * np.trace(np.linalg.inv(K - A * np.identity(4)) * - np.identity(4))
        A = A - r/r1
        error = abs(A - A_old)
    q_a = ((A+sigma) * np.identity(3) - S)
    q_a = np.linalg.inv(q_a)
    q = np.dot(q_a ,Z)
    print(q)
    print(get_dcm(q))


# TESTING CODE
# v_i = np.array([[1,0],
#                 [0,0],
#                 [0,1]])
# s_i = np.array([[0.0099994,0.0049997],
#                 [0.99994,-0.0049997],
#                 [0.0049997,0.99994]])
#
# weight = np.array([0.5, 0.5])
# static_attitude(weight, v_i, s_i)
