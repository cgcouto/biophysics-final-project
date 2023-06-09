import numpy as np

# These are some helper functions used in multiple files!

# 
# points () : 
def find_bond_angle(points):

    A = points[:3]
    B = points[3:6]
    C = points[6:]

    vector_AB = A - B
    vector_BC = C - B

    angle = np.arccos(np.dot(vector_AB, vector_BC) / (np.linalg.norm(vector_AB)*np.linalg.norm(vector_BC)))

    return angle


def get_positions(m):
    inds = []
    for i in range(m-1):
        j = i + 1
        inds.append(m * i + j - ((i + 2) * (i + 1)) // 2)

    return np.array(inds)