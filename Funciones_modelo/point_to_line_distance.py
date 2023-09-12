import numpy as np

def point_to_line_distance(pt, v1, v2):
    # Calculate distance between a point and a line in 2D or 3D.
    # pt is a nx3 matrix with xyz coordinates for n points
    # v1 and v2 are vertices on the line (each 1x3)
    # d is a nx1 vector with the orthogonal distances

    # Check inputs
    if pt.shape[1] not in [2, 3] or pt.size == 0:
        raise ValueError('First input (pt) is not numeric or has an incorrect shape.')

    if v1.size != pt.shape[1]:
        raise ValueError('Second input (v1) is not numeric or has an incorrect size. Expected 1x3 or 1x2.')

    if v2.size != pt.shape[1]:
        raise ValueError('Third input (v2) is not numeric or has an incorrect size. Expected 1x3 or 1x2.')

    # Prepare inputs
    v1 = v1.flatten()
    v2 = v2.flatten()

    if len(v1) == 2:
        v1 = np.append(v1, 0)

    if len(v2) == 2:
        v2 = np.append(v2, 0)

    if pt.shape[1] == 2:
        pt = np.column_stack((pt, np.zeros(pt.shape[0])))

    v1_ = np.tile(v1, (pt.shape[0], 1))
    v2_ = np.tile(v2, (pt.shape[0], 1))

    # Actual calculation
    a = v1_ - v2_
    b = pt - v2_
    distance = np.sqrt(np.sum(np.cross(a, b, axis=1) ** 2, axis=1)) / np.sqrt(np.sum(a ** 2, axis=1))

    return distance