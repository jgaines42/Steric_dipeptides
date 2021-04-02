##############################################################################################################################
# dihedral = calcDihedral(AtomArray, position)
#
# Input:
#  AtomArray: 4x1 array of the 4 atom indexes involved in the dihedral angle
#  position: Nx3 array of atomic coordinates
#
# Returns:
#   Dihedral angle in degrees
#
# New approach is from https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
##############################################################################################################################

def calcDihedral(AtomArray, position):
    import numpy as np

    # sub_ij = position[AtomArray[0], :] - position[AtomArray[1], :]
    # sub_jk = position[AtomArray[1], :] - position[AtomArray[2], :]
    # sub_kl = position[AtomArray[2], :] - position[AtomArray[3], :]

    # cross_ijk = np.cross(sub_ij, sub_jk)
    # cross_jkl = np.cross(sub_jk, sub_kl)

    # ijk_jkl = sum(cross_ijk * cross_jkl)

    # cijk_mag2 = np.sum(np.square(cross_ijk))
    # cjkl_mag2 = np.sum(np.square(cross_jkl))

    # kj_mag = np.sqrt(np.sum(np.square(sub_jk)))

    # cijk_cjlk = np.cross(cross_ijk, cross_jkl)
    # cosphi = ijk_jkl / np.sqrt(cijk_mag2 * cjkl_mag2)
    # sinphi = np.sum(cijk_cjlk * sub_jk) / np.sqrt(cijk_mag2 * cjkl_mag2) / kj_mag
    # #print('sin', sinphi, 'cos', cosphi)
    # if (sinphi < 0):
    #     phi = np.arccos(cosphi)
    # else:
    #     phi = -np.arccos(cosphi)
    # DA = phi / np.pi * 180.0
    # print(DA)
    # #return DA

    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = position[AtomArray[0], :]
    p1 = position[AtomArray[1], :]
    p2 = position[AtomArray[2], :]
    p3 = position[AtomArray[3], :]

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    #print(DA, np.degrees(np.arctan2(y, x)) )
    return np.degrees(np.arctan2(y, x))
