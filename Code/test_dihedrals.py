import numpy as np
from calcDihedral import calcDihedral
from rotate_DA import rotate_DA
A = np.array([[0, 1.0, 0], [0, 0, 0], [1.0, 0, 0], [1.0, 1.0, 0.0]])
B = np.array([0, 1, 2, 3])

X = calcDihedral(B, A)
print(0, X)

InitChi1 = X
setChi = -135
delta_term = np.pi * np.sign(InitChi1) * InitChi1 / 180
moveAtomID = np.array([3])
iChiArray = B
newA = rotate_DA(A, setChi, delta_term, iChiArray, moveAtomID)
print(newA)

X = calcDihedral(B, newA)
print(setChi, X)


setChi = 10
delta_term = np.pi * np.sign(InitChi1) * InitChi1 / 180
moveAtomID = np.array([3])
iChiArray = B
newA = rotate_DA(A, setChi, delta_term, iChiArray, moveAtomID)
print(newA)

X = calcDihedral(B, newA)
print(setChi, X)

# A = np.array([[0, 1.0, 0], [0, 0, 0], [1.0, 0, 0], [1.0, -1.0, -1.0]])
# B = np.array([0, 1, 2, 3])

# X = calcDihedral(B, A)
# print(-135, X)

# A = np.array([[0, 1.0, 0], [0, 0, 0], [1.0, 0, 0], [1.0, 1.0, 1.0]])
# B = np.array([0, 1, 2, 3])

# X = calcDihedral(B, A)
# print(45, X)

# A = np.array([[0, 1.0, 0], [0, 0, 0], [1.0, 0, 0], [1.0, -1.0, 1.0]])
# B = np.array([0, 1, 2, 3])

# X = calcDihedral(B, A)
# print(135, X)