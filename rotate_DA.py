##########################################
# function new_Pos = Rotate_DA(Pos, ang, SA, dt, iCA, mAID)
# 
# Rotates the coordinates by the dihedral angle given
#
# Input:
# Position: Coordinates of the dipeptide
# SetChi:  Angle to rotate to
# delta_term: pi*sign(InitChi)*InitChi/180
# iChiArray: indexes of the atoms that define the dihedral angle
# moveAtomID: Atoms to move during rotation
#
# Return:
# new_Pos: new coordinates of dipeptide
##########################################
def rotate_DA(Position, setChi, delta_term, iChiArray, moveAtomID):
    import numpy as np

    subtract_atom = Position[iChiArray[1], :]

    deltaChi1_F153 = delta_term - setChi * np.pi / 180.0
    TempPosition = Position - subtract_atom 			# Move all atoms so 2nd atom of dihedral is at orgin
    CAtoCB_F153 = -TempPosition[iChiArray[2], :]
    CAtoCB_F153 = CAtoCB_F153 / np.linalg.norm(CAtoCB_F153)
    q0 = np.cos(deltaChi1_F153 / 2)
    sindelta = np.sin(deltaChi1_F153 / 2)
    q1 = CAtoCB_F153[0] * sindelta
    q2 = CAtoCB_F153[1] * sindelta
    q3 = CAtoCB_F153[2] * sindelta
    q02 = q0 * q0
    q12 = q1 * q1
    q22 = q2 * q2
    q32 = q3 * q3

    Q = np.array([[(q02 + q12 - q22 - q32), 2.0 * (q1 * q2 - q0 * q3), 2 * (q0 * q2 + q1 * q3)],
                 [2 * (q1 * q2 + q0 * q3), (q02 - q12 + q22 - q32), 2 * (-q0 * q1 + q2 * q3)],
                 [2 * (-q0 * q2 + q1 * q3), 2 * (q0 * q1 + q2 * q3), (q02 - q12 - q22 + q32)]])

    move_coordinates = TempPosition[moveAtomID, :]
    newPos = np.matmul(Q, move_coordinates.T)               # Move the needed atoms
    TempPosition[moveAtomID, :] = newPos.T        # Put the moved atoms back into the original array
    new_Pos1 = TempPosition + subtract_atom     # Move back to original location
    return new_Pos1
