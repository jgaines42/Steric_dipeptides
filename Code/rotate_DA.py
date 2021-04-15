##########################################
# function new_Pos1 = Rotate_DA(Position, setChi, delta_term, iChiArray, moveAtomID)
#
# Rotates the coordinates by the dihedral angle given
#
# Input:
#   Position: Coordinates of the dipeptide
#   SetChi:  Angle to rotate to
#   delta_term: pi*sign(InitChi)*InitChi/180
#   iChiArray: indexes of the atoms that define the dihedral angle
#   moveAtomID: Atoms to move during rotation
#
# Returns:
#   new_Pos1: new coordinates of dipeptide
##########################################
def rotate_DA(Position, setChi, delta_term, iChiArray, moveAtomID):

    import numpy as np
    from math import sin, cos
    # Calculate how much the dihedral needs to change (in radians)
    deltaChi1_F153 = delta_term - setChi * np.pi / 180.0

    # Get the coordinates of the 2nd atom in the dihedral
    subtract_atom = Position[iChiArray[1], :]

    # Move all atoms so 2nd atom of dihedral is at orgin
    TempPosition = Position - subtract_atom

    # Get the vector from the origin to the 3rd atom of the dihedral.
    # This is the vector that we will rotate around
    CAtoCB_F153 = -TempPosition[iChiArray[2], :]
    CAtoCB_F153 = CAtoCB_F153 / np.linalg.norm(CAtoCB_F153)

    # Do complicated math
    q0 = cos(deltaChi1_F153 / 2.0)
    sindelta = sin(deltaChi1_F153 / 2.0)
    q1 = CAtoCB_F153[0] * sindelta
    q2 = CAtoCB_F153[1] * sindelta
    q3 = CAtoCB_F153[2] * sindelta
    q02 = q0 * q0
    q12 = q1 * q1
    q22 = q2 * q2
    q32 = q3 * q3

    # Q is the rotation vector
    Q = np.array([[(q02 + q12 - q22 - q32), 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q0 * q2 + q1 * q3)],
                 [2.0 * (q1 * q2 + q0 * q3), (q02 - q12 + q22 - q32), 2.0 * (-q0 * q1 + q2 * q3)],
                 [2.0 * (-q0 * q2 + q1 * q3), 2 * (q0 * q1 + q2 * q3), (q02 - q12 - q22 + q32)]])

    # Extract the coordinates that should move
    move_coordinates = TempPosition[moveAtomID, :]

    # Use Q to rotate these atoms
    newPos = np.matmul(Q, move_coordinates.T)

    # Put the moved atoms back into the original array
    TempPosition[moveAtomID, :] = newPos.T

    # Translate atoms back to original location
    new_Pos1 = TempPosition + subtract_atom
    return new_Pos1
