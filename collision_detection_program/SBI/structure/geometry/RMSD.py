"""
    Support Library to calculate RMSD difference between two Structure Objects

                                        jbonet @ boliva's lab        2012
"""
import math
import numpy as np


def superimposition_matrix_vector(structure1, structure2):

    Matrix = kabsch(structure1, structure2)

    # structure1 = structure1.duplicate(hetero = False, water = False, backbone = True)
    # structure1.rotate(Matrix)

    Vector = structure1.geometric_center(structure = True, hetero = False, water = False) - \
        structure2.geometric_center(structure = True, hetero = False, water = False)

    for i in range(3):
        for j in range(3):
            if (Matrix[i][j] > 0 and Matrix[i][j] < 1e-5) or (Matrix[i][j] < 0 and Matrix[i][j] > -1e-5):
                Matrix[i][j] = 0
    return Matrix, Vector


def kabsch(str1 = None, str2 = None):
    '''
    Given two structures, kabsch algorithm tries to identify the rotation matrix to move str1 and optimize the RMSD over str2
    '''

    '''
    STEP1: Center into origin
    '''
    centered_str1 = str1.duplicate(hetero = False, water = False, backbone = True)
    centered_str2 = str2.duplicate(hetero = False, water = False, backbone = True)

    centered_str1.translate_onto_origin()
    centered_str2.translate_onto_origin()

    '''
    STEP2: Retrieve all coordinates
    '''
    centered_str1_coordinates = centered_str1.get_backbone_coordinates()
    centered_str2_coordinates = centered_str2.get_backbone_coordinates()

    V, S, Wt = np.linalg.svd(np.dot(np.transpose(centered_str2_coordinates), centered_str1_coordinates))

    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    U = np.dot(V, Wt)

    return U


# def repositioned_RMSD(str1 = None, str2 = None, c_type = "all", selection = None):
#     '''
#     Checking the RMSD between two IDENTICAL structures implies calculating the difference between all their atoms.
#     Structures are assumed to have been previously positioned as they are expected to be for the comparison

#     This definition accepts that some flexibility might have been applied to the structures
#     '''

#     '''
#     For lack of a better/faster way to check it out, we will consider that two structures are the same if they
#     have the same number of atoms.
#     Else, we will rise a ValueError
#     '''
#     if str1.get_structure_length() != str2.get_structure_length():
#         raise ValueError("Both structures need to be the same for the RMSD of the reposition to be calculated\n")

#     '''
#     Several flags can be called over this function:
#         c_type: the allowed values are:
#             all: it will compare the full RMSD
#             ca: it will compare the RMSD between the CA atoms
#             cb: it will compare the RMSD between the CB atoms -> GLY are ignored
#             cb-ca: CB RMSD - CA RMSD -> kind of compares only how the lateral chains have changed
#             any other option defaults to 'all'
#         selection: an array of positions can be submitted and only those Aa are going to be compared
#                    This can be used, for example, if comparing the Interface RMSD
#     '''
#     allowed_c_type = set(['all', 'ca', 'cb', 'cb-ca'])
#     if c_type not in allowed_c_type:
#         c_type = 'all'

#     if selection is not None:
#         selection_set = set(selection)

#     E = 0
#     atom_count = 0
#     for res in xrange(str1.get_structure_length()):
#         if selection is None or str1.get_residues()[res].get_residue_num() in selection_set:
#             if c_type == 'all':
#                 for at in xrange(str1.get_residues()[res].get_residue_size()):
#                     distance = str1.get_residues()[res].get_atoms()[at].get_distance(atom2=str2.get_residues()[res].get_atoms()[at])
#                     atom_count += 1
#                     E += math.pow(distance, 2)
#             if c_type == 'ca':
#                 distance = str1.get_residues()[res].get_ca().get_distance(atom2=str2.get_residues()[res].get_ca())
#                 atom_count += 1
#                 E += math.pow(distance, 2)
#             if c_type == 'cb' and str1.get_residues()[res].get_cb() is not None:
#                 distance = str1.get_residues()[res].get_cb().get_distance(atom2=str2.get_residues()[res].get_cb())
#                 atom_count += 1
#                 E += math.pow(distance, 2)
#             if c_type == 'cb-ca' and str1.get_residues()[res].get_cb() is not None:
#                 distance_a = str1.get_residues()[res].get_ca().get_distance(atom2=str2.get_residues()[res].get_ca())
#                 distance_b = str1.get_residues()[res].get_cb().get_distance(atom2=str2.get_residues()[res].get_cb())
#                 atom_count += 12
#                 E += math.pow(distance_b - distance_a, 2)

#     E = E / atom_count

#     return math.sqrt(math.fabs(E))


# def repositioned_RMSD_byMatrix(AV = None, AV2 = None, NormL = None, rotation_matrix1 = None, rotation_matrix2 = None,
#                                translation_vector1 = None, translation_vector2 = None):
#     '''
#     This method will try to implement an RMSD calculation based on the rotation matrices and the translation vectors
#     '''

#     x, y, z, V2 = calculateVandV2(translation_vector1 = translation_vector1, translation_vector2 = translation_vector2)
#     V  = [x, y, z]

#     V1 = calculateV1(rotation_matrix1 = rotation_matrix1, rotation_matrix2 = rotation_matrix2, AV2 = AV2)
#     V3 = calculateV3(V = V, rotation_matrix1 = rotation_matrix1, rotation_matrix2 = rotation_matrix2, AV = AV)

#     '''
#     Calculate RMSD
#     '''
#     rmsd = V2 + NormL - V1 + V3
#     return math.sqrt(math.fabs(rmsd))

# """
#     __Suporting Function to repositioned_RMSD_byMatrix__
# """


# def calculateVandV2(translation_vector1 = None, translation_vector2 = None):
#     '''
#     The translation vector clusters the data in order x, y, z
#     This one creates the vectors:
#         V:    x, y, z
#     And the value V2
#     '''
#     x = translation_vector1[0] - translation_vector2[0]
#     y = translation_vector1[1] - translation_vector2[1]
#     z = translation_vector1[2] - translation_vector2[2]

#     v2 = x * x + y * y + z * z

#     return x, y, z, v2


# def calculateV3(V = None, rotation_matrix1 = None, rotation_matrix2 = None, AV = None):
#     '''
#     This one calculates V3 from the arrays V and AV and the rotation matrices
#         AV[0] : average_x
#         AV[1] : average_y
#         AV[2] : average_z
#         V[0]  : x
#         V[1]  : y
#         V[2]  : z
#     '''

#     v3 = 0

#     v3 = v3 + ((rotation_matrix1[0][0] - rotation_matrix2[0][0]) * AV[0] +
#                (rotation_matrix1[0][1] - rotation_matrix2[0][1]) * AV[1] +
#                (rotation_matrix1[0][2] - rotation_matrix2[0][2]) * AV[2]) * V[0]
#     v3 = v3 + ((rotation_matrix1[1][0] - rotation_matrix2[1][0]) * AV[0] +
#                (rotation_matrix1[1][1] - rotation_matrix2[1][1]) * AV[1] +
#                (rotation_matrix1[1][2] - rotation_matrix2[1][2]) * AV[2]) * V[1]
#     v3 = v3 + ((rotation_matrix1[2][0] - rotation_matrix2[2][0]) * AV[0] +
#                (rotation_matrix1[2][1] - rotation_matrix2[2][1]) * AV[1] +
#                (rotation_matrix1[2][2] - rotation_matrix2[2][2]) * AV[2]) * V[2]

#     return 2 * v3


# def calculateV1(rotation_matrix1 = None, rotation_matrix2 = None, AV2 = None):
#     '''
#     This one calculates V1 from the rotation matrices and the AV2 vector
#     '''

#     v1 = 0

#     v1 = v1 + (rotation_matrix1[0][0] * rotation_matrix2[0][0] +
#                rotation_matrix1[1][0] * rotation_matrix2[1][0] +
#                rotation_matrix1[2][0] * rotation_matrix2[2][0]) * AV2[0]  # xx
#     v1 = v1 + (rotation_matrix1[0][1] * rotation_matrix2[0][1] +
#                rotation_matrix1[1][1] * rotation_matrix2[1][1] +
#                rotation_matrix1[2][1] * rotation_matrix2[2][1]) * AV2[1]  # yy
#     v1 = v1 + (rotation_matrix1[0][2] * rotation_matrix2[0][2] +
#                rotation_matrix1[1][2] * rotation_matrix2[1][2] +
#                rotation_matrix1[2][2] * rotation_matrix2[2][2]) * AV2[2]  # zz

#     v1 = v1 + (rotation_matrix1[0][0] * rotation_matrix2[0][1] +
#                rotation_matrix1[0][1] * rotation_matrix2[0][0] +
#                rotation_matrix1[1][0] * rotation_matrix2[1][1] +
#                rotation_matrix1[1][1] * rotation_matrix2[1][0] +
#                rotation_matrix1[2][0] * rotation_matrix2[2][1] +
#                rotation_matrix1[2][1] * rotation_matrix2[2][0]) * AV2[3]  # xy
#     v1 = v1 + (rotation_matrix1[0][0] * rotation_matrix2[0][2] +
#                rotation_matrix1[0][2] * rotation_matrix2[0][0] +
#                rotation_matrix1[1][0] * rotation_matrix2[1][2] +
#                rotation_matrix1[1][2] * rotation_matrix2[1][0] +
#                rotation_matrix1[2][0] * rotation_matrix2[2][2] +
#                rotation_matrix1[2][2] * rotation_matrix2[2][0]) * AV2[4]  # xz
#     v1 = v1 + (rotation_matrix1[0][1] * rotation_matrix2[0][2] +
#                rotation_matrix1[0][2] * rotation_matrix2[0][1] +
#                rotation_matrix1[1][1] * rotation_matrix2[1][2] +
#                rotation_matrix1[1][2] * rotation_matrix2[1][1] +
#                rotation_matrix1[2][1] * rotation_matrix2[2][2] +
#                rotation_matrix1[2][2] * rotation_matrix2[2][1]) * AV2[5]  # yz

#     return 2 * v1


# def calculateAVandAV2andNormL(str1):
#     '''
#     This one creates the vectors:
#         AV:     average_x, average_y, average_z
#         AV2:    xx, yy, zz, xy, xz, yz
#     and the NormL value
#     '''

#     average_x, average_y, average_z, total_atoms = 0, 0, 0, 0
#     xx, yy, zz, xy, xz, yz = 0, 0, 0, 0, 0, 0

#     for residue in str1.get_residues():
#         for atom in residue.get_atoms():
#             average_x = average_x + atom.get_x()
#             average_y = average_y + atom.get_y()
#             average_z = average_z + atom.get_z()

#             xx = xx + atom.get_x() * atom.get_x()
#             yy = yy + atom.get_y() * atom.get_y()
#             zz = zz + atom.get_z() * atom.get_z()
#             xy = xy + atom.get_x() * atom.get_y()
#             xz = xz + atom.get_x() * atom.get_z()
#             yz = yz + atom.get_y() * atom.get_z()

#             total_atoms += 1

#     average_x = average_x / float(total_atoms)
#     average_y = average_y / float(total_atoms)
#     average_z = average_z / float(total_atoms)

#     xx = xx / float(total_atoms)
#     yy = yy / float(total_atoms)
#     zz = zz / float(total_atoms)
#     xy = xy / float(total_atoms)
#     xz = xz / float(total_atoms)
#     yz = yz / float(total_atoms)

#     NormL = 2 * (xx + yy + zz)

#     return average_x, average_y, average_z, xx, yy, zz, xy, xz, yz, NormL
