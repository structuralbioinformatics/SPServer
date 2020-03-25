import math, numpy

def get_points_distance(point1, point2):
    '''
    Get the distance between two points
    '''
    point1 = numpy.array(point1)
    point2 = numpy.array(point2) 
    diff = point1-point2
    return numpy.sqrt(numpy.dot(diff, diff))

def translate_point(point, direction):
    '''
    '''
    x = point[0] + direction[0]
    y = point[1] + direction[1]
    z = point[2] + direction[2]
    return x, y, z

###########
#   HEX   #
###########

_AXES2TUPLE = { 'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
                'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
                'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
                'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
                'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
                'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
                'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
                'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1) }

_TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())

# axis sequences for Euler angles
_NEXT_AXIS = [1, 2, 0, 1]

# epsilon for testing whether a number is close to zero
_EPS = numpy.finfo(float).eps * 4.0

def get_hex_rotation_matrix(angles, axes='szyz'):
    '''
    ai, aj, ak : Euler's roll, pitch and yaw angles
    axes : One of 24 axis sequences as string or encoded tuple (default HEX: szyz)
    '''
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes.lower()]
    except (AttributeError, KeyError):
        _TUPLE2AXES[axes]  # validation
        firstaxis, parity, repetition, frame = axes

    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    ai = -angles[0]
    aj = -angles[1]
    ak = -angles[2]

    if frame:
        ai, ak = ak, ai
    if parity:
        ai, aj, ak = -ai, -aj, -ak

    si, sj, sk = math.sin(ai), math.sin(aj), math.sin(ak)
    ci, cj, ck = math.cos(ai), math.cos(aj), math.cos(ak)
    cc, cs = ci*ck, ci*sk
    sc, ss = si*ck, si*sk

    M = numpy.identity(3)
    if repetition:
        M[i, i] = cj
        M[i, j] = sj*si
        M[i, k] = sj*ci
        M[j, i] = sj*sk
        M[j, j] = -cj*ss+cc
        M[j, k] = -cj*cs-sc
        M[k, i] = -sj*ck
        M[k, j] = cj*sc+cs
        M[k, k] = cj*cc-ss
    else:
        M[i, i] = cj*ck
        M[i, j] = sj*sc-cs
        M[i, k] = sj*cc+ss
        M[j, i] = cj*sk
        M[j, j] = sj*ss+cc
        M[j, k] = sj*cs-sc
        M[k, i] = -sj
        M[k, j] = cj*si
        M[k, k] = cj*ci

    return M

def hex_transform_point(point, direction, angles):
    '''
    '''
    rotation_matrix = get_hex_rotation_matrix(angles)
    rotated_point = numpy.dot(point, rotation_matrix)
    return translate_point(rotated_point, direction)

############
# PACHDOCK #
############

def get_pachdock_rotation_matrix(angles):
    '''
    '''
    rotation_matrix = numpy.identity(3)

    cos_x = math.cos(angles[0])
    sin_x = math.sin(angles[0])
    cos_y = math.cos(angles[1])
    sin_y = math.sin(angles[1])
    cos_z = math.cos(angles[2])
    sin_z = math.sin(angles[2])

    rotation_matrix[0][0] = cos_z * cos_y
    rotation_matrix[1][0] = sin_z * cos_y
    rotation_matrix[2][0] = sin_y   
    rotation_matrix[0][1] = -sin_y * sin_x * cos_z - sin_z * cos_x
    rotation_matrix[1][1] = -sin_y * sin_x * sin_z + cos_x * cos_z
    rotation_matrix[2][1] = cos_y * sin_x
    rotation_matrix[0][2] = -sin_y * cos_x * cos_z + sin_z * sin_x
    rotation_matrix[1][2] = -sin_y * cos_x * sin_z - sin_x * cos_z
    rotation_matrix[2][2] = cos_y * cos_x

    return rotation_matrix

def pachdock_transform_point(point, direction, angles):
    '''
    transformed_point = rotation_matrix*point + direction
    '''
    rotation_matrix = get_pachdock_rotation_matrix(angles)
    rotated_point = numpy.dot(rotation_matrix, point)
    return translate_point(rotated_point, direction)

#def rotate_point_euclidean(point, angles):
#
#    '''
#    Rotate a 3D point using euclidean angles
#    '''
# 
#    rotation_matrix = numpy.zeros([3,3], float)
#
#    cos_x = math.cos(numpy.radians(angles[0]))
#    sin_x = math.sin(numpy.radians(angles[0]))
#    cos_y = math.cos(numpy.radians(angles[1]))
#    sin_y = math.sin(numpy.radians(angles[1]))
#    cos_z = math.cos(numpy.radians(angles[2]))
#    sin_z = math.sin(numpy.radians(angles[2]))
#
#    rotation_matrix[0][0] = cos_y * cos_z
#    rotation_matrix[1][0] = cos_y * sin_z
#    rotation_matrix[2][0] = -sin_y
#
#    rotation_matrix[0][1] = -cos_x * sin_z + sin_x * sin_y * cos_z
#    rotation_matrix[1][1] = cos_x * cos_z + sin_x * sin_y * sin_z
#    rotation_matrix[2][1] = sin_x * cos_y
#
#    rotation_matrix[0][2] = sin_x * sin_z + cos_x * sin_y * cos_z
#    rotation_matrix[1][2] = -sin_x * cos_z + cos_x * sin_y * sin_z
#    rotation_matrix[2][2] = cos_x * cos_y
#
#    x = rotation_matrix[0][0] * point[0] + rotation_matrix[0][1] * point[1] + rotation_matrix[0][2] * point[2]
#    y = rotation_matrix[1][0] * point[0] + rotation_matrix[1][1] * point[1] + rotation_matrix[1][2] * point[2]
#    z = rotation_matrix[2][0] * point[0] + rotation_matrix[2][1] * point[1] + rotation_matrix[2][2] * point[2]
#
#    return x, y, z
#
#def rotate_point_twists(point, angles):
#
#    '''
#    Rotate a 3D point using the Z, theta and phi twists (used by FTDock)
#    '''
#
#    #Z axis twist
#    z_twist = numpy.radians(angles[0])
#    z_twist_x = point[0] * math.cos(z_twist) - point[1] * math.sin(z_twist)
#    z_twist_y = point[0] * math.sin(z_twist) + point[1] * math.cos(z_twist)
#    z_twist_z = point[2]
#    #theta twist along plane of x-z
#    theta = numpy.radians(angles[1])
#    theta_x = z_twist_z * math.sin(theta) + z_twist_x * math.cos(theta)
#    theta_y = z_twist_y
#    theta_z = z_twist_z * math.cos(theta) - z_twist_x * math.sin(theta)
#    #phi twist around z axis
#    phi = numpy.radians(angles[2])
#    phi_x = theta_x * math.cos(phi) - theta_y * math.sin(phi)
#    phi_y = theta_x * math.sin(phi) + theta_y * math.cos(phi)
#    phi_z = theta_z
#
#    x = phi_x
#    y = phi_y
#    z = phi_z
#
#    return x, y, z
