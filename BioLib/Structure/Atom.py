from BioLib.Algebra.Transforms import get_points_distance, hex_transform_point, pachdock_transform_point
from BioLib.Tools.BioExceptions import RelocateProgramError

class Atom(object):
    '''
    Represents an structure atom
    '''
    def __init__(self, num, atom_type, x, y, z):
        '''
        Constructor
        '''
        self.num = num
        self.type = atom_type
        self.x = x
        self.y = y
        self.z = z

    def get_num(self):
        return self.num

    def get_type(self):
        return self.type

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    def get_coords(self):
        return (self.x, self.y, self.z)

    def get_distance(self, atom):
        return get_points_distance(self.get_coords(), atom.get_coords())
 
    def get_distance_coords(self, coords):
        return get_points_distance(self.get_coords(), coords)

    def relocate(self, translation_vector, rotation_vector, docking='HEX'):
        '''
        Relocate the atom position (x, y, z). 
        @translation_vector = Vector used to translate the atom
        @Rotation_vector = Vector used to rotate the structure
        @docking = Docking program used to get the vectors
        '''
        if docking == 'HEX':
            self.x, self.y, self.z = hex_transform_point(self.get_coords(), translation_vector, rotation_vector)
        elif docking == 'PATCHDOCK':
            self.x, self.y, self.z = pachdock_transform_point(self.get_coords(), translation_vector, rotation_vector)
        else:
            raise RelocateProgramError(docking)

    def __str__(self):
        '''
        Return the atom information
        '''
        return 'ATOM%s  %s%s%s%s' % (str(self.get_num()).rjust(7), 
                                     self.get_type().ljust(3), 
                                     ('%.3f' % float(self.get_x())).rjust(12), 
                                     ('%.3f' % float(self.get_y())).rjust(8), 
                                     ('%.3f' % float(self.get_z())).rjust(8))
