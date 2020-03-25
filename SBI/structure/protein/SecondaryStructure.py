"""
jbonet @ boliva's lab 2013
"""
import  numpy               as np
from    numpy import linalg as LA
import  copy, sys

class SecondaryStructure(object):

    min_ss_length     = {'H':4, 'E':2, 'G':3}
    max_ini_distance  = {'H':8, 'E':2, 'G':6}
    structure_regions = set(['H','G','E'])

    def __init__(self, sstype, initposition):

        self._sstype                                        = sstype
        self._inip                                          = initposition
        self._struct                                        = None
        self._endp                                          = None
        self._length                                        = None
        self._f11, self._cmf11, self._eigf11, self._tsf11   = None, None, None, None
        self._f44, self._cmf44, self._eigf44, self._tsf44   = None, None, None, None

    @property
    def structure(self):
        return self._struct
        
    def get_moment_of_inertia_length(self, workp):
        if   workp == 'f11':    (p1, p2) = (self._get_coordinate(self._f11),  self._get_coordinate(self._endp))
        elif workp == 'f44':    (p1, p2) = (self._get_coordinate(self._inip), self._get_coordinate(self._f44) )
        return p2 - p1 + 1

    def _get_coordinate(self, identifier):
        return self._struct._get_structure_array_coordinate(identifier)

    def calculate_center_of_masses(self):
        _struct = self._struct.duplicate(backbone = True)
        if self.max_ini_distance[self._sstype] >= self._length:
            self._f11   = self._inip
            self._f44   = self._endp
            self._cmf11 = _struct.geometric_center()
            self._cmf44 = self._cmf11
        else:
            self._f11   = _struct.aminoacids[ len(_struct.aminoacids) - self.max_ini_distance[self._sstype] ].identifier
            self._f44   = _struct.aminoacids[ self.max_ini_distance[self._sstype] - 1 ].identifier
            self._cmf11 = _struct.extract(self._f11, self._endp).geometric_center()
            self._cmf44 = _struct.extract(self._inip, self._f44).geometric_center()

        self._jacobi_angles('f11')
        if self._f11 != self._inip and self._f44 != self._endp:
            self._jacobi_angles('f44')
            if self._sstype == 'E': self._process_betas()
        else:
            self._eigf44 = self._eigf11 

        self._orientvectors()

    def _orientvectors(self):
        first_halfdif   = np.subtract(self._struct.first_aminoacid.ca.coordinates, self._cmf44)
        second_halfdif  = np.subtract(self._struct.last_aminoacid.ca.coordinates,  self._cmf11)
    
        first_sign      = -1 * np.sign(np.dot(first_halfdif,  self._eigf44))
        second_sign     = np.sign(np.dot(second_halfdif, self._eigf11))

        self._eigf44    = np.multiply(self._eigf44, first_sign)
        self._eigf11    = np.multiply(self._eigf11, second_sign)

        first_lambda    = np.dot(first_halfdif,  self._eigf44)
        second_lambda   = np.dot(second_halfdif, self._eigf11)

        self._tsf44     = np.add(np.multiply(first_lambda,  self._eigf44), self._cmf44)
        self._tsf11     = np.add(np.multiply(second_lambda, self._eigf11), self._cmf11)

    def _jacobi_angles(self, workp):
        _struct = self._struct.duplicate(backbone = True)
        if workp == 'f11':
            moving_point = self._f11
            fixed_point  = self._endp
            ini_coord    = _struct.extract(self._f11, self._endp)._all_atoms_coordinates()
            distance     = len(ini_coord)/3
            cm           = self._cmf11
        elif workp == 'f44':
            moving_point = self._f44
            ini_coord    = _struct.extract(self._inip, self._f44)._all_atoms_coordinates()
            distance     = len(ini_coord)/3
            cm           = self._cmf44

        new_coord = np.subtract(ini_coord, cm)
        x2, y2, z2, xy, yz, zx = 0, 0, 0, 0, 0, 0
        if self._sstype != 'E':
            for row in new_coord:
                x, y, z = row
                x2 += np.power(x, 2)
                y2 += np.power(y, 2)
                z2 += np.power(z, 2)
                xy += np.multiply(x, y)
                yz += np.multiply(y, z)
                zx += np.multiply(z, x)
        else:
            for i in range(0, len(ini_coord) - 4, 3):
                j=i+3
                x0, y0, z0 = new_coord[i:i+3][:,0], new_coord[i:i+3][:,1], new_coord[i:i+3][:,2]
                x1, y1, z1 = new_coord[j:j+3][:,0], new_coord[j:j+3][:,1], new_coord[j:j+3][:,2]

                xm, ym, zm = np.mean([x0, x1], axis = 0), np.mean([y0, y1], axis = 0), np.mean([z0, z1], axis = 0)

                x2 += np.sum(np.power(xm, 2))
                y2 += np.sum(np.power(ym, 2))
                z2 += np.sum(np.power(zm, 2))
                xy += np.sum(np.multiply(xm, ym))
                yz += np.sum(np.multiply(ym, zm))
                zx += np.sum(np.multiply(zm, xm))

        a = np.matrix([[y2 + z2,     -xy,     -zx],
                       [    -xy, x2 + z2,     -yz],
                       [    -zx,     -yz, x2 + y2]])

        eigenVal, eigenVec = LA.eig(a)
        #sort eigenVal / eigenVec by descending eigenVal value
        #get the smallest value as represents the axis that actually follows the structure direction
        idx = eigenVal.argsort()[::-1]
        eigenVal = eigenVal[idx]
        eigenVec = eigenVec[:,idx]

        if   workp == 'f11':    self._eigf11 = np.asarray(eigenVec[:,2]).reshape(-1)
        elif workp == 'f44':    self._eigf44 = np.asarray(eigenVec[:,2]).reshape(-1)

    def _process_betas(self):
        advance_limit        = 2 #plus the two already selected makes a total of 4
        degree_dif_threshold = 10

        end, advance        = False, 1
        original_eigienf11  = self._eigf11
        while not end:
            new_me          = copy.deepcopy(self)
            _struct         = new_me._struct.duplicate(backbone = True)
            new_me._f11     = new_me._struct.aminoacids[ len(new_me._struct.aminoacids) - (self.max_ini_distance[new_me._sstype]  + advance)].identifier
            new_me._cmf11   = _struct.extract(new_me._f11, new_me._endp).geometric_center()

            new_me._jacobi_angles('f11')
            difference = np.degrees(np.arccos(np.absolute(np.dot(new_me._eigf11, original_eigienf11))))
            if new_me._f11 == new_me._inip or advance >= advance_limit: end = True
            advance += 1
            if difference < degree_dif_threshold:
                self._f11    = new_me._f11
                self._cmf11  = new_me._cmf11
                self._eigf11 = new_me._eigf11

        end, advance        = False, 1
        original_eigienf44  = self._eigf44
        while not end:
            new_me          = copy.deepcopy(self)
            _struct         = new_me._struct.duplicate(backbone = True)
            new_me._f44     = new_me._struct.aminoacids[ self.max_ini_distance[new_me._sstype] - 1 + advance].identifier
            new_me._cmf44   = _struct.extract(new_me._inip, new_me._f44).geometric_center()

            new_me._jacobi_angles('f44')
            difference = np.degrees(np.arccos(np.absolute(np.dot(original_eigienf44, new_me._eigf44))))
            if new_me._f44 == new_me._endp or advance >= advance_limit: end = True
            advance += 1
            if difference < degree_dif_threshold:
                self._f44    = new_me._f44
                self._cmf44  = new_me._cmf44
                self._eigf44 = new_me._eigf44

    def headerformat(self, ssnum = 1, ssidentifier = 'A'):
        inires = self._struct.aminoacids[0]
        endres = self._struct.aminoacids[-1]
        if self._sstype == 'H' or self._sstype == 'G':
            ini = '{1.type} {0._struct.chain}{1.number:>5d}{1.version}'.format(self, inires)
            end = '{1.type} {0._struct.chain}{1.number:>5d}{1.version}'.format(self, endres)
            return 'HELIX {0:>4d} {1:>3s} {2} {3} 1 {4:>35d}'.format(ssnum, ssidentifier, ini, end, self._length)
        elif self._sstype == 'E':
            ini = '{1.type} {0._struct.chain}{1.number:>4d}{1.version}'.format(self, inires)
            end = '{1.type} {0._struct.chain}{1.number:>4d}{1.version}'.format(self, endres)
            return 'SHEET {0:>4d} {1:>3s} 1 {2} {3} 0'.format(ssnum, ssidentifier, ini, end)
        else:
            raise NotImplementedError('Unknown secondary structure type!')

    def strdata(self, workp):
        data = []
        data.append("( {0._sstype} ) {0._inip:>4s} <-- {0._length:>2d} --> {0._endp:>4s}".format(self))
        if workp == 'f11':
            data.append("\tf11: {0._f11:>4s} cmf11: {0._cmf11:>50s} eigf11: {0._eigf11:>50s}".format(self))
        elif workp == 'f44':
            data.append("\tf44: {0._f44:>4s} cmf44: {0._cmf44:>50s} eigf44: {0._eigf44:>50s}".format(self))
        return "\n".join(data)

    def __repr__(self):
        data = []
        data.append("( {0._sstype} ) {0._inip:>4s} <-- {0._length:>2d} --> {0._endp:>4s}".format(self))
        data.append("\tf11: {0._f11:>4s} cmf11: {0._cmf11:>50s} eigf11: {0._eigf11:>50s}".format(self))
        data.append("\tf44: {0._f44:>4s} cmf44: {0._cmf44:>50s} eigf44: {0._eigf44:>50s}".format(self))
        # data.append(repr(self._struct))
        return "\n".join(data)
