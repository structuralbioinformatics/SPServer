"""
jbonet @ boliva's lab 2013
"""
from SBI.beans import StorableObject
from SBI.beans import JSONer
from SBI.beans import File

import numpy as np
from collections import Counter


class Arch(StorableObject, JSONer):

    def __init__(self, source, ss1, ss2, internalss, distanceAA, structure,
                 torsionsCA, PsiPhi, order):
        self._source = source
        self._ss1    = ss1  # N-term
        self._ss2    = ss2  # C-term
        self._intss  = internalss  # number of internal secondary structures
        self._inttxt = []          # header format of internal ss (if any)
        self._inttyp = ''          # internal structures type (if any)
        self._distAA = distanceAA
        self._str    = structure
        self._vlo    = self._ss2._tsf44 - self._ss1._tsf11
        self._dist   = np.sqrt(np.dot(self._vlo, self._vlo))
        self._angls  = [None, None, None]  # theta, rho, delta
        self._torCA  = torsionsCA
        self._psiphi = PsiPhi
        self._order  = order
        self._tordis = self._distance_torsion_vector()
        self._angles()
        # theta, rho, delta
        self._torgls = [self._angle_torsion_vector(self._angls[0], 180),
                        self._angle_torsion_vector(self._angls[1], 360),
                        self._angle_torsion_vector(self._angls[2], 180)]
        self._htotal = 0
        self._hss    = 0
        self._hloop  = 0
        self._hdim   = 0
        self._hmxini = 0
        self._hmatix = self._calculate_hydrogen_bond_matrix()

        distribution = Counter(list(self.structure_sequence))
        self._gaps   = distribution['x'] if 'x' in distribution else 0

    @property
    def type(self):
        if self._ss1._sstype == 'E' and self._ss2._sstype == 'E':
            if self._hss < 1:
                # If this condition was == 0 some BK of length 1 can become BN
                return 'BK'
            else:
                if self.aminoacid_distance == 1 and self.theta < 160:
                    return 'BK'
                if self.theta < 90:
                    return 'BK'
                return 'BN'
        else:
            return self._ss1._sstype + self._ss2._sstype

    @property
    def boundaries(self):
        return self.initial_position + '-' + self.end_position

    @property
    def initial_position(self):
        return str(self._ss1._inip).strip()

    @property
    def end_position(self):
        return str(self._ss2._endp).strip()

    @property
    def selfboundaries(self):
        return '0-' + str(len(self._torCA) - 1)

    @property
    def access_surface(self):
        sec = []
        for residues in self._str.aminoacids:
            sec.append(residues.accessibilitycoded)
        return ''.join(sec)

    @property
    def psiphi(self):
        data = ['', '', '']
        for i in range(0, len(self._psiphi)-4, 3):
            data[0] += self._psiphi[i]
            data[1] += self._psiphi[i+1]
            data[2] += self._psiphi[i+2]
        return data

    @property
    def cartesian_distance(self):
        return self._dist

    @property
    def aminoacid_distance(self):
        return self._distAA

    @property
    def length(self):
        return self._distAA

    @property
    def internal_structures(self):
        return self._intss

    @property
    def gaps(self):
        return self._gaps

    @property
    def identifier(self):
        ident = [self._source, str(self._ss1._inip).strip()]
        if self.is_superarch:
            ident.append(str(self._intss))
        return "_".join(ident)

    @property
    def theta(self):
        return self._angls[0]

    @property
    def rho(self):
        return self._angls[1]

    @property
    def delta(self):
        return self._angls[2]

    @property
    def aminoacid_sequence(self):
        return self._str.gapped_protein_sequence

    @property
    def structure_sequence(self):
        return self._str.gapped_protein_secondary_structure

    @property
    def structure(self):
        return self._str

    @property
    def length_ss1(self):
        return self._ss1._length

    @property
    def length_ss2(self):
        return self._ss2._length

    @property
    def has_gaps(self):
        return self._gaps > 0

    @property
    def is_superarch(self):
        return self._intss > 0

    def unbound_ss(self, percentage):
        if not self.has_gaps:
            return False
        return (self._gaps/self._distAA)*100 >= percentage

    def _calculate_hydrogen_bond_matrix(self):
        matrix_min_border_index = self._str._get_structure_array_coordinate(self._ss1._endp)
        matrix_max_border_index = self._str._get_structure_array_coordinate(self._ss2._inip)
        ss1_end_index           = self._str._get_structure_array_coordinate(self._ss1._endp)
        ss2_ini_index           = self._str._get_structure_array_coordinate(self._ss2._inip)

        if matrix_min_border_index > 1:
            matrix_min_border_index -= 2
        else:
            matrix_min_border_index  = 0
        if matrix_max_border_index < len(self._str.aminoacids) - 2:
            matrix_max_border_index += 2
        else:
            matrix_max_border_index = len(self._str.aminoacids) - 1

        matrix_size = matrix_max_border_index - matrix_min_border_index + 1
        matrix = np.zeros(shape=(matrix_size, matrix_size))
        ii = 0
        for i in range(matrix_min_border_index, matrix_min_border_index + matrix_size):
            jj = 0
            for j in range(matrix_min_border_index, matrix_min_border_index + matrix_size):
                diff = j-i
                if self._str.aminoacids[i].dssp.nhoa[0] != diff and self._str.aminoacids[i].dssp.nhob[0] != diff:
                    matrix[ii][jj] = 0
                elif self._str.aminoacids[i].dssp.nhoa[0] == diff and self._str.aminoacids[i].dssp.nhob[0] == diff:
                    matrix[ii][jj]  = self._str.aminoacids[i].dssp.nhoa[1]
                    matrix[ii][jj] += self._str.aminoacids[i].dssp.nhob[1]
                elif self._str.aminoacids[i].dssp.nhoa[0] == diff and self._str.aminoacids[i].dssp.nhob[0] != diff:
                    matrix[ii][jj]  = self._str.aminoacids[i].dssp.nhoa[1]
                else:
                    matrix[ii][jj]  = self._str.aminoacids[i].dssp.nhob[1]
                jj += 1
            ii += 1
        jj = 0
        for j in range(matrix_min_border_index, matrix_min_border_index + matrix_size):
            ii = 0
            for i in range(matrix_min_border_index, matrix_min_border_index + matrix_size):
                diff = i-j
                if self._str.aminoacids[j].dssp.ohna[0] != diff and self._str.aminoacids[j].dssp.ohnb[0] != diff:
                    matrix[ii][jj] += 0
                elif self._str.aminoacids[j].dssp.ohna[0] == diff and self._str.aminoacids[j].dssp.ohnb[0] == diff:
                    matrix[ii][jj] += self._str.aminoacids[j].dssp.ohna[1]
                    matrix[ii][jj] += self._str.aminoacids[j].dssp.ohnb[1]
                    matrix[ii][jj]  = matrix[ii][jj]/2
                elif self._str.aminoacids[j].dssp.ohna[0] == diff and self._str.aminoacids[j].dssp.ohnb[0] != diff:
                    matrix[ii][jj] += self._str.aminoacids[j].dssp.ohna[1]
                    matrix[ii][jj]  = matrix[ii][jj]/2
                else:
                    matrix[ii][jj] += self._str.aminoacids[j].dssp.ohnb[1]
                    matrix[ii][jj]  = matrix[ii][jj]/2
                ii += 1
            jj += 1
        hbond_total, hbond_ss, hbond_loop = 0, 0, 0
        for j in range(matrix_size):
            for i in range(matrix_size):
                if matrix[i][j] <= -0.01:
                    hbond_total += 1
                if matrix[i][j] <= -0.6:
                    position_i = matrix_min_border_index + i
                    position_j = matrix_min_border_index + j
                    if (position_i <= ss1_end_index and position_j >= ss2_ini_index) or \
                       (position_i >= ss2_ini_index and position_j <= ss1_end_index):
                        hbond_ss += 1
                    elif (position_i > ss1_end_index and position_i < ss2_ini_index) or \
                         (position_j > ss1_end_index and position_j < ss2_ini_index):
                            hbond_loop += 1

        self._htotal = hbond_total
        self._hss    = hbond_ss
        self._hloop  = hbond_loop
        self._hdim   = matrix_size
        self._hmxini = matrix_min_border_index

        return matrix

    def _distance_torsion_vector(self):

        data    = ["-" for i in range(len(self._torCA))]
        border1 = self._ss1._length - 1
        border2 = border1 + self._distAA + 1
        if self._dist   <= 7.5:
            data[border1], data[border2] = 's', 's'
        elif self._dist <= 10:
            data[border1], data[border2] = 'm', 'm'
        elif self._dist <= 20:
            data[border1], data[border2] = 'l', 'l'
        else:
            data[border1], data[border2] = 'L', 'L'
        return ''.join(data)

    def _angle_torsion_vector(self, angle, max_range):
        data     = ["-" for i in range(len(self._torCA))]
        border1  = self._ss1._length - 1
        border2  = border1 + self._distAA + 1
        accuracy = 36
        definitions = 'ABCDEFGHIJKLMNOPQRSTUVXYZabcdefghijk'
        for i in range(accuracy):
            if angle >= i * (max_range/accuracy) and angle < (i+1) * (max_range/accuracy):
                break
        data[border1] = definitions[i]
        data[border2] = definitions[i]
        return ''.join(data)

    def _angles(self):

        normal = self._vectorproduct(self._vlo, self._ss1._eigf11)
        target = self._vectorproduct(self._ss1._eigf11, normal)
        const  = np.dot(target, target)

        if const != 0:
            self._angls[0] = np.degrees(np.arccos(np.dot(self._ss1._eigf11, self._ss2._eigf44)))

            const2 = np.dot(self._ss2._eigf44, self._ss1._eigf11)
            proj   = np.subtract(self._ss2._eigf44, np.multiply(const2, self._ss1._eigf11))
            dproj  = np.sqrt(np.dot(proj, proj))

            if dproj != 0:
                angle  = np.degrees(np.arccos(np.dot(proj, normal)/dproj))
                const3 = np.dot(proj, target)
                if const3 < 0:
                    angle = 360 - angle
                self._angls[1] = angle
            else:
                self._angls[1] = 0

        else:
            angle = np.degrees(np.arccos(np.dot(self._ss1._eigf11, self._ss2._eigf44)))
            self._angls[0] = angle
            self._angls[1] = 500

        self._angls[2] = np.degrees(np.arccos(np.dot(self._ss1._eigf11,self._vlo)/self._dist))

    def _vectorproduct(self, v, n):

        retvec    = [0, 0, 0]
        retvec[0] = v[1]*n[2] - v[2]*n[1]
        retvec[1] = v[2]*n[0] - v[0]*n[2]
        retvec[2] = v[0]*n[1] - v[1]*n[0]

        return np.divide(retvec, np.sqrt(np.dot(retvec,retvec)))

    def _format_hbond_matrix(self):
        data = []
        data.append('______|')
        data.append('      |')
        ii = 0
        for i in range(self._hmxini, self._hmxini + self._hdim):
            data[0] = data[0] + '{0:>5d}|'.format(i)
            data[1] = data[1] + '-----+'
            data.append('{0:>4d}  |'.format(i))
            data.append('      |')
            jj = 0
            for j in range(self._hmxini, self._hmxini + self._hdim):
                data[-2] = data[-2] + '{0:>5.1f}|'.format(self._hmatix[ii][jj])
                data[-1] = data[-1] + '-----+'
                jj += 1
            ii += 1
        return "\n".join(data)

    def _format_border_hbond_matrix(self):
        data = []
        data.append('______|')
        data.append('      |')
        ii = 0
        for i in range(self._hmxini, self._hmxini + self._hdim):
            if i < self._ss1._length or i >= self._ss1._length + self._distAA:
                data[0] = data[0] + '{0:>5d}|'.format(i)
                data[1] = data[1] + '-----+'
                data.append('{0:>4d}  |'.format(i))
                data.append('      |')
                jj = 0
                for j in range(self._hmxini, self._hmxini + self._hdim):
                    if j < self._ss1._length or j >= self._ss1._length + self._distAA:
                        data[-2] = data[-2] + '{0:>5.1f}|'.format(self._hmatix[ii][jj])
                        data[-1] = data[-1] + '-----+'
                    jj += 1
            ii += 1
        return "\n".join(data)

    def format2file(self, filename, extension = 'pdb', center = False):
        if extension not in ('pdb', 'js'):
            raise AttributeError('Not accepted extension')
        structure = File('.'.join([filename, extension]), 'w')
        if extension == 'pdb':
            structure.write(self.pdb_format(center = center))
        elif extension == 'js':
            structure.write(self.js_format(center = center))
        structure.close()

    def pdb_format(self, center = False):
        data = []
        data.append(self._ss1.headerformat())
        for internalss in self._inttxt:
            data.append(internalss)
        data.append(self._ss2.headerformat(self._intss + 2, chr(65 + self._intss + 1)))
        tempstr = self._str.duplicate(hetero = False, water = False)
        if center:
            tempstr.translate_onto_origin()
        data.append(tempstr.PDB_format())
        return "\n".join(data)

    def js_format(self, center = False):
        init = 'var {0}="'.format('pdb_' + self.identifier)
        return init + self.pdb_format(center = center).replace('\n', '\\n') + '"'

    def as_dict(self):
        return {'id'         : self.identifier,
                'type'       : self.type,
                'aadistance' : self.aminoacid_distance,
                'distance'   : self.cartesian_distance,
                'inipos'     : self.initial_position,
                'endpos'     : self.end_position,
                'theta'      : self.theta,
                'rho'        : self.rho,
                'delta'      : self.delta,
                'gaps'       : self.gaps,
                'ss1L'       : self.length_ss1,
                'ss2L'       : self.length_ss2,
                'sequence'   : self.aminoacid_sequence,
                'structure'  : self.structure_sequence}

    def json_format(self, as_string=True):
        return self.json() if as_string else self.as_dict()

    def archtype_header(self):
        data = []
        data.append('CODE: {0._source}'.format(self))
        data.append('CHAIN: {0}'.format(self._source.split('_')[1]))
        data.append('TYPE: {0.type}'.format(self))
        data.append('ORDER: {0._order}'.format(self))
        data.append('ID_AA: {0.selfboundaries}'.format(self))
        data.append('REAL_AA: {0.boundaries}'.format(self))
        if self.is_superarch:
            data.append('INT_SS: {0._intss}'.format(self))
            data.append('INT_SS_TYPE: {0._inttyp}'.format(self))
        return "\n".join(data)

    def archtype_loop(self):
        loopinfo = []
        loopinfo.append('LOOP:')
        loopinfo.append('{0._distAA:>3d}'.format(self))
        loopinfo.append('{0._ss1._length:>3d}'.format(self))
        loopinfo.append('{0._ss2._length:>3d}'.format(self))
        loopinfo.append('{0:>3d}'.format(self._ss1.get_moment_of_inertia_length('f11')))
        loopinfo.append('{0:>3d}'.format(self._ss2.get_moment_of_inertia_length('f44')))
        loopinfo.append('{0:>.6f}'.format(self._dist))
        loopinfo.append('{0:>.6f}'.format(self._angls[0]))
        loopinfo.append('{0:>.6f}'.format(self._angls[1]))
        loopinfo.append('{0:>.6f}'.format(self._angls[2]))
        return '\t'.join(loopinfo)

    def archtype_hbond(self, matrix = True):
        data = []
        data.append('HBOND_DIM:  {0:>3d}'.format(self._hdim))
        data.append('HBOND_LOOP: {0:>3d}'.format(self._distAA))
        data.append('HBOND_NUM:  \t{0:>3d}\t{1:>3d}\t{2:>3d}\n'.format(self._htotal, self._hss, self._hloop))
        if matrix:
            data.append('HBOND_MATRIX\n\n' + self._format_hbond_matrix())
        return "\n".join(data)

    def archtype_psiphi(self):
        data = []
        psiphi = self.psiphi
        data.append('PSI_PHIA: {0}'.format(psiphi[0]))
        data.append('PSI_PHIB: {0}'.format(psiphi[1]))
        data.append('PSI_PHIC: {0}'.format(psiphi[2]))
        return "\n".join(data)

    def archtype_format(self, matrix = True):
        data = []
        data.append('\n******\n')
        data.append(self.archtype_header())
        data.append(self.archtype_loop())
        data.append('\nSEQUENCE: {0}'.format(self.aminoacid_sequence))
        data.append('TOR_CA  : {0}'.format(self._torCA))
        data.append('ACCES_SU: {0}'.format(self.access_surface))
        data.append(self.archtype_psiphi())
        data.append('SEC_STR : {0}'.format(self.structure_sequence))
        data.append('TOR_THET: {0}'.format(self._torgls[0]))
        data.append('TOR_RHO : {0}'.format(self._torgls[1]))
        data.append('ANG_DELT: {0}'.format(self._torgls[2]))
        data.append('DISTANCE: {0}\n'.format(self._tordis))
        data.append(self.archtype_hbond(matrix))
        return "\n".join(data)

    def line_format(self, sequence = False):
        data = [self.identifier, self.type, str(self.aminoacid_distance)]
        data.append('{0:07.3f}'.format(self.cartesian_distance))
        data.append('{0:07.3f}'.format(self.theta))
        data.append('{0:07.3f}'.format(self.rho))
        data.append('{0:07.3f}'.format(self.delta))
        data.append(str(self.gaps))
        if sequence:
            data.extend([self.aminoacid_sequence, self.structure_sequence])
        return "\t".join(data)

    def __lt__(self, other):
        return self._distAA < other._distAA

    def __le__(self, other):
        return self.__lt__(other) or self.__eq__(other)

    def __eq__(self, other):
        return self._distAA == other._distAA

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        return self._distAA > other._distAA

    def __ge__(self, other):
        return self.__gt__(other) or self.__eq__(other)

    def __repr__(self):
        data = []
        data.append('**Secondary Structure Relation:')
        data.append('STR1:\t{0}'.format(self._ss1.strdata('f11')))
        data.append('STR2:\t{0}'.format(self._ss2.strdata('f44')))
        data.append('ARCH TYPE: {0.type}'.format(self))
        data.append('INTERNAL SS: {0._intss}'.format(self))
        data.append('DISTANCE: {0._dist}'.format(self))
        data.append('THETA: {0.theta:<.5f} RHO: {0.rho:<.5f} DELTA: {0.delta:<.5f}'.format(self))
        data.append('**\n\n')
        return "\n".join(data)
