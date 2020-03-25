import numpy as np
from SBI.beans.JSONer import JSONer


class BioMolecule(JSONer):
    def __init__(self, identifier):
        self._identifier = identifier
        self._matrices   = []
        self._chains     = set()

    #
    # ATTRIBUTES
    #
    @property
    def identifier(self):
        return self._identifier

    @property
    def chains(self):
        return self._chains

    @chains.setter
    def chains(self, value):
        if isinstance(value, list):
            self._chains.update(set(value))
        else:
            self._chains.add(value)

    @property
    def matrices(self):
        if len(self._matrices) > 0:
            return self._matrices
        else:
            return [MatrixAndVector.stillmatrix(), ]

    #
    # FUNCTIONS
    #
    def new_matrix(self):
        self._matrices.append(MatrixAndVector.stillmatrix())

    def update_last_matrix(self, row, mx, my, mz, v):
        self._matrices[-1].update_matrix(row, mx, my, mz)
        self._matrices[-1].update_vector(row, v if v != '' else 0)

    #
    # OUTPUT FUNCTIONS
    #
    def as_dict(self):
        data = {'id':       self.identifier,
                'chains':   list(self.chains),
                'matrices': [x.as_dict() for x in self.matrices]}

        return data

    #
    # OVERWRITE DEFAULT FUNCTIONS
    #
    def __len__(self):
        return len(self._matrices)


class MatrixAndVector(object):
    def __init__(self, matrix, vector):
        self._matrix = matrix
        self._vector = vector

    #
    # ATTRIBUTES
    #
    @property
    def matrix(self):
        return self._matrix

    @property
    def vector(self):
        return self._vector

    #
    # FUNCTIONS
    #
    @staticmethod
    def stillmatrix():
        return MatrixAndVector(np.identity(3, float),
                               np.zeros(3, float))

    def update_matrix(self, row, mx, my, mz):
        try:  # 3DPA
            self._matrix[int(row)-1][0] = float(mx)
        except ValueError:
            self._matrix[int(row)-1][0] = float(mx.split()[-1])
        try:
            self._matrix[int(row)-1][1] = float(my)
        except ValueError:
            self._matrix[int(row)-1][1] = float(my.split()[-1])
        try:
            self._matrix[int(row)-1][2] = float(mz)
        except ValueError:
            self._matrix[int(row)-1][2] = float(mz.split()[-1])

    def update_vector(self, row, v):
        self._vector[int(row)-1] = float(v)

    #
    # OUTPUT FUNCTIONS
    #
    def as_dict(self):
        data = {'matrix': self.matrix.tolist(),
                'vector': self.vector.tolist()}
        return data
