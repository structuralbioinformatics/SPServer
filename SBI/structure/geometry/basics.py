import numpy as np


def vectors_angle(a, b, degrees = False):
    def _vdot(a, b):
        if len(a.shape) == 1:
            return np.sum(a*b)
        elif len(a.shape) == 2:
            return np.sum(a*b, 1)

    def _vabs(a):
        if len(a.shape) == 1:
            return np.sqrt(np.sum(a**2))
        elif len(a.shape) == 2:
            return np.sqrt(np.sum(a**2, 1))

    a = np.asarray(a)
    b = np.asarray(b)

    if a.shape != b.shape:
        raise 'vectors have different shape'

    if not degrees:
        return np.arccos(_vdot(a, b) / (_vabs(a)*_vabs(b)))
    else:
        return np.degrees(np.arccos(_vdot(a, b) / (_vabs(a)*_vabs(b))))
