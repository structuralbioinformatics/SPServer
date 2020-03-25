from SBI.data import crystal_method_has_resolution
from . import process_EXPERIMENT_line, process_RESOLUTION_line
from . import process_RFACTOR_line, process_FREER_line
from SBI.beans.JSONer import JSONer


class Experiment(JSONer):
    def __init__(self, line):
        self._xpdta      = process_EXPERIMENT_line(line)
        self._res        = None
        self._rfactor    = None
        self._freeR      = None

    #
    # ATTRIBUTES
    #
    @property
    def xpdta(self):
        return self._xpdta

    @property
    def resolution(self):
        return float(self._res) if self._res is not None else -1

    @property
    def rfactor(self):
        if self._rfactor is None:
            if self.has_resolution:
                return float(0)
            else:
                return 'NULL'
        return float(self._rfactor)

    @property
    def freeR(self):
        if self._freeR is None:
            if self.has_resolution:
                return float(0)
            else:
                return 'NULL'
        return float(self._freeR)

    #
    # BOOLEANS
    #
    @property
    def has_resolution(self):
        return self.xpdta in crystal_method_has_resolution

    #
    # READ FUNCTIONS
    #
    def update_experiment(self, line):
        self._xpdta += ' ' + process_EXPERIMENT_line(line)

    def add_resolution(self, line):
        self._res = process_RESOLUTION_line(line)

    def add_rfactor(self, line):
        if self.has_resolution and self._rfactor is None:
            self._rfactor = process_RFACTOR_line(line)

    def add_freer(self, line):
        if self.has_resolution and self._freeR is None:
            self._freeR = process_FREER_line(line)

    #
    # OUT FUNCTIONS
    #
    def as_dict(self):
        nobj = {}
        nobj['xpdta']      = self.xpdta
        nobj['resolution'] = self.resolution
        nobj['rfactor']    = self.rfactor
        nobj['freer']      = self.freeR
        return nobj

    def __repr__(self):
        return self.json()
