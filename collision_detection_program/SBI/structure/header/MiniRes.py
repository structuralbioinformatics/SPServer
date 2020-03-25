from SBI.beans.IndexedNum import IndexedNum
from SBI.beans.JSONer     import JSONer


class MiniResidue(JSONer):
    def __init__(self, restype, position):
        self._type  = restype
        self._chain = None
        self._pos   = None
        self._parse_position(position)

    #
    # ATTRIBUTES
    #
    @property
    def restype(self):
        return self._type

    @property
    def chain(self):
        return self._chain

    @property
    def position(self):
        return int(self._pos)

    @property
    def idxp(self):
        return self._pos.index

    def _parse_position(self, position):
        try:
            int(position)
            self._chain = 'A'
            self._pos   = IndexedNum(position)
        except ValueError:
            if position != '':  # 1DLK
                self._chain = position[0]
                self._pos   = IndexedNum(position[1:].strip())

    #
    # OUT FUNCTION
    #
    def as_dict(self):
        nobj = {}
        nobj['type']  = self.restype
        nobj['chain'] = self.chain
        nobj['pos']   = self.position
        nobj['idxp']  = self.idxp
        return nobj
