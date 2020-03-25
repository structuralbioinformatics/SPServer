from SBI.beans.IndexedNum import IndexedNum
from SBI.beans.JSONer     import JSONer


class SecondaryStructureInfo(JSONer):
    # http://www.wwpdb.org/documentation/format32/sect5.html

    _sstype = '-'

    def __init__(self, line):
        self._line = line

    #
    # ATTRIBUTES
    #
    @property
    def line(self):
        return self._line

    @property
    def chain(self):
        raise NotImplemented

    @property
    def init(self):
        raise NotImplemented

    @property
    def end(self):
        raise NotImplemented

    @property
    def sstype(self):
        return self._sstype

    #
    # FUNCTIONS
    #
    def change_chain(self, chain):
        raise NotImplemented

    def __str__(self):
        return self.line

    #
    # OUT FUNCTIONS
    #
    def as_dict(self):
        nobj = {}
        nobj['type']  = self.sstype
        nobj['ini']   = int(self.init)
        nobj['idxs']  = self.init.index
        nobj['end']   = int(self.end)
        nobj['idxe']  = self.end.index
        nobj['chain'] = self.chain
        nobj['line']  = self.line
        return nobj


class HelixInfo(SecondaryStructureInfo):

    _sstype = 'H'
    # pdb2v7b -> I've been unable to find the correct def for code 0
    pdb_helix_def = {'0' : 'Right-handed alpha',
                     '1' : 'Right-handed alpha',
                     '2' : 'Right-handed omega',
                     '3' : 'Right-handed pi',
                     '4' : 'Right-handed gamma',
                     '5' : 'Right-handed 310',
                     '6' : 'Left-handed alpha',
                     '7' : 'Left-handed omega',
                     '8' : 'Left-handed gamma',
                     '9' : '27 ribbon/helix',
                     '10': 'Polyproline'}

    def __init__(self, line):
        super(HelixInfo, self).__init__(line)
        self._code = line[38:40].strip()

    #
    # ATTRIBUTES
    #
    @property
    def chain(self):
        return self.line[19]

    @property
    def init(self):
        return IndexedNum(self.line[20:27].strip())

    @property
    def end(self):
        return IndexedNum(self.line[32:39].strip())

    @property
    def helix_type(self):
        return self.pdb_helix_def[str(self._code)]

    #
    # FUNCTIONS
    #
    def change_chain(self, chain):
        '''
                 1         2         3         4...7
        1234567890123456789012345678901234567890...0123456
        HELIX    1  HA GLY A   86  GLY A   94  1...      9
        HELIX    2  HB GLY B   86  GLY B   94  1...      9
        '''
        string_line     = list(self.line)
        string_line[19] = chain
        string_line[31] = chain
        self._line      = "".join(string_line)

    def defines(self, chain, pos1, pos2):
        defpos1  = self.line[20:27].strip()
        defpos2  = self.line[32:38].strip()

        defchain = self.line[19:21].strip()
        if chain != defchain:
            return False

        try:
            defpos1int = int(defpos1)
        except:
            defpos1int = int(defpos1[:-1])
        try:
            defpos2int = int(defpos2)
        except:
            defpos2int = int(defpos2[:-1])
        try:
            pos1int = int(pos1)
        except:
            pos1int = int(pos1[:-1])
        try:
            pos2int = int(pos2)
        except:
            pos2int = int(pos2[:-1])

        return (pos1int >= defpos1int and pos1int <= defpos2int) or \
               (pos2int >= defpos1int and pos2int <= defpos2int) or \
               (pos1int <= defpos1int and pos2int >= defpos2int)


class SheetInfo(SecondaryStructureInfo):

    _sstype = 'E'

    def __init__(self, line):
        super(SheetInfo, self).__init__(line)

    #
    # ATTRIBUTES
    #
    @property
    def chain(self):
        return self.line[21]

    @property
    def init(self):
        return IndexedNum(self.line[22:28].strip())

    @property
    def end(self):
        return IndexedNum(self.line[33:38].strip())

    @property
    def beta_direction(self):
        return self.line[38:40].strip()

    def change_chain(self, chain):
        '''
                 1         2         3         4         5         6         7
        1234567890123456789012345678901234567890123456789012345678901234567890
        SHEET    1   A 5 THR A 107  ARG A 110  0
        SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
        SHEET    3   A 5 ARG A  87  SER A  91 -1  N  LEU A  89   O  TYR A  97
        SHEET    4   A 5 TRP A  71  ASP A  75 -1  N  ALA A  74   O  ILE A  88
        SHEET    5   A 5 GLY A  52  PHE A  56 -1  N  PHE A  56   O  TRP A  71
        '''
        string_line     = list(self.line)
        string_line[21] = chain
        string_line[32] = chain

        if (len(string_line) > 45 and string_line[49] != ' '):
            string_line[49] = chain
            string_line[64] = chain

        self._line      = "".join(string_line)

    def defines(self, chain, pos1, pos2):
        defpos1 = self.line[22:28].strip()
        defpos2 = self.line[33:38].strip()

        defchain = self.line[21:22].strip()
        if chain != defchain:
            return False

        try:
            defpos1int = int(defpos1)
        except:
            defpos1int = int(defpos1[:-1])
        try:
            defpos2int = int(defpos2)
        except:
            defpos2int = int(defpos2[:-1])
        try:
            pos1int = int(pos1)
        except:
            pos1int = int(pos1[:-1])
        try:
            pos2int = int(pos2)
        except:
            pos2int = int(pos2[:-1])

        decission = (pos1int >= defpos1int and pos1int <= defpos2int) or \
                    (pos2int >= defpos1int and pos2int <= defpos2int) or \
                    (pos1int <= defpos1int and pos2int >= defpos2int)

        double = (len(self.line.strip()) > 45 and self.line[39] == '1')

        if (double and not decission):
            defpos1 = self.line[50:56].strip()
            defpos2 = self.line[65:].strip()
            try:
                defpos1int = int(defpos1)
            except:
                defpos1int = int(defpos1[:-1])
            try:
                defpos2int = int(defpos2)
            except:
                defpos2int = int(defpos2[:-1])

            decission = (pos1int >= defpos1int and pos1int <= defpos2int) or \
                        (pos2int >= defpos1int and pos2int <= defpos2int) or \
                        (pos1int <= defpos1int and pos2int >= defpos2int)

        return decission


class TurnInfo(SecondaryStructureInfo):

    _sstype = 'C'

    def __init__(self, line):
        super(TurnInfo, self).__init__(line)

    #
    # ATTRIBUTES
    #
    @property
    def chain(self):
        return self.line[19]

    @property
    def init(self):
        return IndexedNum(self.line[20:26].strip())

    @property
    def end(self):
        return IndexedNum(self.line[31:39].strip())

    @property
    def turn_type(self):
        return self.line[40:].strip()

    #
    # FUNCTIONS
    #
    def change_chain(self, chain):
        '''
                 1         2         3         4         5
        12345678901234567890123456789012345678901234567890
        TURN     1 S1A GLY A  16  GLN A  18     SURFACE
        TURN     2 FLA ILE A  50  GLY A  52     FLAP
        TURN     3 S2A ILE A  66  HIS A  69     SURFACE
        TURN     4 S1B GLY B  16  GLN B  18     SURFACE
        TURN     5 FLB ILE B  50  GLY B  52     FLAP
        TURN     6 S2B ILE B  66  HIS B  69     SURFACE
        '''
        string_line     = list(self.line)
        string_line[19] = chain
        string_line[30] = chain
        self._line      = "".join(string_line)

    def defines(self, chain, pos1, pos2):
        defpos1 = self.line[22:26].strip()
        defpos2 = self.line[31:37].strip()

        defchain = self.line[19:21].strip()
        if chain != defchain:
            return False

        try:
            defpos1int = int(defpos1)
        except:
            defpos1int = int(defpos1[:-1])
        try:
            defpos2int = int(defpos2)
        except:
            defpos2int = int(defpos2[:-1])
        try:
            pos1int = int(pos1)
        except:
            pos1int = int(pos1[:-1])
        try:
            pos2int = int(pos2)
        except:
            pos2int = int(pos2[:-1])

        return (pos1int >= defpos1int and pos1int <= defpos2int) or \
               (pos2int >= defpos1int and pos2int <= defpos2int) or \
               (pos1int <= defpos1int and pos2int >= defpos2int)
