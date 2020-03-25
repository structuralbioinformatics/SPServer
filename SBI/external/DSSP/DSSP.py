'''
@file: DSSP.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   03/2013

@ [oliva's lab](http://sbi.imim.es)

@class: DSSP
'''
from SBI.data import aminoacids3to1     as a3to1
from SBI.data import aminoacids_surface as asurf


class DSSP(object):
    '''
    Stores a DSSP prediction for a single amino acid.

    '''
    _EXPOSED_THRESHOLD = 2.5
    # Why? DSSP exposition refers only to secondary chain (aprox.) see PMID: 15768403

    def __init__(self, secondary_structure, accessibility, amino):
        '''
        @param:    secondary_structure
        @pdef:     1 letter code of secondary structure
        @ptype:    {String}

        @param:    accessibility
        @pdef:     numeric value of accessibility of the residue
        @ptype:    {Integer}

        @param:    amino
        @pdef:     aminoacid name (1 or 3 letter code)
        @ptype:    {String}

        '''
        self._ss         = secondary_structure
        self._access     = int(accessibility)
        self._type       = amino if len(amino) == 1 else a3to1[amino]
        if amino in asurf:
            self._access10   = int((10 * self._access) / asurf[amino])
            self._accesscode = DSSP._codifyaccess((float(self._access) / asurf[amino]) * 100)
        else:
            self._access10   = 1
            self._accesscode = 1
        self._exposed    = self._access10 > DSSP.exposition_threshold
        self._rnhoa      = None
        self._enhoa      = None
        self._rohna      = None
        self._eohna      = None
        self._rnhob      = None
        self._enhob      = None
        self._rohnb      = None
        self._eohnb      = None

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def secondary_structure(self):
        '''
        @return: {String}
        '''
        return self._ss

    @property
    def aminoacid(self):
        '''
        @return: {String}
        '''
        return self._type

    @property
    def accessibility(self):
        '''
        @return: {Integer}
        '''
        return self._access

    @property
    def accessibility10(self):
        '''
        @return: {Integer}
        '''
        return self._access10

    @property
    def accesscode(self):
        '''
        @return: {String}
        '''
        return str(self._accesscode)

    @staticmethod
    @property
    def exposition_threshold():
        '''
        @return: {Float}
        '''
        return DSSP._EXPOSED_THRESHOLD

    @property
    def exposed(self):
        '''
        @return: {Boolean}
        '''
        return self._exposed

    @property
    def nhoa(self):
        '''
        @return: {List} of {Integer} and {Float}
        '''
        return self._rnhoa, self._enhoa

    @property
    def ohna(self):
        '''
        @return: {List} of {Integer} and {Float}
        '''
        return self._rohna, self._eohna

    @property
    def nhob(self):
        '''
        @return: {List} of {Integer} and {Float}
        '''
        return self._rnhob, self._enhob

    @property
    def ohnb(self):
        '''
        @return: {List} of {Integer} and {Float}
        '''
        return self._rohnb, self._eohnb

    ###########
    # METHODS #
    ###########
    def add_hydrogen_links(self, nhoa, ohna, nhob, ohnb):
        '''
        Hydrogen link data from the DSSP file

        @param:    nhoa
        @pdef:     N-HO link CA
        @ptype:    {List}

        @param:    ohna
        @pdef:     OH-N link CA
        @ptype:    {List}

        @param:    nhob
        @pdef:     N-HO link CB
        @ptype:    {List}

        @param:    ohnb
        @pdef:     OH-N link CB
        @ptype:    {List}
        '''
        self._rnhoa = int(nhoa.split(',')[0].strip())
        self._enhoa = float(nhoa.split(',')[1].strip())
        self._rohna = int(ohna.split(',')[0].strip())
        self._eohna = float(ohna.split(',')[1].strip())
        self._rnhob = int(nhob.split(',')[0].strip())
        self._enhob = float(nhob.split(',')[1].strip())
        self._rohnb = int(ohnb.split(',')[0].strip())
        self._eohnb = float(ohnb.split(',')[1].strip())

    ###################
    # PRIVATE METHODS #
    ###################
    @staticmethod
    def _codifyaccess(value):
        '''
        Codifies the access value into a numeric code.

        @param:    value
        @pdef:     actual accessibility
        @ptype:    {Float}

        @return: {String}
        '''
        value = float(value)
        if value == 0:
            return '*'  # Buried
        elif value > 100:
            return '?'  # PTM that make the Aa bigger
        else:
            if value > 0  and value <= 10:
                return '1'
            elif value > 10 and value <= 20:
                return '2'
            elif value > 20 and value <= 30:
                return '3'
            elif value > 30 and value <= 40:
                return '4'
            elif value > 40 and value <= 50:
                return '5'
            elif value > 50 and value <= 60:
                return '6'
            elif value > 60 and value <= 70:
                return '7'
            elif value > 70 and value <= 80:
                return '8'
            elif value > 80 and value <= 90:
                return '9'
            elif value > 90 and value <= 100:
                return '#'
