'''
@file: HmmHit.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2014

@ [oliva's lab](http://sbi.imim.es)

@class: HmmHit
'''


class HmmHit(object):
    '''
    Stores the data of each hit output of HMMER
    '''
    def __init__(self, name, accession, length, query_lim, hit_lim,
                 values, repeats):
        '''
        @param:    name
        @pdef:     hit name
        @ptype:    {String}

        @param:    accession
        @pdef:     hit identifier
        @ptype:    {String}

        @param:    length
        @pdef:     length of the query protein
        @ptype:    {Integer}

        @param:    query_lim
        @pdef:     limits of the assignation in the query sequence
        @ptype:    {List} of {Integer}

        @param:    hit_lim
        @pdef:     limits of the assignation in the hit sequence
        @ptype:    {List} of {Integer}

        @param:    values
        @pdef:     evalue and ievalue
        @ptype:    {List} of {Float}

        @param:    repeats
        @pdef:     number of apparitions of this hit in the sequence
        @ptype:    {Integer}
        '''
        self._accession = accession if accession != '-' else name.split(':')[1]
        self._name      = name if accession != '-'  else name.split(':')[0]
        self._length    = length
        self._qini      = query_lim[0]
        self._qend      = query_lim[1]
        self._alength   = int(self._qend) - int(self._qini) + 1
        self._hini      = hit_lim[0]
        self._hend      = hit_lim[1]
        self._evalue    = values[0]
        self._repeats   = repeats
        self._ievalue   = values[1]
        self._coverage  = self._set_coverage()

        self._used      = False

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def accession(self):
        '''
        @return: {String}
        '''
        return self._accession

    @property
    def name(self):
        '''
        @return: {String}
        '''
        return self._name

    @property
    def length(self):
        '''
        @return: {Integer}
        '''
        return int(self._length)

    @property
    def applied_length(self):
        '''
        @return: {Integer}
        '''
        return int(self._alength)

    @property
    def query_ini(self):
        '''
        @return: {Integer}
        '''
        return int(self._qini)

    @property
    def query_end(self):
        '''
        @return: {Integer}
        '''
        return int(self._qend)

    @property
    def hit_ini(self):
        '''
        @return: {Integer}
        '''
        return int(self._hini)

    @property
    def hit_end(self):
        '''
        @return: {Integer}
        '''
        return int(self._hend)

    @property
    def evalue(self):
        return float(self._evalue)

    @property
    def repeats(self):
        '''
        @return: {Integer}
        '''
        return int(self._repeats)

    @property
    def ievalue(self):
        '''
        @return: {Float}
        '''
        return float(self._ievalue)

    @property
    def coverage(self):
        '''
        @return: {Float}
        '''
        return float(self._coverage)

    @property
    def used(self):
        '''
        @return: {Boolean}
        '''
        return self._used

    @used.setter
    def used(self, value):
        '''
        @param:    value
        @pdef:     new used status
        @ptype:    {Boolean}
        '''
        self._used = value

    ###########
    # METHODS #
    ###########
    def overlap(self, domain):
        '''
        calculates the overlap of two domains over the same sequence

        @param:    domain
        @pdef:     domain to which compare
        @ptype:    {HmmHit}

        @return: {Float}
        '''
        if domain.query_ini >= self.query_ini and domain.query_end <= self.query_end:
            return 1
        if domain.query_end <= self.query_ini or domain.query_ini > self.query_end:
            return 0

        start = max(domain.query_ini, self.query_ini)
        end   = min(domain.query_end, self.query_end)

        return float(end-start)/float(domain.applied_length)

    ###################
    # PRIVATE METHODS #
    ###################
    def _set_coverage(self):
        '''
        calculates how much of the domain is assigned to the protein

        @return: {Float}
        '''
        return float(self.hit_end - self.hit_ini + 1) / float(self.length)
