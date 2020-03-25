'''
@file: HmmResult.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2014

@ [oliva's lab](http://sbi.imim.es)

@class: HmmResult
'''
import re

from SBI.beans import StorableObject


class HmmResult(StorableObject):
    '''
    Contains the results of a HMMER output file.

    '''
    _PFAM = re.compile('pfam', re.IGNORECASE)
    _CATH = re.compile('cath', re.IGNORECASE)

    def __init__(self, query_name, query_sequence, database = None):
        '''
        @param:    query_name
        @pdef:     name of the query protein
        @ptype:    {String}

        @param:    query_sequence
        @pdef:     sequence of the query protein
        @ptype:    {String}

        @param:    database
        @pdef:     name of the queried database
        @ptype:    {String}
        '''
        self._qname = query_name
        self._qseq  = query_sequence
        self._db    = database
        self._hits  = []

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def query_name(self):
        '''
        @return: {String}
        '''
        return self._qname

    @property
    def query_sequence(self):
        '''
        @return: {String}
        '''
        return self._qseq

    @property
    def database(self):
        '''
        @return: {String}
        '''
        return self._db

    @property
    def hits(self):
        '''
        @return: {List} of {HmmHit}
        '''
        return self._hits

    ############
    # BOOLEANS #
    ############
    @property
    def is_pfam(self):
        '''
        Check if we queried against a PFAM database

        @return: {Boolean}
        '''
        return bool(self._PFAM.search(self._db))

    @property
    def is_cath(self):
        '''
        Check if we queried against a CATH database

        @return: {Boolean}
        '''
        return bool(self._CATH.search(self._db))

    @property
    def has_domains(self):
        '''
        Check if any hit was found.

        @return:{Boolean}
        '''
        return len(self.hits) > 0

    @property
    def is_empty(self):
        '''
        Check if no hit was found.

        @return:{Boolean}
        '''
        return not self.has_domains

    ###########
    # METHODS #
    ###########
    def add_hit(self, hit):
        '''
        append a new hmm hit result to the object.

        @param:    hit
        @pdef:     new hit to add
        @ptype:    {HmmHit}
        '''
        self._hits.append(hit)

    def get_hits(self, evalue = None, coverage = None, overlap = None):
        '''
        Retrieve hits that fulfill different filters.

        @param:    evalue
        @pdef:     maximum evalue allowed
        @ptype:    {Float}

        @param:    coverage
        @pdef:     minimum coverage of the assigned hit required
        @ptype:    {Float}

        @param:    overlap
        @pdef:     maximum overlap to previous retrieved hits allowed
        @ptype:    {Float}

        @return: {List} of {HmmHit}
        '''
        if evalue is None and coverage is None and overlap is None:
            return self.hits
        else:
            hits = []
            coverage = float(0) if coverage is None else float(coverage)
            evalue   = float(1000) if evalue is None else float(evalue)
            for domain in self.hits:
                if evalue >= domain.evalue and coverage <= domain.coverage:
                    if overlap is None:
                        hits.append(domain)
                    else:
                        addDomain = True
                        for previous_domain in hits:
                            if previous_domain.overlap(domain) > overlap:
                                addDomain = False
                        if addDomain:
                            hits.append(domain)
            return hits

    def get_used_hits(self):
        '''
        Retrieve those hits marked as used.

        @return: {List} of {HmmHit}
        '''
        hits = []
        for domain in self.hits:
            if domain.used:
                hits.append(domain)
        return hits
