class CDhitHomolog(object):
    '''
    Representative of a sequence in the cd-hit cluster.

    '''
    def __init__(self, name, length, homology):
        '''
        @param:    name
        @pdef:     name of the sequence
        @ptype:    {String}

        @param:    length
        @pdef:     length of the sequence
        @ptype:    {String}

        @param:    homology
        @pdef:     level of homology inside the cluster
        @ptype:    {String}
        '''
        self._name     = name.lstrip('>').rstrip('.')
        self._length   = int(length.replace('aa', '').rstrip(','))
        self._homology = homology.replace('%', '')

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def name(self):
        '''
        name of the sequence

        @return: {String}
        '''
        return self._name

    @property
    def length(self):
        '''
        length of the protein.

        @return: {Integer}
        '''
        return self._length

    @property
    def homology(self):
        '''
        level of homology inside the cluster

        @return: {Integer} or {String}
        '''
        if self.is_master:
            return self._homology
        else:
            return int(self._homology)

    @homology.setter
    def homology(self, value):
        '''
        @param:    value
        @pdef:     level of homology inside the cluster
        @ptype:    {String}
        '''
        self._homology = int(value)

    ############
    # BOOLEANS #
    ############
    @property
    def is_master(self):
        '''
        checks if protein is representative of the cluster.

        @return: {Boolean}
        '''
        return self._homology == '*'

    def __repr__(self):
        if not self.is_master:
            return '{0.name}: {0.length:0004d} Aa with {0.homology:003d}%'.format(self)
        else:
            return '{0.name}: {0.length:0004d} Aa.'.format(self)
