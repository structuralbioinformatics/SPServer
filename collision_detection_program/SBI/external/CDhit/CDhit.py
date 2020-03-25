'''
@file: CDhit.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2013

@ [oliva's lab](http://sbi.imim.es)

@class: CDhit.py
'''


class CDhit(object):
    '''
    Represents a cd-hit cluster.

    '''
    def __init__(self, cluster_id):
        '''
        @param:    cluster_id
        @pdef:     identifier of the cd-hit cluster
        @ptype:    {String}
        '''
        self._id        = cluster_id
        self._master    = None
        self._sequences = {}

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def identifier(self):
        '''
        cluster identifier.

        @return: {String}
        '''
        return self._id

    @property
    def master(self):
        '''
        Name of the master sequence of the cluster.

        @return: {String}
        '''
        return self._master

    @property
    def sequences(self):
        '''
        list of the proteins belonging to the cluster.

        @return: {List}
        '''
        return self._sequences

    ############
    # BOOLEANS #
    ############
    def is_master(self, seq):
        '''
        Checks if a given sequence is the master protein.

        @param:    seq
        @pdef:     query sequence id
        @ptype:    {String}

        @return: {Boolean}
        '''
        return self._master.name == seq

    def is_sequence(self, seq):
        '''
        Checks if a given sequence belongs to the cluster

        @param:    seq
        @pdef:     query sequence id
        @ptype:    {String}

        @return: {Boolean}
        '''
        return seq in self._sequences

    ###########
    # METHODS #
    ###########
    def add_sequence(self, cdhit_homolog):
        '''
        Adds a new cd-hit homolog to the cd-hit cluster.

        @param:    cdhit_homolog
        @pdef:     new homolog to the cluster
        @ptype:    {CDhitHomolog}
        '''
        if cdhit_homolog.is_master:
            self._master = cdhit_homolog
        else:
            self._sequences[cdhit_homolog.name] = cdhit_homolog

    def __repr__(self):
        text = []
        text.append('Cluster {0.identifier}:'.format(self))
        text.append('\tMaster Sequence: {0.master}'.format(self))
        for s in self.sequences:
            text.append('\t\t{0}'.format(self.sequences[s]))
        return '\n'.join(text)
