'''
@file: CDhitList.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2013

@ [oliva's lab](http://sbi.imim.es)

@class: CDhitList
'''
import re

from SBI.beans      import StorableObject
from SBI.beans      import File
from .CDhit         import CDhit
from .CDhitHomolog  import CDhitHomolog


class CDhitList(StorableObject):
    '''
    List of cd-hit clusters.

    '''
    def __init__(self, cdhit_file = None):
        '''
        @param:    cdhit_file
        @pdef:     name of the cd-hit output file
        @pdefault: _None_. Create an empty list
        @ptype:    {String}

        '''
        self._clusters  = []
        self._allseqids = {}
        if cdhit_file is not None:
            self._file  = File(file_name = cdhit_file)
        else:
            self._file  = None

        if self._file is not None:
            self._parse_file()

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def clusters(self):
        '''
        List of cd-hit clusters.

        @return: {List} of {CDhit}
        '''
        return self._clusters

    ###########
    # METHODS #
    ###########
    def get_cluster4sequence(self, sequence):
        '''
        Retrieve a cluster for a given sequence. _None_ if the sequence is not
        found.

        @param:    sequence
        @pdef:     name of the query sequence
        @ptype:    {String}

        @return: {CDhit}
        '''
        if sequence in self._allseqids:
            return self._clusters[self._allseqids[sequence]]
        else:
            return None

    def is_in_cluster(self, sequence):
        '''
        Evaluate if the sequence is in a cluster.

        @param:    sequence
        @pdef:     name of the query sequence
        @ptype:    {String}

        @return: {String} as 'N' if no, 'H' if yes and 'M' if cluster master
        '''
        c = self.get_cluster4sequence(sequence)
        if c is None:
            return 'N'
        else:
            return 'M' if c.is_master(sequence) else 'H'

    def add_cluster(self, cluster):
        '''
        Add a cd-hit cluster to the object.

        @param:    cluster
        @pdef:     new cd-hit cluster to add
        @ptype:    {CDhit}
        '''
        self._clusters.append(cluster)

    def add_sequence2cluster(self, sequence, cluster_id = None):
        '''
        Add a new sequence to a given cluster.

        @param:    sequence
        @pdef:     name of the query sequence
        @ptype:    {String}

        @param:    cluster_id
        @pdef:     identifier of the cluster
        @pdefault: _None_. Refers to the last added cluster.
        @ptype:    {String}
        '''
        if cluster_id is None:
            self.clusters[-1].add_sequence(sequence)
            self._allseqids[sequence.name] = len(self.clusters) - 1
        else:
            for x in range(len(self._clusters)):
                if self._clusters[x].identifier == cluster_id:
                    self._clusters[x].add_sequence(sequence)
                    self._allseqids[sequence.name] = x
                    break

    def dictionary_role_summary(self):
        '''
        Creates a dictionary separating master sequences and homolog sequences.

        @return: {Dictionary}
        '''
        data = {'master': [], 'homolog': []}
        for c in self.clusters:
            data['master'].append(c.master.name)
            for s in c.sequences:
                data['homolog'].append(s)
        return data

    def merge_clusters(self, cluster_file):
        '''
        When using an intermediate state to cluster by homology,
        the result of the second clustering is a clustering of clusters.
        We need to transform this into the original sequences

        @param:    cluster_file
        @pdef:     name of the second-step cluster output
        @ptype:    {String}
        '''
        clustlist  = CDhitList(cluster_file)
        newlist    = CDhitList()
        cluster_re = re.compile('Cluster\s+(\d+)')
        for cl in clustlist.clusters:
            c = CDhit(cluster_id = cl.identifier)
            newlist.add_cluster(c)
            cnum = int(cluster_re.search(cl.master.name).group(1))
            oldclust = self.clusters[cnum]
            newlist.add_sequence2cluster(sequence = oldclust.master)
            for s in oldclust.sequences:
                newlist.add_sequence2cluster(sequence = oldclust.sequences[s])
            for s in cl.sequences:
                idclust = cl.sequences[s]
                cnum = int(cluster_re.search(idclust.name).group(1))
                oldclust = self.clusters[cnum]
                master   = oldclust.master
                master.homology = idclust.homology
                newlist.add_sequence2cluster(sequence = master)
                for s in oldclust.sequences:
                    h = oldclust.sequences[s]
                    h.homology = int(h.homology * float(idclust.homology)/10)
                    newlist.add_sequence2cluster(sequence = h)

        self._clusters  = newlist._clusters
        self._allseqids = newlist._allseqids

    ###################
    # PRIVATE METHODS #
    ###################
    def _parse_file(self):
        '''
        Read the cd-hit output file into a {CDhitList}

        '''
        homolog_re = re.compile('(\d+)aa,\s+\>([\s\w]+)\.{3}')
        for line in self._file.read():
            if line.startswith('>'):
                c = CDhit(cluster_id = line.split()[-1].strip())
                self.add_cluster(c)
            else:
                data = homolog_re.search(line)
                d = line.split()
                h    = CDhitHomolog(name     = data.group(2),
                                    length   = data.group(1),
                                    homology = d[-1])
                self.add_sequence2cluster(sequence = h)
        self._file.close()

    def __len__(self):
        return len(self._clusters)

    def __repr__(self):
        text = []
        for c in self.clusters:
            text.append('{0}'.format(c))
        return '\n'.join(text)
