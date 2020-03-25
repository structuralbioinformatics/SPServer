'''
@file: BlastResult.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2011

@ [oliva's lab](http://sbi.imim.es)

@class: BlastResult
@class: BlastHeader
@class: HitFilter
'''
import os
import copy

import BlastExe
from SBI.beans              import StorableObject
from SBI.beans              import File
from SBI.sequence           import Sequence
from SBI.sequence.alignment import Rost
from SBI                    import SBIglobals as SBIg


class BlastResult(StorableObject):
    '''
    Contains all the information provided by a BLAST execution.
    Basically is a container of {BlastHit}.

    Theoretically, it will always be created through the {BlastExe} object.

    The length of the {BlastResult} correspond to the number of hits in the
    iteration of interest (last iteration by default)

    '''
    def __init__(self, query_name, query_sequence, header = None):
        '''
        @param:    query_name
        @pdef:     name of the query sequence
        @ptype:    {String}

        @param:    query_sequence
        @pdef:     query sequence
        @ptype:    {String}

        @param:    header
        @pdef:     main data about the blast execution
        @pdefault: _None_
        @ptype:    {BlastHeader}
        '''
        self._query  = Sequence(query_name, query_sequence)
        self._header = header
        self._filter = HitFilter()

        self._lastiteration = 0                  # Keep the last iteration
        self._hits          = []                 # List of BlastHit objects
        self._correctedHits = False

        self._filter_hits   = None
        self._query_index   = None

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def query(self):
        '''
        Name of the query sequence

        @return: {String}
        '''
        return self._query.full_identifier

    @property
    def query_length(self):
        '''
        Length of the query sequence

        @return: {Integer}
        '''
        return len(self._query)

    @property
    def query_sequence(self):
        '''
        Query sequence

        @return: {String}
        '''
        return self._query.sequence

    @property
    def query_object(self):
        '''
        {Sequence} query

        @return: {Sequence}
        '''
        return self._query

    @property
    def version(self):
        '''
        Version of the blast executable.

        @return: {String}
        '''
        return self._header.version

    @property
    def matrix(self):
        '''
        Similarity matrix used.

        @return: {String}
        '''
        return self._header.matrix

    @property
    def gap_open(self):
        '''
        Gap open penalty

        @return: {Integer}
        '''
        return self._header.gap_open

    @property
    def gap_extend(self):
        '''
        Gap extension penalty

        @return: {Integer}
        '''
        return self._header.gap_extend

    @property
    def database(self):
        '''
        Name of the database searched

        @return: {String}
        '''
        return self._header.database

    @property
    def lastiteration(self):
        '''
        Number of the last iteration.
        0 means that there are no hits found.

        @return: {Integer}
        '''
        return int(self._lastiteration)

    @property
    def hits(self):
        '''
        Blast hits contained in the {BlastResult}. According to certain
        selection parameters.

        @return: {List} of {BlastHit}
        '''
        if self._filter_hits is not None:
            return self._filter_hits

        self._filter_hits = []

        i = self._filter.iteration
        e = self._filter.evalue
        r = self._filter.tz_type
        p = self._filter.tz_parameter
        o = self._filter.overlap

        for hit in self.raw_hits:
            ros = (r is None or hit.evaluate_Rost_twilight_zone(r, p))
            ite = hit.iteration == i
            eva = (e is None or hit.evalue <= float(e))
            if ite and eva and ros:
                self._filter_hits.append(hit)
            if hit.iteration > i:  # ignore bigger iterations, if any
                break

        if o == 1 or len(self._filter_hits) == 0:
            return self._filter_hits
        else:
            return self._get_nonoverlaping_hits(self._filter_hits, o)

    @property
    def raw_hits(self):
        '''
        All {BlastHit} contained in the object.

        @return: {List} of {BlastHit}
        '''
        return self._hits

    ############
    # BOOLEANS #
    ############
    @property
    def are_hits_corrected(self):
        '''
        Checks if hits have already been corrected.

        @return: {Boolean}
        '''
        return self._correctedHits

    @property
    def is_selfhit_ignored(self):
        '''
        Checks if self hit has been ignored (if any)

        @return: {Boolean}
        '''
        return not self._header.self_hit

    @property
    def has_hits(self):
        '''
        Checks that there is some hit from the blast result

        @return: {Boolean}
        '''
        return not self.lastiteration == 0

    @property
    def is_empty(self):
        '''
        Checks that there are no hits from the blast result (no homologs)

        @return: {Boolean}
        '''
        return not self.has_hits

    ###########
    # METHODS #
    ###########
    def add_hit(self, hit):
        '''
        Add a new blast hit to the container.

        '''
        self._hits.append(hit)

    def set_hit_filter(self, iteration = None, evalue = None,
                       tz_type = None, tz_parameter = None, overlap = 1):
        '''
        Set parameters to exclude some hits that do not fulfill the filter.
        The parameters can be passed to the different functions that call the hits,
        but doing it here implies not having to do it every time.
        _IMPORTANT:_ Any non-specified parameter will be set to defaults.

        @param:    iteration
        @pdef:     query iteration. if a iteration bigger than the last iteration
                   is given, it defaults to last iteration. if 0 is given, the first
                   iteration (1) is assumed.
        @pdefault: last iteration
        @ptype:    {Integer}

        @param:    evalue
        @pdef:     maximum e-value considered.
        @pdefault: _None_, no filter by e-value
        @ptype:    {Float}

        @param:    tz_type
        @pdef:     Rost curve.
        @pdefault: _None_, no filter by Rost
        @poptions: 'ID', 'identity', 'SIM', 'similarity', 'hssp'
        @ptype:    {String}

        @param:    tz_parameter
        @pdef:     threshold parameter for Rost's curve
        @pdefault: 5 for 'ID', 'identity'; 12 for 'SIM', 'similarity';
                   8 for 'hssp'
        @pclash:   if tz_type is _None_ is ignored
        @ptype:    {Integer}

        @param:    overlap
        @pdef:     maximum amount of sequence overlap (over the query) allowed
                   between results. from 0 (no overlap) to 1 (total overlap).
                   Overlap values bigger than 1 default to 1, those smaller than
                   0 default to 0.
        @pdefault: 1, no filter by overlap.
        @ptype:    {Float}
        '''
        # TODO review overlap
        self.set_hit_filter_default()

        if iteration is not None:
            self._filter.iteration = iteration
            SBIg.alert('debug', self, 'Set iteration filter to {0}'.format(self._filter.iteration))
        if evalue is not None:
            self._filter.evalue = evalue
            SBIg.alert('debug', self, 'Set evalue filter to {0}'.format(self._filter.evalue))
        if tz_type is not None:
            self._filter.tz_type = tz_type
            SBIg.alert('debug', self, 'Activating {0} Rost'.format(self._filter.tz_type))
            if tz_parameter is not None:
                self._filter.tz_parameter = tz_parameter
                SBIg.alert('debug', self, 'Setting Rost at {0}'.format(self._filter.tz_parameter))
        if overlap != 1:
            self._filter.overlap = overlap
            SBIg.alert('debug', self, 'Set overlap filter to {0}'.format(self._filter.overlap))

    def set_hit_filter_default(self):
        '''
        Return the hit filter to default parameters: NO FILTER
        '''
        self._filter_hits = None
        self._filter.lastiteration = self._lastiteration
        self._filter.set_default()

    def get_hit(self, hit_ID):
        '''
        Get a hit by hit sequence ID.

        @param:    hit_ID
        @pdef:     name of the hit
        @ptype:    {String}
        '''
        for h in self.hits:
            if h.sequenceID == hit_ID:
                return h
        return None

    def get_used_hits(self):
        '''
        Retrieves all the {BlastHit} that have been marked as used.

        @return: {List} of {BlastHit}
        '''
        hits_of_interest = []
        for hit in self._hits:
            if hit.is_used:
                hits_of_interest.append(hit)
        return hits_of_interest

    def correct_hit_count(self, count_hit_file = None, count_query_file = None,
                          return_correction_dict = False):
        '''
        Corrects the starting point of the hits and the query, if needed.
        Why?
        When blasting vs. PDB (for example), sometimes the hit positions given
        by blast are wrong, as the blast always consider the first position of
        the hit sequence as 1 and PDB does not.
        Even more, the position reference doesn't even need to be a number.
        As the specific location in the PDB is important, we need to adapt our
        blasts so than we can read that data.
        Keep in mind that hits and query must be corrected together in this step,
        as this function cannot be called twice for a same instance.

        @param:    count_hit_file
        @pdef:     file containing the idex data for the query database
                   each sequence in this file will have a format such as:
                   >3K2K_A -7 ;-6 ;-5 ;-4 ;-3 ;-2 ;-1 ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ...
        @ptype:    {String}

        @param:    count_query_file
        @pdef:     sometimes we might also need to correct the query (if PDB vs.
                   PDB). Same format as count_hit_file. They might be the same file.
        @ptype:    {String}

        @param:    return_correction_dict
        @pdef:     instead of actually executing the correction, it only returns
                   the dictionary for further use.
        @pdefault: _False_
        @ptype:    {Boolean}

        @raises: {IOError} if the correction index file does not exist.
        @raises: {AttributeError} if the BlastResult does not contain any BlastHit.
        @raises: {BlastError} if it has been called before for this instance.

        '''
        if not self.has_hits:
            SBIg.warn(self, "BlastResult of {0} has no hits to correct".format(self.query))
            return

        if self.are_hits_corrected:
            be = BlastExe.BlastError()
            raise be.corrected_hits()

        SBIg.alert('debug', self, 'Correcting indexes for {0}'.format(self.query))
        cfile = File(count_hit_file)
        cq    = False

        codes_of_interest = set([hit.sequenceID for hit in self.raw_hits])
        if count_query_file == count_hit_file:
            codes_of_interest.add(self.query)
            count_query_file = None
            cq               = True

        start_index_dic = {}
        for line in cfile.read():
            if len(line.strip()) > 0:
                k = line.split('\t')
                if k[0].lstrip('>') in codes_of_interest:
                    start_index_dic[k[0].lstrip('>')] = k[1].strip().split(';')
        cfile.close()

        if count_query_file is not None:
            cfile = File(count_query_file)
            for line in cfile.read():
                if len(line.strip()) > 0:
                    k = line.split('\t')
                    if k[0].lstrip('>') == self.query:
                        start_index_dic[k[0].lstrip('>')] = k[1].strip().split(';')
            cfile.read().close()
            cq = True

        if cq:
            SBIg.alert('debug', self, '\tFixing Query {0}'.format(self.query))
            self._query_index = start_index_dic[self.query]

        if return_correction_dict:
            return start_index_dic

        for hit in self._hits:
            # This tests between the options PDB/PDB_ID or PDB_ID in case
            # the TAB file has different codification
            h      = hit.sequenceID
            hit_ID = h if h in start_index_dic else h.split("/")[-1]
            SBIg.alert('debug', self, '\tFixing {0}'.format(hit_ID))
            hit.correct_hit_count(new_index = start_index_dic[hit_ID])
            if cq:
                SBIg.alert('debug', self, '\tFixing Query {0}'.format(self.query))
                hit.correct_query_count(new_index = start_index_dic[self.query])

        self._correctedHits = True

    def set_last_iteration(self):
        '''
        Modifies the value of the last iteration with the value of iteration
        of the last BlastHit found

        '''
        self._lastiteration = 0 if self._hits == [] else self._hits[-1].iteration
        self._header.lastiteration = self._lastiteration

    def same_blast_version(self, blastObject):
        '''
        Are the two blasts the same executable?

        @param:    blastObject
        @pdef:     object against which compare
        @ptype:    {BlastResult}
        '''
        return self._header.same_version(blastObject._header)

    def same_blast_database(self, blastObject):
        '''
        Do the two blast use the same database?

        @param:    blastObject
        @pdef:     object against which compare
        @ptype:    {BlastResult}
        '''
        return self._header.same_database(blastObject._header)

    def same_blast_matrix(self, blastObject):
        '''
        Have two blast the same similarity matrix?

        @param:    blastObject
        @pdef:     object against which compare
        @ptype:    {BlastResult}
        '''
        return self._header.same_matrix(blastObject._header)

    def same_blast_iterations(self, blastObject):
        '''
        Have two blast the same number of iterations?

        @param:    blastObject
        @pdef:     object against which compare
        @ptype:    {BlastResult}
        '''
        return self.lastiteration == blastObject.lastiteration

    def same_blast_gap_properties(self, blastObject):
        '''
        Have two blast the same gap opening/extension policy?

        @param:    blastObject
        @pdef:     object against which compare
        @ptype:    {BlastResult}
        '''
        return self._header.same_gap_properties(blastObject._header)

    def comparable(self, blastObject):
        '''
        Are two blasts comparable?

        @param:    blastObject
        @pdef:     object against which compare
        @ptype:    {BlastResult}
        '''
        return self._header == blastObject._header

    def duplicate(self, new_query_sequence = None):
        '''
        Returns a copy of the object

        @param:    new_query_sequence
        @pdef:     New sequence to substitute the source query
        @ptype:    {Sequence}

        @return: {BlastResult}
        '''
        b = copy.deepcopy(self)
        if new_query_sequence is not None:
            b._query = new_query_sequence
        return b

    @staticmethod
    def build_fileobject_name(query_sequence, directory):
        '''
        Constructs a standardized name for the file to store the fileobject in.

        @param:    query_sequence
        @pdef:     query sequence
        @ptype:    {Sequence}

        @param:    directory
        @pdef:     global output directory
        @ptype:    {String}

        @return: {String}
        '''
        blastid  = '.'.join([query_sequence.md5, query_sequence.id, 'blastObj.gz'])
        blastdir = os.path.join(directory, os.path.join('blastOut', blastid[2:5]))
        return os.path.join(blastdir, blastid)

    def str_compacted_blast(self):
        '''
        Create a compacted representation of the alignment.

        @return: {String}
        '''
        lines = self.str_blast_header()
        lines.append("#Query Sequence:        {0}".format(self.query_sequence))
        for hit in self.hits:
            lines.append("\t".join((self.query, str(self.query_length), str(hit))))
        return "\n".join(lines)

    def str_PIR(self, result = 0):
        '''
        Create a specific hit as a PIR alignment.

        @param:    result
        @pdef:     numer-identifier of the {BlastHit} of interest
        @pdefault: 0. First hit.
        @ptype:    {Integer}
        '''
        lines = []
        hit = self.hits[result]

        #Counting in the possibility that the ID of the PDB hit is:
        #    PDB/2ERD_A
        hit_name  = hit.sequenceID.replace("/", "|").upper()
        hit_chain = hit_name.split("_")[-1]
        hit_id    = hit_name.split("|")[-1].split("_")[0]
        hit_sq    = hit.hit_seq + "*"
        q_sq      = hit.query_seq + "*"

        hs      = '>P1;{0}\nsequence:{0}:{1}:{2}:{3}:{2}:.:.:.:.'
        header1 = hs.format(self.query, hit.query_pos[0], '.', hit.query_pos[-1])
        header2 = hs.format(hit_id, hit.hit_pos[0], hit_chain, hit.hit_pos[-1])

        lines.append(header1)
        lines.extend([q_sq[i:i+60] for i in range(0, len(q_sq), 60)])
        lines.append(header2)
        lines.extend([hit_sq[i:i+60] for i in range(0, len(hit_sq), 60)])

        return "\n".join(lines)

    def str_blast_header(self, as_string = False):
        '''
        Returns a readable account of the blast execution main properties.

        @param:    as_string
        @pdef:     instead of returning each line as a position in an array,
                   it already returns everything as a string
        @pdefault: _False_
        @ptype:    {Boolean}

        @returns: {List} or {String}
        '''
        return self._header.represent(self.query, as_string)

    def str_representation(self, line_split = 160):
        '''
        Creates a representation of the coverage of each hit over the
        query sequence.

        @param:    line_split
        @pdef:     number of characters per line
        @pdefault: 160
        @ptype:    {Integer}
        '''
        names = []
        seqs  = []

        n = line_split

        names.append(self.query)
        line = ''.join(['*'] * self.query_length)
        seqs.append([line[i:i+n] for i in range(0, len(line), n)])
        for hit in self.hits:
            names.append(hit.sequenceID)
            sequence = []
            symbol   = ' '
            ini, end = hit.query_threshold_coord
            for i in range(self.query_length):
                if self._query_index is not None:
                    position = self._query_index[i]
                else:
                    position = i + 1
                if position == ini:
                    symbol = '*'
                sequence.append(symbol)
                if position == end:
                    symbol = ' '
            line = ''.join(sequence)
            seqs.append([line[i:i+n] for i in range(0, len(line), n)])

        data = []
        for j in range(len(seqs[0])):
            for i in range(len(names)):
                data.append('{0:<10}'.format(names[i]) + ' ' + seqs[i][j])
            data.append('')
        return "\n".join(data)

    def print_query_fasta(self):
        '''
        Gives the query in FASTA format.

        @return: {String}
        '''
        return self._query.format('FASTA')

    def print_compacted_blast(self, out_file = None):
        '''
        Print the compacted format of the blast hit.

        @param:    out_file
        @pdef:     file to print the blast data into.
        @pdefault: _None_
        @ptype:    {String}
        '''
        if out_file is not None:
            output = File(out_file, 'w')
            output.write("%s\n" % self.str_compacted_blast())
            output.close()
        else:
            print self.str_compacted_blast()

    def print_representation(self, line_split = 160, out_file = None):
        '''
        Print the alignment representation of the blast hit.

        @param:    line_split
        @pdef:     number of characters per line
        @pdefault: 160
        @ptype:    {Integer}

        @param:    out_file
        @pdef:     file to print the blast data into.
        @pdefault: _None_
        @ptype:    {String}
        '''
        if out_file is not None:
            output = File(out_file, 'w')
            output.write("%s\n" % self.str_representation(line_split))
            output.close()
        else:
            print self.str_representation(line_split)

    @staticmethod
    def read_compacted_blast(compacted_blast_file):
        '''
        Read data from a printed compacted blast into {BlastResult}.
        Not all options will be available in that new object.

        @param:    compacted_blast_file
        @pdef:     file of the compacted blast print
        @ptype:    {String}

        @return: {BlastResult}
        '''
        from BlastHit import BlastHit
        query_name, query_sequence     = None, None
        version, matrix, database      = None, None, None
        gap_open, gap_extend, self_hit = None, None, None

        br = None

        cbf = File(compacted_blast_file)
        for line in cbf.read():
            if line.startswith('#'):
                if line.startswith('#Query:'):
                    query_name = line.strip().split()[-1]
                if line.startswith('#Query Sequence:'):
                    query_sequence = line.strip().split()[-1]
                if line.startswith('#Blast Version:'):
                    version = line.strip().split()[-1]
                if line.startswith('#Search on matrix:'):
                    matrix = line.strip().split()[-1]
                if line.startswith('#Gap open penalty:'):
                    gap_open = line.strip().split()[-1]
                if line.startswith('#Gap extension penalty:'):
                    gap_extend = line.strip().split()[-1]
                if line.startswith('#Database searched:'):
                    database = line.strip().split()[-1]
                if line.startswith('#Self Hit is omitted:'):
                    self_hit = line.strip().split()[-1]
            else:
                if br is None:
                    if version is None:
                        bh = None
                    else:
                        bh = BlastHeader(version, matrix, gap_open, gap_extend,
                                         database, self_hit)
                    br = BlastResult(query_name, query_sequence, bh)
                d = line.strip().split()
                hit = BlastHit([d[2], d[3]], [d[8], d[9]],
                               [int(x) for x in d[10].split(',')[0].split(':')],
                               1, [d[4], d[5], d[6], d[7]])
                br.add_hit(hit)
        cbf.close()

        return br

    ###################
    # PRIVATE METHODS #
    ###################
    def _get_nonoverlaping_hits(self, hitlist, overlap):
        '''
        Filters a list of selected hits according to their overlap, prioritizing
        the order of the list.

        @param:    hitlist
        @pdef:     list of selected hits
        @ptype:    {List} of {BlastHit}

        @param:    overlap
        @pdef:     maximum amount of sequence overlap (over the query) allowed
                   between results. from 0 (no overlap) to 1 (total overlap).
        @ptype:    {Float}

        @return: {List} of {BlastHit}
        '''
        final_hits = []
        final_hits.append(hitlist[0])
        for x in xrange(1, len(hitlist)):
            accepted = True
            for y in xrange(0, len(final_hits)):
                if hitlist[x].overlap(final_hits[y]) > overlap:
                    accepted = False
                    break
            if accepted:
                final_hits.append(hitlist[x])

        self._filter_hits = final_hits
        return self._filter_hits

    #################
    # CLASS METHODS #
    #################
    def __len__(self):
        c = 0
        for hit in self.raw_hits:
            if hit.iteration == self._filter.iteration:
                c += 1
            if hit.iteration > self._filter.iteration:
                break
        return c

    def __str__(self):
        lines = self.str_blast_header()
        for hit in self._hits:
            lines.append("\t".join((self.query, str(self.query_length), str(hit))))
        return "\n".join(lines)


class BlastHeader(object):
    '''
    Stores the main data that identifies the blast execution. For comparisons and
    as a record.

    '''
    def __init__(self, version, matrix, gap_open, gap_extend, database, self_hit):
        '''
        @param:    version
        @pdef:     version of the BLAST executable
        @ptype:    {String}

        @param:    matrix
        @pdef:     name of the similarity matrix
        @ptype:    {String}

        @param:    gap_open
        @pdef:     cost of opening a gap
        @ptype:    {Integer}

        @param:    gap_extend
        @pdef:     cost of extending a gap
        @ptype:    {Integer}

        @param:    database
        @pdef:     file name of the database
        @ptype:    {String}

        @param:    self_hit
        @pdef:     match with itself has not been captured (if any)
        @ptype:    {Boolean} or {String} ('TRUE', 'FALSE')
        '''
        self.version    = version
        self.matrix     = matrix
        self.gap_open   = int(gap_open)
        self.gap_extend = int(gap_extend)
        self.database   = database
        if isinstance(self_hit, basestring):
            if self_hit.upper() == 'TRUE':
                self.self_hit = True
            else:
                self.self_hit = False
        else:
            self.self_hit     = self_hit

    ###########
    # METHODS #
    ###########
    def represent(self, query_name, as_string = False):
        '''
        Returns a readable account of the blast execution main properties.

        @param:    query_name
        @pdef:     name of the query sequence
        @ptype:    {String}

        @param:    as_string
        @pdef:     instead of returning each line as a position in an array,
                   it already returns everything as a string
        @pdefault: _False_
        @ptype:    {Boolean}

        @returns: {List} or {String}
        '''
        lines = []

        lines.append("#Query:                 {0}".format(query_name))
        lines.append("#Blast Version:         {0}".format(self.version))
        lines.append("#Search on matrix:      {0}".format(self.matrix))
        lines.append("#Gap open penalty:      {0}".format(self.gap_open))
        lines.append("#Gap extension penalty: {0}".format(self.gap_extend))
        lines.append("#Database searched:     {0}".format(self.database))
        lines.append("#Self Hit is omitted:   {0}".format(str(not self.self_hit).upper()))

        return '\n'.join(lines) if as_string else lines

    def same_version(self, other):
        '''
        Compares the version with another header.

        @param:    other
        @pdef:     header against which compare
        @ptype:    {BlastHeader}

        @return: {Boolean}
        '''
        return self.version == other.version

    def same_matrix(self, other):
        '''
        Compares the matrix with another header.

        @param:    other
        @pdef:     header against which compare
        @ptype:    {BlastHeader}

        @return: {Boolean}
        '''
        return self.matrix == other.matrix

    def same_gap_open(self, other):
        '''
        Compares the gap_open with another header.

        @param:    other
        @pdef:     header against which compare
        @ptype:    {BlastHeader}

        @return: {Boolean}
        '''
        return self.gap_open == other.gap_open

    def same_gap_extend(self, other):
        '''
        Compares the gap_extend with another header.

        @param:    other
        @pdef:     header against which compare
        @ptype:    {BlastHeader}

        @return: {Boolean}
        '''
        return self.gap_extend == other.gap_extend

    def same_gap_properties(self, other):
        '''
        Compares all gap properties with another header.

        @param:    other
        @pdef:     header against which compare
        @ptype:    {BlastHeader}

        @return: {Boolean}
        '''
        return self.same_gap_open(other) and self.same_gap_extend(other)

    def same_database(self, other):
        '''
        Compares the database with another header.

        @param:    other
        @pdef:     header against which compare
        @ptype:    {BlastHeader}

        @return: {Boolean}
        '''
        return self.database == other.database

    #################
    # CLASS METHODS #
    #################
    def __eq__(self, other):
        a = self.same_version(other) and self.same_matrix(other)
        b = self.same_gap_properties(other) and self.same_database(other)
        return a and b

    def __ne__(self, other):
        return not self.__eq__(other)


class HitFilter(object):
    '''
    Stores different parameters that can be used to filter among all the {BlastHit}
    contained inside the {BlastResult}.
    It is initialized with all the parameters set to values that will not produce
    any kind of filter.
    '''
    _EVALUE     = None
    _ROST_CURVE = None
    _ROST_PARAM = None
    _OVERLAP    = 1

    def __init__(self):
        self._iteration, self._evalue          = None, None
        self._rost_curve, self._rost_parameter = None, None
        self._overlap                          = None

        self._lastiteration = 1

        self.set_default()

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def iteration(self):
        return self._iteration

    @iteration.setter
    def iteration(self, value):
        try:
            iteration = int(value)
        except:
            iteration = self._lastiteration

        if iteration > self._lastiteration:
            iteration = self._lastiteration
        elif iteration == 0:
            iteration = 1

        self._iteration = iteration

    @property
    def lastiteration(self):
        '''
        Last iteration of the blast

        @return: {Integer}
        '''
        return int(self._lastiteration)

    @lastiteration.setter
    def lastiteration(self, value):
        self._lastiteration = int(value)
        self._iteration     = int(value)

    @property
    def evalue(self):
        return self._evalue

    @evalue.setter
    def evalue(self, value):
        self._evalue = float(value) if value is not None else None

    @property
    def tz_type(self):
        return self._rost_curve

    @tz_type.setter
    def tz_type(self, value):
        value = value.upper() if value is not None else None
        self._rost_curve = value if Rost.available_call(value) else None

    @property
    def tz_parameter(self):
        return self._rost_parameter

    @tz_parameter.setter
    def tz_parameter(self, value):
        self._rost_parameter = int(value) if value is not None else None

    @property
    def overlap(self):
        return self._overlap

    @overlap.setter
    def overlap(self, value):
        value = value if value >= 0 else 0
        self._overlap = value if value <= 1 else 1

    ###########
    # METHODS #
    ###########
    def set_default(self):
        '''
        Returns all the parameters to default, non-filtering mode.
        '''
        self._iteration      = self._lastiteration
        self._evalue         = HitFilter._EVALUE
        self._rost_curve     = HitFilter._ROST_CURVE
        self._rost_parameter = HitFilter._ROST_PARAM
        self._overlap        = HitFilter._OVERLAP
