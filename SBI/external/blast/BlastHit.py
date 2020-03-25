"""
BlastHit

The BlastHit object stores all the data related to a blast XML output hit.

That means the data inside each <Hit> tag of the blast XML file which corresponds to a single blast match to the query.

WARNING!: It has to be a blastpgp from a single sequence; multi-fasta blasts are not compatible
          with xml.dom.minidom parsing capabilities.

                                   ############################################
                                        jbonet @ boliva's lab        2011
"""
import re
from collections import Counter

from SBI.sequence.alignment import SeqAli


class BlastHit(SeqAli):
    '''
    A sequence alignment that stores the data related to a single query-hit
    relation from a Blast output.
    All the query-hit ({BlastHit}) pairs are grouped into a {BlastResult}.

    Normally, the object will be created directly through the {BlastParser}

    '''

    _VAST = {'query': 0, 0: 0, 'hit': 1, 1: 1}

    def __init__(self, hit, sequences, sequence_inits, iteration, stats):
        '''

        @param:    sequences
        @pdef:     list with the aligned query and hit sequences (in this order)
                   plus the alignment pattern representation derived from blast
                   (in the last position of the list).
        @ptype:    {List} of {String}

        @param:    sequence_inits
        @pdef:     list of the first position number for the query and hit
                   sequences, in this order.
        @ptype:    {List} of {Integer}

        @param:    iteration
        @pdef:     iteration of the blast output to which the {BlastHit} belongs
        @ptype:    {Integer}

        @param:    stats
        @pdef:     list alignment statistics. Namely: identities, positives, gaps
                   and e-value
        @ptype:    {List} as [{Integer}, {Integer}, {Integer}, {Float}]
        '''
        align_pattern = sequences.pop() if len(sequences) > 2 else None
        super(BlastHit, self).__init__(sequences      = sequences,
                                       sequence_inits = sequence_inits,
                                       identities     = stats[0],
                                       positives      = stats[1],
                                       gaps           = stats[2])

        if align_pattern is not None:
            self.add_alignment_pattern(align_pattern = align_pattern, method = 'blast')

        self._sequenceID   = str(hit[0])        # Hit Name
        self._length       = int(hit[1])        # Hit Length
        self._iteration    = int(iteration)     # Iteration to which the hit belongs
        self._e_value      = stats[3]           # Alignment e-value

        self._used         = False

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def query_seq(self):
        '''
        The aligned query sequence

        @return: {String}
        '''
        return str(self._seq[0].sequence)

    @property
    def ungapped_query_seq(self):
        '''
        The aligned query sequence. without gaps

        @return: {String}
        '''
        return str(re.sub('-', '', self.query_seq))

    @property
    def hit_seq(self):
        '''
        The aligned hit sequence

        @return: {String}
        '''
        return str(self._seq[1].sequence)

    @property
    def ungapped_hit_seq(self):
        '''
        The aligned hit sequence. without gaps

        @return: {String}
        '''
        return str(re.sub('-', '', self.hit_seq))

    @property
    def query_pos(self):
        return self._segment[0]

    @property
    def hit_pos(self):
        return self._segment[1]

    '''ATTRIBUTES'''
    @property
    def sequenceID(self):
        return self._sequenceID

    @property
    def length(self):
        return self._length

    @property
    def iteration(self):
        return self._iteration

    @property
    def evalue(self):
        return float(self._e_value) if self._e_value is not None else None

    @property
    def is_used(self):
        return self._used

    @is_used.setter
    def is_used(self, value):
        self._used = value

    '''MIX ATTRIBUTES'''
    @property
    def query_threshold_coord(self):
        return [self._first_segment_value(0), self._last_segment_value(0)]

    @property
    def query_align_segment_length(self):
        return self._aligned_sequence_length(0)

    @property
    def hit_threshold_coord(self):
        return [self._first_segment_value(1), self._last_segment_value(1)]

    @property
    def hit_align_segment_length(self):
        return self._aligned_sequence_length(1)

    '''OVERWRITTE PARENT GETTERS TO ADAPT TO refseq=["query","hit"]'''

    def get_hit_to_query_coordinate(self, coordinate):
        refseq = self._check_valid_alignment_seq_type(1)
        return super(BlastHit, self).get_respective_coordinate(refseq, int(not refseq), coordinate)

    def get_query_to_hit_coordinate(self, coordinate):
        refseq = self._check_valid_alignment_seq_type(0)
        return super(BlastHit, self).get_respective_coordinate(refseq, int(not refseq), coordinate)

    def get_coverage_of_sequence(self, refseq, fulllength = None, ini = None, end = None):
        """
        METHOD: get_coverage_of_sequence(refseq =["query","hit", fulllength = None, ini = None, end = None))

            __DEVEL_INFO__
            => Overwrittes from parent (SeqAli)

            __EVALUATOR__
            =>Returns coverage ratio (0<->1) of the quey/hit over the full sequence (for query fulllength is required)

            =>If a *range* is given, returns __the sequence coverage for the given range__
             (adjusted to actually aligned regions)

            =>[TIP] To obtain __the sequence coverage for the aligned region over itself__:
                ini = self.query_threshold_coord[0] or self.hit_threshold_coord[0]
                AND
                end = self.query_threshold_coord[1] or self.hit_threshold_coord[1]
        """
        refseq = self._check_valid_alignment_seq_type(refseq)
        if bool(refseq):  # hit
            fulllength = self.length
        else:  # query
            if fulllength is None:
                raise AttributeError('To check the query coverage the query length must be given')

        return super(BlastHit, self).get_coverage_of_sequence(BlastHit._VAST[refseq], fulllength, ini, end)

    def get_min_coverage_of_sequence(self, fulllength):
        return min(self.get_coverage_of_sequence('hit'), self.get_coverage_of_sequence('query', fulllength))

    def get_coverage_of_full_sequence_segment(self, refseq):
        refseq = self._check_valid_alignment_seq_type(refseq)

        if bool(refseq):  # hit
            length = len(self.hit_seq.replace('-', ''))
        else:  # query
            length = len(self.query_seq.replace('-', ''))

        return super(BlastHit, self).get_coverage_of_sequence(BlastHit._VAST[refseq], length)

    def get_min_coverage_of_full_sequence_segment(self):
        return min(self.get_coverage_of_full_sequence_segment('hit'), self.get_coverage_of_full_sequence_segment('query'))

    def get_section_from_sequence_position(self, refseq, ini, end):
        refseq = self._check_valid_alignment_seq_type(refseq)
        return super(BlastHit, self).get_section_from_sequence_position(refseq, ini, end)

    def correct_hit_count(self, new_index):
        if isinstance(self._idx[1], list):
            raise AttributeError('Complex index of hit cannot be modified')

        if not isinstance(new_index, list):
            self.increment_sequence_index(1, new_index) ##
        else:
            self.add_complex_index(1, new_index)

    def correct_query_count(self, new_index):
        if isinstance(self._idx[0], list):
            raise AttributeError('Complex index of querycannot be modified')

        if not isinstance(new_index, list):
            self.increment_sequence_index(0, new_index)
        else:
            self.add_complex_index(0, new_index)

    def _check_valid_alignment_seq_type(self, seq_type):
        """
        METHOD: _check_valid_alignment_seq_type()

            __STATUS_CHECKER__
            =>Returns BOOLEAN:
                                True if the seq_type is accepted
                                False otherwise
        """
        if not seq_type in BlastHit._VAST:
            raise ValueError("'ref_seq' must be in %s. Recieved 'ref_seq' is %s" % (repr(BlastHit._VAST), seq_type))
        return BlastHit._VAST[seq_type]

    def overlap(self, blastHit):
        # It returns a range between 0<->1 of how much the actual blastHit and the one passed
        # overlap in the query sequence.
        # The overlap is measured OVER THE SHORTEST lenght
        return super(BlastHit, self).overlap(1, 0, blastHit)

    """
        toString
    """
    def _build_slice_alignment(self, sequences, index):
        '''
        This function helps the slicing of the alignment through slice objects as
        an iterable. It is required to support inheritance of the slicing.
        It covers special requirements of child objects when building themselves.

        @param:    sequences
        @pdef:     array of sequences for the new alignment
        @ptype:    {List}

        @param:    index
        @pdef:     transformation and adaptation of the slice object
        @ptype:    {List}

        @return: {IndexedSeqAli}
        '''
        seqInits  = ['' for x in range(self.number_of_sequences)]
        complex_indexes = {}
        for x in range(self.number_of_sequences):
            seqInits[x] = self._sequence_position_from_alignment_position(x, index[0])
            complex_indexes[x] = list(filter(None, self._seq2ali[x]))
            complex_indexes[x] = complex_indexes[x][complex_indexes[x].index(seqInits[x]):]
        patt   = self._alipatt[index[0] - 1:index[1] - 1:index[2]]
        newali = self.__class__(hit = [self._sequenceID, self._length],
                                sequences = [sequences[0], sequences[1], patt],
                                sequence_inits = [seqInits[0], seqInits[1]],
                                iteration      = self._iteration,
                                stats          = [None, None, None, None])
        for refident in complex_indexes:
            newali._update_index(refident, complex_indexes[refident])
            newali._search_segments()
        return newali

    def __str__(self):
        return "\t".join((self.sequenceID, str(self.length), str(self.identities), str(self.positives), str(self.gaps),
                          "{:.3}".format(self.evalue), self.query_seq, self.hit_seq, self.format_positions()))
