'''
@file: SeqAli.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2013

@ [oliva's lab](http://sbi.imim.es)

@class: SeqAli
@class: Rost
@class: SeqAliError
'''
import re
import decimal
import math

from SBI.sequence import Sequence
from SBI.beans    import IndexedNum


class SeqAli(object):
    '''
    Sequence alignment container.
    Capacity for multiple sequence alignment.
    Iterable. Iteration refers to the *alignment position*, which starts in *1*
    and counts also the gaps.
    A slice from a {SeqAli} also returns a {SeqAli} of the requested section.
    Capacity for *sequence positions* defined as {Integer} or {IndexedNum},
    which is useful for PDB sequence alignment.

    __DEVELOPER ADVISE:__ there is an intended difference between *alignment
    position* (that which refers to global alignment) and *sequence position*,
    which refers to the position of a given sequence according to where it
    begins. As a rule, everything works internally through *alignment position*.
    A request for a *sequence position* is internally transformed.

    '''
    _PATTERN_FORMATS = frozenset(['BLAST'])

    def __init__(self, sequences, sequence_inits,
                 identities = None, positives = None, gaps = None):
        '''
        @param:    sequences
        @pdef:     sequences to add to the alignment.
        @ptype:    {List} of {String} or {Sequence}

        @param:    sequence_inits
        @pdef:     initial number of each sequence of the alignment. (not all
                   sequences start alignment at 1)
        @ptype:    {List} of {Integer}

        @param:    identities
        @pdef:     number of identities in the alignment
        @pdefault: _None_
        @ptype:    {Integer}

        @param:    positives
        @pdef:     number of positives in the alignment
        @pdefault: _None_
        @ptype:    {Integer}

        @param:    gaps
        @pdef:     number of gaps in the alignment
        @pdefault: _None_
        @ptype:    {Integer}

        @raises: {AttributeError} if sequence or sequence_inits is not a {List}
                 or if they are not {List}s of the same length.
                 Also, if the sequence list contains something different than
                 {String} or {Sequence}
        @raises: {SeqAliError} if the fragmentation of the alignment encounters
                 some problem.
        '''
        if not isinstance(sequences, list):
            raise AttributeError('Sequences must be added in a list\n')
        if not isinstance(sequence_inits, list):
            raise AttributeError('Sequence inits must be added in a list\n')
        if len(sequences) != len(sequence_inits):
            raise AttributeError('One init is required for each sequence\n')

        self._error = SeqAliError()
        self._seq   = []

        for aliseq in sequences:
            if isinstance(aliseq, Sequence):
                self._seq.append(aliseq)
            elif isinstance(aliseq, basestring):
                self._seq.append(Sequence(sequence = aliseq))
            else:
                raise AttributeError('sequences must be specified as strings or Sequence objects.')

        self._num_seq = len(self._seq)
        self._segment = []

        self._idx = sequence_inits

        self._seq2ali  = self._built_seq2ali(sequence_inits)
        self._staticSA = self._seq2ali

        self._segmentation_ok = self._search_segments()

        if not self._segmentation_ok:
            raise self._error.wrong_segmentation()

        self._identities = identities
        self._positives  = positives
        self._gaps       = gaps
        self._alipatt    = None
        self._aliptmeth  = None

        self._aligned_aa = self._get_aligned_aa()

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def sequences(self):
        '''
        Sequences in the alignment.

        @return: {List} of {Sequence}
        '''
        return self._seq

    @property
    def identities(self):
        '''
        Identities in the alignment.

        @raises: {TypeError} if identities is _None_
        @return: {Integer}
        '''
        return int(self._identities)

    @property
    def precentage_identity(self):
        '''
        Percentage of identity.

        @raises: {TypeError} if identities is _None_
        @return: {Float}
        '''
        return 100 * float(self.identities) / self.aligned_aminoacids

    @property
    def positives(self):
        '''
        Positives in the alignment.

        @raises: {TypeError} if positives is _None_
        @return: {Integer}
        '''
        return int(self._positives)

    @property
    def precentage_positives(self):
        '''
        Percentage of positives (similarity).

        @raises: {TypeError} if positives is _None_
        @return: {Float}
        '''
        return 100 * float(self.positives) / self.aligned_aminoacids

    @property
    def gaps(self):
        '''
        Gaps in the alignment.

        @raises: {TypeError} if gaps is _None_
        @return: {Integer}
        '''
        return int(self._gaps)

    @property
    def alignment_pattern(self):
        '''
        Pattern of alignment. Summarizes the score of the alignment.

        @return: {String}
        '''
        return self._alipatt

    @property
    def aligned_aminoacids(self):
        '''
        Number of aligned amino acids.

        @return: {Integer}
        '''
        return self._aligned_aa

    @property
    def number_of_sequences(self):
        '''
        Number of sequences in the alignment.

        @return: {Integer}
        '''
        return self._num_seq

    ############
    # BOOLEANS #
    ############
    @property
    def is_multiple(self):
        '''
        Queries if we are dealing with a multiple alignment.

        @return: {Boolean}
        '''
        return self._num_seq > 2

    @property
    def are_segments_ok(self):
        '''
        Checks that the definition of the segments in the alignment is correct.

        @return: {Boolean}
        '''
        return self._segmentation_ok

    ###########
    # METHODS #
    ###########
    def add_complex_index(self, refseq, index, pregaps = 0):
        '''
        Manages the addition of a complex index to one of the sequences in the
        alignment.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    index
        @pdef:     a complex index
        @ptype:    {List}

        @param:    pregaps
        @pdef:     a measure of the gap previous to the position
        @ptype:    {Integer}

        @raises: {AttributeError} if the index is not a list or if it
                 is a list smaller than the number of aligned residues in the
                 sequence.
        '''
        if not isinstance(index, (list, tuple)):
            raise AttributeError('A complex index is a list of IndexedNum or integers')
        if len(index) < len(re.sub('-', '', self._seq[refseq].sequence)):
            raise AttributeError('The given index is not long enough to cover the assigned sequence')

        data        = self._process_index(index)
        part_to_add = data[self._idx[refseq] - 1:]
        self._update_index(refseq, part_to_add)
        self._search_segments()

    def overlap(self, qrefseq, mappedrefseq, otherAli):
        '''
        Measures the overlap between the qrefseq of two alignments when both sequences
        are aligned over mappedrefseq.
        Overlap is measured over the shortest sequence and ranges between 0 and 1.

        mappedrefseq    ------------------------|-------|----------------------
        self.qrefseq     -----------------------|-------|
        otherAli.qrefseq                        |-------|--------------------
                                                 *******

        @param:    qrefseq
        @pdef:     reference sequence (remove lateral gaps from all the alignment).
                   those are the sequence which overlap we evaluate when mapped over
                   the mappedrefseq.
        @ptype:    {Integer}

        @param:    mappedrefseq
        @pdef:     sequence to evaluate overlap upon.
        @ptype:    {Integer}

        @param:    otherAli
        @pdef:     other alignment to which we are comparing. Is mandatory that
                   the alignment has at least as many sequences as the requested.
        @ptype:    {SeqAli}

        @return: {Float}
        '''
        length = min(len(self.sequences[mappedrefseq].ungapped()),
                     len(otherAli.sequences[mappedrefseq].ungapped()))
        mf = self._first_segment_value(qrefseq)
        ml = self._last_segment_value(qrefseq)
        me = self.get_section_from_sequence_position(qrefseq, mf, ml)

        hf = otherAli._first_segment_value(qrefseq)
        hl = otherAli._last_segment_value(qrefseq)
        he = otherAli.get_section_from_sequence_position(qrefseq, hf, hl)

        mf = me._first_segment_value(mappedrefseq)
        ml = me._last_segment_value(mappedrefseq)
        hf = he._first_segment_value(mappedrefseq)
        hl = he._last_segment_value(mappedrefseq)

        imax = self._max_position(mappedrefseq, mf, hl)
        emin = self._min_position(mappedrefseq, ml, hf)
        if (imax == mf and not mf == hl) or \
           (emin == ml and not ml == hf):
            return 0

        fmax = self._max_position(mappedrefseq, mf, hf)
        lmin = self._min_position(mappedrefseq, ml, hl)

        me = self.get_section_from_sequence_position(mappedrefseq, fmax, lmin)
        he = otherAli.get_section_from_sequence_position(mappedrefseq, fmax, lmin)

        l = self.get_sequence_section_from_sequence_position(mappedrefseq,
                                                             fmax, lmin,
                                                             gapped = False)

        return float(len(l)) / float(length)

    def get_sequence_section_from_sequence_position(self, refseq, ini, end, gapped = True):
        '''
        Retrieve a sequence of interest between two points of that sequence.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    ini
        @pdef:     initial *sequence position*
        @ptype:    {Integer}

        @param:    end
        @pdef:     end *sequence position*
        @ptype:    {Integer}

        @param:    gapped
        @pdef:     decide whether to return or not the gaps in the sequence
        @pdefault: _True_
        @ptype:    {Boolean}

        @return: {Sequence}
        '''
        qseq = self.get_section_from_sequence_position(refseq, ini, end)._seq[refseq]

        if gapped:
            return qseq
        else:
            qseq.do_ungap()
            return qseq

    def get_segment_section_from_sequence_position(self, refseq, ini, end):
        '''
        Retrieves the section definition for a defined segment of the alignemnt.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    ini
        @pdef:     initial *sequence position*
        @ptype:    {Integer}

        @param:    end
        @pdef:     end *sequence position*
        @ptype:    {Integer}

        @return: {List}
        '''
        return self.get_section_from_sequence_position(refseq, ini, end)._segment[refseq]

    def increment_sequence_index(self, refseq, new_index):
        '''
        If the start index is corrected _a posteriori_, (for example because the
        sequences are PDB sequences and they did not start in 1 in the first place)
        this correction needs to be extended to the section indexes.

        @param:    refseq
        @pdef:     identifier of the reference sequence
        @ptype:    {Integer}

        @param:    new_index
        @pdef:     new index by which to correct the rest.
        @ptype:    {Integer}
        '''
        self._idx[refseq] += (new_index - 1)
        for x in xrange(len(self._segment[refseq])):
            self._segment[refseq][x] = self._segment[refseq][x] - 1 + self._idx[refseq]

    def get_respective_coordinate(self, refseq, destseq, pos):
        '''
        From the *sequence position* of one sequence, obtain the *sequence position*
        of another.

        @param:    refseq
        @pdef:     identifier of the reference sequence
        @ptype:    {Integer}

        @param:    destseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    pos
        @pdef:     *sequence position* of the reference sequence
        @ptype:    {Integer}

        @return: {Integer}
        '''
        pos = self._alignment_position_from_sequence_position(refseq, pos)
        return self._sequence_position_from_alignment_position(destseq, pos)

    def get_coverage_of_sequence(self, refseq, refseq_full_length = None,
                                 ini = None, end = None):
        '''
        If refseq_full_length is specified, calculates the coverage of the alignment
        over the complete sequence (including the part that does not appear in
        the alignment)
        If ini and/or end are specified, it calculates how much that section of
        sequence is covered by the alignment. In multiple sequence alignments this
        implies that at least is aligned with one other sequence. They supersede
        refseq_full_length.

        @param:    refseq
        @pdef:     identifier of the reference sequence
        @ptype:    {Integer}

        @param:    refseq_full_length
        @pdef:     real length of the sequence, not only of the aligned region.
        @pclash:   ini, end. If specified any of them this is ignored
        @ptype:    {Integer}

        @param:    ini
        @pdef:     initial *sequence position* if we want to calculate only a section
        @ptype:    {Integer}

        @param:    end
        @pdef:     end *sequence position* if we want to calculate only a section
        @ptype:    {Integer}

        @return: {Float}
        '''
        is_section = (ini is not None or end is not None)

        ini = self._first_segment_value(refseq) if ini is None else ini
        end = self._last_segment_value(refseq) if end is None else end

        if ini is None and end is None:  # this sequence is all gaps
            return 0

        section = self.get_section_from_sequence_position(refseq, ini, end)

        if is_section:
            refseq_full_length = section._aligned_sequence_length(refseq)

        coverage = 0
        tk_seqs  = section._tokenize()
        for x in range(len(tk_seqs[refseq])):
            if bool(tk_seqs[refseq][x]) and sum(map(lambda n: int(n[x]), tk_seqs)) > 1:
                coverage += 1

        return float(coverage) / int(refseq_full_length)

    def get_section_from_sequence_position(self, refseq, ini, end):
        '''
        Extracts a section of the alignment according to the given positions
        of the sequence of interest.

        @param:    refseq
        @pdef:     identifier of the reference sequence
        @ptype:    {Integer}

        @param:    ini
        @pdef:     initial *sequence position*
        @ptype:    {Integer}

        @param:    end
        @pdef:     end *sequence position*
        @ptype:    {Integer}

        @return: {SeqAli}
        '''

        align_ini = self._alignment_position_from_sequence_position(refseq, ini)
        align_end = self._alignment_position_from_sequence_position(refseq, end)

        if align_ini > len(self) or align_end > len(self):
            raise IndexError
        if align_ini <= 0 or align_end <= 0:
            raise IndexError
        if align_ini == align_end:
            align_end = None if align_ini == len(self) else (align_end + 1)

        return self.__getitem__(slice(align_ini, align_end, 1))

    def add_alignment_pattern(self, align_pattern, method = 'blast'):
        '''
        Add a summary line that summarizes the alignment, as it exists in the
        blast output.
        If identities, positives and gaps have not been defined before, they
        are filled, as long as the summary format is known and its identification
        provided through the method attribute.

        @param:    align_pattern
        @pdef:     summary line of the alignment.
        @ptype:    {String}


        @param:    method
        @pdef:     format of the summary pattern
        @pdefault: 'blast'
        @ptype:    {String}
        '''
        self._alipatt   = align_pattern
        self._aliptmeth = method

        if self._positives is None and self._identities is None and self._gaps is None:
            if self._aliptmeth == 'blast' and not self.is_multiple:
                self._positives  = 0
                self._identities = 0
                self._gaps       = 0
                for x in range(len(self._alipatt)):
                    position = self._alipatt[x]
                    if position == '+':
                        self._positives  += 1
                    elif position == ' ':
                        if bool(re.search('-', "".join(self[x + 1]))):
                            self._gaps += 1
                    else:
                        self._identities += 1
                        self._positives  += 1
            else:
                return NotImplemented

    def trim_sequence_from(self, source_id, query_id, gapped = False):
        '''
        Extracts the part of the query sequence that aligns with the source
        sequence trimming the unaligned lateral sections.

        @param:    source_id
        @pdef:     identifier of the source sequence
        @ptype:    {Integer}

        @param:    query_id
        @pdef:     identifier of the query sequence that will be returned
        @ptype:    {Integer}

        @param:    gapped
        @pdef:     decide whether to return or not internal gaps in the sequence
        @pdefault: _False_
        @ptype:    {Boolean}

        @return: {String}
        '''
        seq = ''
        letter_found = False
        source = self.sequences[source_id]
        target = self.sequences[query_id]
        gap    = Sequence.GAP_DEFINITION
        for i in range(len(source)):
            if not letter_found:
                if not bool(re.search(gap, source[i])):
                    letter_found = True
                    if not bool(re.search(gap, target[i])) or gapped:
                        seq += target[i]
            else:
                if not bool(re.search(gap, source[i])):
                    if not bool(re.search(gap, target[i])) or gapped:
                        seq += target[i]
                else:
                    if len(re.sub(gap, '', source[i:])) > 0:
                        if not bool(re.search(gap, target[i])) or gapped:
                            seq += target[i]
                    else:
                        break
        return seq

    def evaluate_Rost_twilight_zone(self, equation = 'ID', parameter = None):
        '''
        Rost curves for evaluation homology.
        Use them to decide whether the aligned segments are homologous.

        Rost's Twilight Zone can only be applyed to alignments between
        TWO (2) SEQUENCES.
        This means is not allowed for MULTIPLE SEQUENCE ALIGNMENT.

        @param:    equation
        @pdef:     type of curve to evaluate homology
        @pdefault: 'identity'
        @poptions: 'id', 'identity', 'sim', 'similarity', 'hssp'
        @ptype:    {String}

        @@param:   parameter
        @pdef:     value of the corrector of the curve. The curve gets more
                   restrictive as the value gets bigger.
        @pdefault: 5 for 'identity', 12 for 'similarity' and 8 for 'hssp'
        @ptype:    {Integer}

        @raises: {SeqAliError} if it is multiple alignment.
        @raises: {AttributeError} if requested equation is not available.

        @return: {Boolean}
        '''
        if self.is_multiple:
            raise self._error.two_sequence_method()

        if equation.upper() in ['ID', 'IDENTITY'] :
            r = Rost.identity_curve(self.aligned_aminoacids, parameter)
            isHomolog = self.precentage_identity >= r
        elif equation.upper() in ['SIM', 'SIMILARITY']:
            r = Rost.similarity_curve(self.aligned_aminoacids, parameter)
            isHomolog = self.precentage_positives  >= r
        elif equation.upper() == 'HSSP':
            r = Rost.hssp_curve(self.aligned_aminoacids, parameter)
            isHomolog = self.precentage_identity >= r
        else:
            raise AttributeError('Selected equation not available.')

        return isHomolog

    def format_positions(self, human_readable = False):
        '''
        Creates a visual representation of the segments of the alignment.

        @param:    human_readable
        @pdef:     select if the visual representation is expected to be read
                   by humans.
        @pdefault: _False_
        @ptype:    {Boolean}

        @return: {String}
        '''
        if not human_readable:
            return self._compressed_segment_format()
        else:
            return self._human_segment_format()

    ###################
    # PRIVATE METHODS #
    ###################
    def _max_position(self, refseq, pos1, pos2):
        '''
        Gives whichever position is the biggest according to the sequence.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    pos1
        @pdef:     position 1 to compare
        @ptype:    {IndexedNum}

        @param:    pos2
        @pdef:     position 2 to compare
        @ptype:    {IndexedNum}

        @return: {IndexedNum}
        '''
        p1, p2 = None, None
        try:
            p1 = self._alignment_position_from_sequence_position(refseq, pos1, static = True)
        except IndexError:
            pass
        try:
            p2 = self._alignment_position_from_sequence_position(refseq, pos2, static = True)
        except IndexError:
            pass
        
        if p1 is None:
            return pos2
        if p2 is None:
            return pos1

        if p1 >= p2:
            return pos1
        else:
            return pos2

    def _min_position(self, refseq, pos1, pos2):
        '''
        Gives whichever position is the smallest according to the sequence.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    pos1
        @pdef:     position 1 to compare
        @ptype:    {IndexedNum}

        @param:    pos2
        @pdef:     position 2 to compare
        @ptype:    {IndexedNum}

        @return: {IndexedNum}
        '''
        maxim = self._max_position(refseq, pos1, pos2)
        if maxim == pos1:
            return pos2
        else:
            return pos1

    def _sequence_position_from_alignment_position(self, refseq, pos):
        '''
        Retrieve the position of a given sequence with respect to the *alignment
        position*. If the alignment position corresponds with a gap, it moves
        up to the next position.
        Remember that *alignment position* starts with 1 not 0.
        Inverse of _alignment_position_from_sequence_position

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    pos
        @pdef:     *alignment position* of interest
        @ptype:    {Integer}

        @return: {Integer}
        '''
        c = 0
        while self._seq2ali[refseq][pos - 1 + c] is None:
            c += 1
        return self._seq2ali[refseq][pos - 1 + c]

    def _alignment_position_from_sequence_position(self, refseq, pos, static = False):
        '''
        From a *sequence position* obtain the corresponding *alignment position*.
        Inverse of _sequence_position_from_alignment_position

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    pos
        @pdef:     *sequence position* of interest
        @ptype:    {Integer}

        @param:    static
        @pdef:     When _True_ calls the staticSA that contains the seq2ali before splits
        @pdefault: _False_
        @ptype:    {Boolean}

        @raise: {IndexError} if position is outside of aligned sequence
        @return: {Integer}
        '''
        if not static:
            data = self._seq2ali[refseq]
        else:
            data = self._staticSA[refseq]
        pos = IndexedNum(pos)
        try:
            return data.index(pos) + 1
        except ValueError:
            raise IndexError

    def _sequence_position_id(self, refseq, pos):
        '''
        retrieves the direct position of a requested sequence according to its
        defined init. Considers as if there were no gaps.

        Almost exclusively used for _search_segments()

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    pos
        @pdef:     *sequence position* of interest
        @ptype:    {Integer}

        @return: {Integer}
        '''
        return self._seq2ali[refseq][pos]

    def _first_segment_value(self, refseq):
        '''
        Report the starting values of a given sequence.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @return: {Integer}
        '''
        for i in range(len(self)):
            if self._seq2ali[refseq][i] is not None:
                return self._seq2ali[refseq][i]

    def _last_segment_value(self, refseq):
        '''
        Report the finishing values of a given sequence.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @return: {Integer}
        '''
        for i in reversed(range(len(self))):
            if self._seq2ali[refseq][i] is not None:
                return self._seq2ali[refseq][i]

    def _compressed_segment_format(self):
        '''
        Creates a compressed format to represent the segments of the alignment.
        The format is such as:
        qpos[0]:hpos[0],qpos[1]:hpos[1];qpos[2]:hpos[2],qpos[3]:hpos[3]; [...]
        ^       ^       ^       ^       ^       ^       ^       ^
        |       |       |       |       |       |       |       |Hit segment 2 end
        |       |       |       |       |       |       |Query segment 2 end
        |       |       |       |       |       |Hit segment 2 start
        |       |       |       |       |Query segment 2 start
        |       |       |       |Hit segment 1 end
        |       |       |Query segment 1 end
        |       |Hit segment 1 start
        |Query segment 1 start

        @return: {String}
        '''
        outdata = []
        for x in range(len(self._segment[0])):
            text = []
            for i in range(len(self._segment)):
                text.append("{0}".format(self._segment[i][x]))
            outdata.append(":".join(text))
            if x % 2:
                if x != len(self._segment[0]) - 1:
                    outdata.append(";")
            else:
                outdata.append(",")
        return "".join(outdata)

    def _human_segment_format(self):
        '''
        Creates a easy human readable format to represent the segments of the
        alignment.

        @return: {String}
        '''
        outdata = []
        for i in range(len(self._segment)):
            text = []
            for x in range(len(self._segment[i])):
                text.append("{0:>4}".format(self._segment[i][x]))
                if not x % 2:
                    text.append(" :")
                elif x != len(self._segment[0]) - 1:
                    text.append(" |")
            text = "".join(text)
            outdata.append(text)
        return "\n".join(outdata)

    def _tokenize(self):
        '''
        creates a binary representation of the alignment in an array.
        0 represents gaps while 1 represent aligned amino acids.

        @return: {List} of {String}
        '''
        token_seqs = []
        for seq in self._seq:
            token_seqs.append(seq.tokenize('binary'))
        return token_seqs

    def _get_aligned_aa(self):
        '''
        Determines the total number of aligned amino acids.
        At least two sequences need to be aligned to consider the amino acids
        aligned.

        @return: {Integer}
        '''
        this_len = 0
        tk_seqs = self._tokenize()

        for x in range(len(self._seq[0])):
            profile_str = "".join(map(lambda n: n[x], tk_seqs))
            profile_vle = sum([int(i) for i in profile_str])
            if profile_vle > 1:
                this_len += 1
        return this_len

    def _aligned_sequence_length(self, refseq):
        '''
        Calculates the lenght of a specified sequence of the alignment,
        that means without gaps.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @return: {Integer}
        '''
        return len(re.sub('-', '', self.sequences[refseq].sequence))

    def _check_segments(self):
        '''
        Checks that the number of segment separators make sense and, thus,
        no problem has been found when analyzing the segments.

        @return: {Boolean}
        '''
        static_length = len(self._segment[0])
        if static_length % 2 != 0:
            return False
        for i in range(len(self._segment)):
            if len(self._segment[i]) != static_length:
                return False
        return True

    def _search_segments(self):
        '''
        Revises the alignment and separates the different aligned blocks.

        @return: {Boolean}
        '''
        previous_profile = ''
        self._segment = []
        for i in range(self._num_seq):
            previous_profile += '0'
            self._segment.append([])
        tk_seqs = self._tokenize()
        for x in range(len(self)):
            profile_str = "".join(map(lambda n: n[x], tk_seqs))

            for i in range(self._num_seq):
                p0 = self._sequence_position_id(i, x)
                p1 = self._sequence_position_id(i, (x - 1))
                b0 = bool(int(profile_str[i]))
                b1 = bool(int(previous_profile[i]))
                if x == 0:
                    self._segment[i].append(p0 if b0 else '-')
                else:
                    if profile_str != previous_profile:
                        self._segment[i].append(p1 if b1 else '-')
                        self._segment[i].append(p0 if b0 else '-')
                    if x == (len(self) - 1):
                        self._segment[i].append(p0 if b0 else '-')

            previous_profile = profile_str

        x = 0
        while (x < len(self._segment[0])):
            counter = 0
            for y in range(self._num_seq):
                counter = counter + (1 if self._segment[y][x] != '-' else 0)
            if counter == 1:
                for y in range(self._num_seq):
                    self._segment[y].pop(x)
                    self._segment[y].pop(x)
                x = x - 2
            x += 2
        if len(self) == 1:
            for y in range(self._num_seq):
                self._segment[y].append(self._segment[y][0])
        return self._check_segments()

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
        newali = self.__class__(sequences, seqInits)
        for refident in complex_indexes:
            newali._update_index(refident, complex_indexes[refident])
            newali._search_segments()
        if not self._alipatt is None:
            newali.add_alignment_pattern(self._alipatt[index[0]-1:index[1]-1:index[2]], self._aliptmeth)
        return newali


    def _process_index(self, index):
        '''
        Enumerates the 'X' of the index

        @param:    index
        @pdef:     transformation and adaptation of the slice object
        @ptype:    {List}

        @return: {List}
        '''
        a = []
        x = 1
        for i in range(len(index)):
            add = index[i] if index[i] != 'X' else 'X' + str(x)
            a.append(IndexedNum(add))
            x   = x + 1 if index[i] == 'X' else x
        return a

    def _update_index(self, refseq, index):
        '''
        Updates the seq2ali list that allows moving from *sequence* to *alignment*
        position and back.

        @param:    refseq
        @pdef:     identifier of the sequence of interest
        @ptype:    {Integer}

        @param:    index
        @pdef:     transformation and adaptation of the slice object
        @ptype:    {List}

        '''
        a = list(self._seq2ali[refseq])
        c = 0
        x = 1
        for i in range(len(a)):
            if a[i] is not None:
                add  = index[c] if index[c] != 'X' else 'X' + str(x)
                a[i] = IndexedNum(add)
                x    = x + 1 if index[c] == 'X' else x
                c   += 1
        self._seq2ali[refseq] = tuple(a)

        self._idx[refseq] = self._first_segment_value(refseq)

    def _built_seq2ali(self, sequence_inits):
        '''
        builts the seq2ali vector required to move from *sequence* to *alignemnt*
        positions.

        @param:    sequence_inits
        @pdef:     initial number of each sequence of the alignment. (not all
                   sequences start alignment at 1)
        @ptype:    {List} of {Integer}

        @return: {List} of {Tuples}
        '''
        r = []
        for i in range(self.number_of_sequences):
            c, a = 0, []
            for j in range(len(self._seq[i])):
                s = self._seq[i][j]
                if s == '-':
                    a.append(None)
                else:
                    a.append(sequence_inits[i] + c)
                    c += 1
            r.append(tuple(a))
        return r

    #################
    # CLASS METHODS #
    #################
    def __len__(self):
        return len(self._seq[0])

    def __getitem__(self, key):
        # We correct the count as one would try to ask for the first position
        # as 1, not as 0...
        # Thus, we must consider that alignments start always in 1 (as global,
        # not for each particular sequence)
        try:
            int(key)
        except:
            if not isinstance(key, slice):
                raise TypeError()
        else:
            if int(key) > len(self):
                raise IndexError
            return map(lambda n: n[int(key) - 1], self._seq)

        index = key.indices(len(self) + 1)

        if key.stop is not None:
            #Remember! we go according to align position, not array!!
            index = (index[0], index[1] + 1, index[2])
        else:
            index = (index[0], index[1], index[2])
        sequences = ['' for x in range(self.number_of_sequences)]
        for i in range(*index):
            data = self[i]
            for x in range(len(data)):
                sequences[x] += data[x]
            if i == len(self):
                break
        ali = self._build_slice_alignment(sequences, index)
        ali._staticSA = self._staticSA
        return ali

    def __iter__(self):
        for x in range(1, len(self) + 1):
            yield self[x]

    def __repr__(self):
        data = []
        for x in range(self._num_seq):
            data.append('{0}\t{1}\t{2}'.format(self._first_segment_value(
                x), self._seq[x].sequence, self._last_segment_value(x)))
        if not self._alipatt is None:
            data.append('\t{0}'.format(self._alipatt))
        if not self._identities is None and not self._positives is None and not self._gaps is None:
            data.append(
                'I: {0._identities:<4} P: {0._positives:<4} G: {0._gaps:<4}'.format(self))
        data.append(self.format_positions(True))
        return "\n".join(data)


class Rost(object):
    '''
    FORMULA DEFINITIONS FOR ROST TWILIGHT ZONE:

    Rost, B. (1999). Twilight zone of protein sequence alignments.
    Protein Engineering, 12(2), 85 94.

    MODIFIED FROM J.Planas-Iglesias IMPLEMENTATION

    '''
    E = decimal.Decimal(repr(math.e))
    P = {'identity':   decimal.Decimal('-0.32'),
         'similarity': decimal.Decimal('-0.335'),
         'hssp':       decimal.Decimal('-0.562')}
    V = {'identity': 480,  'similarity': 420, 'hssp': 290.15}
    D = {'identity': 1000, 'similarity': 2000}
    N = {'identity':   decimal.Decimal(5),
         'similarity': decimal.Decimal(12),
         'hssp':       decimal.Decimal(8)}

    CALLS = frozenset(['ID', 'IDENTITY', 'SIMILARITY', 'SIM', 'HSSP'])

    @staticmethod
    def available_call(call):
        '''
        Checks if the way to refer to Rost curves is correct.

        @param:    call
        @pdef:     identifier of the curve of interest
        @ptype:    {String}

        @return: {Boolean}
        '''
        return call in Rost.CALLS

    @staticmethod
    def hssp_curve(aligned_aminoacids, curve_parameter = None):
        '''
        HSSP identity curve. As defined in Rost's article.
        Fulfilling the curve implies that the percentage of identity in the
        alignment is greater or equal than the returned value of the curve for
        a given number of aligned residues.
        According to the article, Rost's curves are more accurate, at the
        expense of coverage.
        Original description of the curve in:

        Sander, C., & Schneider, R. (1991). Database of homology-derived protein
        structures and the structural meaning of sequence alignment.
        Proteins, 9(1), 56-68.

        @param:    aligned_aminoacids
        @pdef:     total number of aligned amino acids.
        @ptype:    {Integer}

        @param:    curve_parameter
        @pdef:     (N). This parameters alters the high of the curve.
                   The bigger the number, more restrictive the curve is.
        @pdefault: 8. As defined in Rost's article.
        @ptype:    {Integer}

        @return: {Decimal}
        '''
        if curve_parameter is None:
            curve_parameter = Rost.N['hssp']

        if aligned_aminoacids >= 80:
            p = 25
        else:
            p = decimal.Decimal(repr(Rost.V['hssp']))
            p = p * pow(aligned_aminoacids, Rost.P['hssp'])

        return curve_parameter + p

    @staticmethod
    def identity_curve(aligned_aminoacids, curve_parameter = None):
        '''
        Rost identity curve.
        Fulfilling the curve implies that the percentage of identity in the
        alignment is greater or equal than the returned value of the curve for
        a given number of aligned residues.

        @param:    aligned_aminoacids
        @pdef:     total number of aligned amino acids.
        @ptype:    {Integer}

        @param:    curve_parameter
        @pdef:     (N). This parameters alters the high of the curve.
                   The bigger the number, more restrictive the curve is.
        @pdefault: 5. As defined in Rost's article.
        @ptype:    {Integer}

        @return: {Decimal}
        '''
        return Rost._curve(aligned_aminoacids, curve_parameter, 'identity')

    @staticmethod
    def similarity_curve(aligned_aminoacids, curve_parameter = None):
        '''
        Rost similarity curve.
        Fulfilling the curve implies that the percentage of similarity in the
        alignment is greater or equal than the returned value of the curve for
        a given number of aligned residues.

        @param:    aligned_aminoacids
        @pdef:     total number of aligned amino acids.
        @ptype:    {Integer}

        @param:    curve_parameter
        @pdef:     (N). This parameters alters the high of the curve.
                   The bigger the number, more restrictive the curve is.
        @pdefault: 12. As defined in Rost's article.
        @ptype:    {Integer}

        @return: {Decimal}
        '''
        return Rost._curve(aligned_aminoacids, curve_parameter, 'similarity')

    @staticmethod
    def _curve(aligned_aminoacids, curve_parameter, curve):
        '''
        Executes the requested curve.

        @param:    aligned_aminoacids
        @pdef:     total number of aligned amino acids.
        @ptype:    {Integer}

        @param:    curve_parameter
        @pdef:     (N). This parameters alters the high of the curve.
                   The bigger the number, more restrictive the curve is.
        @ptype:    {Integer}

        @param:    curve
        @pdef:     Which curve is requested.
        @poptions: 'similarity' or 'identity'
        @ptype:    {String}

        formula : N + (V * pow(L, P * (1 + pow(E, (-L / D)))))

        @return: {Decimal}
        '''
        if curve_parameter is None:
            curve_parameter = Rost.N[curve]
        else:
            curve_parameter = decimal.Decimal(curve_parameter)
        m = decimal.Decimal(repr(float(-aligned_aminoacids) / Rost.D[curve]))
        e = pow(Rost.E, m)
        k = Rost.P[curve] * (1 + e)
        curve_shape = Rost.V[curve] * pow(aligned_aminoacids, k)
        return curve_parameter + curve_shape


class SeqAliError(Exception):
    '''
    Manages different errors that can occur during the parsing of an alignment.

    '''

    _MSG = ''

    def __init__(self):
        pass

    def wrong_segmentation(self):
        self._MSG  = 'There has been a problem when parsing the different segments '
        self._MSG += 'of the alignment.'
        return self

    def two_sequence_method(self):
        self._MSG  = 'The required function can only be applied to two sequence'
        self._MSG += ' alignments.'
        return self

    def __str__(self):
        return SeqAliError._MSG
