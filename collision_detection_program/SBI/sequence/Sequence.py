'''
@file: Sequence.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   05/2013

@ [oliva's lab](http://sbi.imim.es)

@class: Sequence
'''
import re
import copy
import hashlib
import difflib
from collections import Counter


class Sequence(object):
    '''
    A simple container for sequences.
    Includes the sequence identifier, the sequence information and the sequence
    itself.
    Default comparison performed between sequences by length.
    Iterable.

    '''
    AVAILABLE_FORMATS          = set(['TAB', 'FASTA', 'PIR'])
    TOKEN_FORMATS              = set(['BINARY'])
    GAP_DEFINITION             = '[-x]'
    NEEDLEMAN_WUNSCH           = set(['NEEDLEMAN-WUNSCH', 'NW'])
    AVAILABLE_ALIGN_ALGORITHMS = set()

    AVAILABLE_ALIGN_ALGORITHMS.update(NEEDLEMAN_WUNSCH)

    def __init__(self, sequence_id = '', sequence = ''):
        '''
        @param:    sequence_id
        @pdef:     Name of the sequence.
        @pdefault: ''
        @ptype:    {String}

        @param:    sequence
        @pdef:     Sequence.
        @pdefault: ''
        @ptype:    {String} or {List}

        @raises: {AttributeError} if sequence is not {String} or {List}.
        '''
        self._seqID = ''
        self._info  = ''
        self.id     = sequence_id

        self._sequence = ''
        self._gapped   = False
        self.sequence  = sequence

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def id(self):
        '''
        Sequence identifier.

        @return: {String}
        '''
        return self._seqID

    @id.setter
    def id(self, value):
        '''
        @param:    id
        @pdef:     Name of the sequence.
        @pdefault: ''
        @ptype:    {String}
        '''
        if bool(re.search(r'\s', value)):
            self._seqID = value.split()[0]
            self._info  = ' '.join(value.split()[1:])
        else:
            self._seqID = value
            self._info  = ''

    @property
    def sequence(self):
        '''
        Sequence.

        @return: {String}
        '''
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        '''
        @param:    sequence
        @pdef:     Sequence.
        @pdefault: ''
        @ptype:    {String} or {List}
        '''
        if isinstance(value, str):
            self._sequence = value
        elif isinstance(value, list):
            self._sequence = ''.join(value)
        else:
            raise AttributeError('sequence attribute must be string or list.')

        self._gapped = bool(re.search(self.GAP_DEFINITION, self.sequence))

    @property
    def info(self):
        '''
        Sequence description.

        @return: {String}
        '''
        return self._info

    @info.setter
    def info(self, value):
        '''
        @param:    info
        @pdef:     Description of the sequence.
        @pdefault: ''
        @ptype:    {String}
        '''
        self._info = value

    @property
    def full_identifier(self):
        '''
        Full name of the sequence.

        @return: {String}
        '''
        return ' '.join([self.id, self.info]).strip()

    @property
    def md5(self):
        '''
        Creates an md5 of the sequence as a unique identifier.

        @return: {String}
        '''
        pre = self.sequence[:3]
        md5 = hashlib.sha224(self.sequence).hexdigest()
        pos = self.sequence[-3:]
        return pre + md5 + pos

    ############
    # BOOLEANS #
    ############
    @property
    def is_gapped(self):
        '''
        Determines it the sequence contains gaps.

        @return: {Boolean}
        '''
        return self._gapped

    ###########
    # METHODS #
    ###########
    def chunk(self, size = 60):
        '''
        Split the sequence in segments of up to 'size' elements.

        @param:    size
        @pdef:     Max size of the chunks.
        @ptype:    {Integer}
        '''
        return [self.sequence[i:i+size] for i in range(0, len(self), size)]

    def contains(self, sequence):
        '''
        Checks if a given sequence is contained inside.

        @param:    sequence
        @pdef:     Sequence to compare with.
        @ptype:    {String}, {List} or {Sequence}

        @raises: {AttributeError} if sequence is the wrong type.
        @return: {Boolean}
        '''
        if isinstance(sequence, basestring):
            return bool(re.search(sequence, self.sequence))
        elif isinstance(sequence, list):
            return bool(re.search(''.join(sequence), self.sequence))
        elif isinstance(sequence, Sequence):
            return bool(re.search(sequence.sequence, self.sequence))
        else:
            return AttributeError('sequence must be string, list or Sequence')

    def contained(self, sequence):
        '''
        Checks if the sequence is contained inside a given one.

        @param:    sequence
        @pdef:     Sequence to compare with.
        @ptype:    {String}, {List} or {Sequence}

        @raises: {AttributeError} if sequence is the wrong type.
        @return: {Boolean}
        '''
        if isinstance(sequence, basestring):
            return bool(re.search(self.sequence, sequence))
        elif isinstance(sequence, list):
            return bool(re.search(self.sequence, ''.join(sequence)))
        elif isinstance(sequence, Sequence):
            return bool(re.search(self.sequence, sequence.sequence))
        else:
            return AttributeError('sequence must be string, list or Sequence')

    def get_matching_blocks(self, sequence):
        '''
        Returns each identical block between the sequence and a given one

        @param:    sequence
        @pdef:     Sequence to compare with.
        @ptype:    {String}, {List} or {Sequence}

        @raises: {AttributeError} if sequence is the wrong type.
        @return: {List} of {String}
        '''
        if isinstance(sequence, basestring):
            pass
        elif isinstance(sequence, list):
            sequence = ''.join(sequence)
        elif isinstance(sequence, Sequence):
            sequence = sequence.sequence
        else:
            return AttributeError('sequence must be string, list or Sequence')

        m = difflib.SequenceMatcher(None, self.sequence, sequence)
        b = []
        for match in m.get_matching_blocks():
            i, j, n = match
            b.append(self[i:i+n])
        return b

    def align(self, sequence, algorithm = 'NW', similarity_matrix = 'BLOSUM62',
              gap_init = None, gap_penalty = None):
        '''
        Aligns the sequence with a given one.

        @param:    sequence
        @pdef:     Sequence to compare with.
        @ptype:    {String}, {List} or {Sequence}

        @@param:   algorithm
        @pdef:     Algorithm to use for the alignment
        @pdefault: 'needleman_wunsch'
        @ptype:    {String}

        @param:    similarity_matrix
        @pdef:     Scoring matrix for the alignment.
        @pdefault: 'BLOSUM62'
        @ptype:    {String}

        @param:    gap_init
        @pdef:     Value of initiating a gap
        @pdefault: Gap value of the similarity_matrix
        @ptype:    {Integer}

        @param:    gap_penalty
        @pdef:     Value of extending a gap
        @pdefault: Gap value of the similarity_matrix
        @ptype:    {Integer}

        @raises: {AttributeError} if sequence is the wrong type or if a non-
                 available alignment method is requested.
        @return: {SeqAli}
        '''
        if algorithm.upper() not in Sequence.AVAILABLE_ALIGN_ALGORITHMS:
            raise AttributeError('Alignment algorithm not available')
        if isinstance(sequence, basestring):
            pass
        elif isinstance(sequence, list):
            sequence = ''.join(sequence)
        elif isinstance(sequence, Sequence):
            sequence = sequence.sequence
        else:
            return AttributeError('sequence must be string, list or Sequence')

        if algorithm.upper() in Sequence.NEEDLEMAN_WUNSCH:
            from alignment.Needleman_Wunsch import needleman_wunsch as nw
            return nw(self.sequence, sequence, similarity_matrix,
                      gap_init, gap_penalty)

    def format(self, format = 'TAB'):
        '''
        Formats the sequence.

        @param:    format
        @pdef:     Output format for the sequence
        @pdefault: 'TAB'
        @ptype:    {String}

        @raises: {AttributeError} for unknown formats
        @return: {String}
        '''
        if format.upper() not in Sequence.AVAILABLE_FORMATS:
            raise AttributeError('format option not available')

        if format.upper() == 'TAB':
            return '{0.id}\t{0.sequence}'.format(self)
        elif format.upper() == 'FASTA':
            return '>{0.id}\n{0.sequence}'.format(self)
        elif format.upper() == 'PIR':
            data = ['>P1;{0.id}'.format(self),
                    'sequence:{0.id}:::::::0.00: 0.00'.format(self)]
            data.extend(self.chunk())
            return '\n'.join(data) + '*'

    def tokenize(self, token_coding = 'BINARY'):
        '''
        Re-codifies the sequence according to a given codification.

        @param:    token_coding
        @pdef:     codification to tokenize.
        @pdefault: 'BINARY'
        @ptype:    {String}
        '''
        if token_coding.upper() not in Sequence.TOKEN_FORMATS:
            raise AttributeError('Unrecognized token coding.')
        if token_coding.upper() == 'BINARY':
            tc = ['\D', '1', '-', '0']
            return re.sub(tc[0], tc[1], re.sub(tc[2], tc[3], self._sequence))

    def items_frequency(self):
        '''
        Displays the frequency of the different items in the sequence.

        @return: {Dictionary}
        '''
        c = Counter(re.sub(Sequence.GAP_DEFINITION, '', self._sequence))
        l = float(len(self))
        return dict([(x, y/l) for x, y in c.iteritems()])

    def duplicate(self, new_id = None):
        '''
        Creates a new copy of the sequence.

        @param:    new_id
        @pdef:     a new identifier can be set in the copy
        @ptype:    {String}

        @return: {Sequence}
        '''
        s = copy.deepcopy(self)
        if new_id is not None:
            s.id = new_id
        return s

    def do_ungap(self):
        '''
        Removes the gaps in the sequence.
        '''
        if self._gapped:
            self.sequence = re.sub(self.GAP_DEFINITION, '', self.sequence)

    def ungapped(self):
        '''
        Returns a copy of the sequence without gaps.
        Does not alter the sequence.

        @return: {String}
        '''
        return re.sub(self.GAP_DEFINITION, '', self.sequence)

    def append(self, sequence):
        '''
        Extends the sequence.

        @param:    sequence
        @pdef:     Sequence to compare with.
        @ptype:    {String}, {List} or {Sequence}

        @raises: {AttributeError} if sequence is the wrong type.
        '''
        if isinstance(sequence, basestring):
            pass
        elif isinstance(sequence, list):
            sequence = ''.join(sequence)
        elif isinstance(sequence, Sequence):
            sequence = sequence.sequence
        else:
            return AttributeError('sequence must be string, list or Sequence')

        self.sequence = self.sequence + sequence

    #################
    # CLASS METHODS #
    #################
    def __len__(self):
        return len(self._sequence)

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.sequence == other.sequence
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __lt__(self, other):
        if isinstance(other, Sequence):
            return len(self) < len(other)
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, Sequence):
            return len(self) > len(other)
        return NotImplemented

    def __getitem__(self, key):
        try:
            int(key)
            return self._sequence[int(key)]
        except:
            if not isinstance(key, slice):
                raise TypeError
            else:
                return self._sequence[key]

    def __iter__(self):
        for s in self._sequence:
            yield s

    def __repr__(self):
        return '<{0.__class__.__name__}: [{0.id}, {0.sequence}]>'.format(self)

    def __str__(self):
        return '{0.id}\t{0.sequence}'.format(self)
