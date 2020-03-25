'''
@file: Fasta.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   05/2013

@ [oliva's lab](http://sbi.imim.es)

@class: Fasta
'''
import os

from SBI.beans    import File
from SBI.sequence import Sequence
from SBI          import SBIglobals as SBIg


class Fasta(object):
    '''
    Manages the most common format for sequences.
    It can directly load to memory all the sequences (fast) or retrieve them
    when needed by reading (memory optimized).
    All the identifiers are loaded directly.

    '''
    def __init__(self, fasta_file, auto_load = 10):
        '''
        @param:    fasta_file
        @pdef:     name of the FASTA file.
        @ptype:    {String} or {File}

        @@param:   auto_load
        @pdef:     maximum number of sequences to autoload.
        @pdefault: 10
        @ptype:    {Integer}
        '''
        if isinstance(fasta_file, basestring):
            self._file = File(file_name = fasta_file, action = 'r')
        elif isinstance(fasta_file, File):
            self._file = File(file_name = fasta_file.full, action = 'r')
        else:
            raise AttributeError('Check the input of the Fasta object')

        self._sequences  = []
        self._sequenceID = {}

        self._total_sequences = 0
        self._loaded          = False
        self._auto_load       = auto_load
        self._check_multifasta()

        self._index_file = None
        self._check_index()

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def file(self):
        '''
        Source file of the sequences.

        @return: {File}
        '''
        return self._file

    @property
    def sequences(self):
        '''
        Returns all the sequences loaded in the object.

        @yields: {Sequence}
        '''
        for s in self._sequences:
            if s is not None:
                yield s

    @property
    def sequence(self):
        '''
        Returns the first sequence. Specially designed for single FASTA files.

        @return: {Sequence}
        '''
        return self._sequences[0]

    @property
    def sequence_identifiers(self):
        '''
        Lists all the identifiers in the FASTA file.

        @return: {List}
        '''
        if self.is_loaded:
            return self._sequenceID.keys()
        else:
            keys = []
            for sequence in self.live_show():
                keys.append(sequence.id)
            return keys

    @property
    def index_file(self):
        '''
        Name of the related index file

        @return: {String}
        '''
        return self._index_file

    @index_file.setter
    def index_file(self, value):
        '''
        @param:    value
        @pdef:     name of the related index file
        @ptype:    {String}
        '''
        if os.path.isfile(value):
            self._index_file = value
        else:
            SBIg.warn(self, '{0} does not exist.'.format(value))

    ############
    # BOOLEANS #
    ############
    @property
    def is_multifasta(self):
        '''
        Checks the type of FASTA we are working with.

        @return: {Boolean}
        '''
        return self._total_sequences > 1

    @property
    def is_loaded(self):
        '''
        Checks if all the sequences have been loaded to memory.

        @return: {Boolean}
        '''
        return self._loaded

    @property
    def has_index(self):
        '''
        Checks if has a index associated.

        @return: {Boolean}
        '''
        return self._index_file is not None

    ###########
    # METHODS #
    ###########
    def split(self, sections):
        '''
        @param:    sections
        @pdef:     snumber of segments to split the file into
        @ptype:    {Integer}

        @return: {List} of {File}
        '''
        return self._file.split(ini_separator = '>',
                                sections = sections)

    def load(self):
        '''
        Uploads to memory all the sequences from the file.

        '''
        if self.is_loaded:
            return

        self._sequences = [None] * self._total_sequences
        self.file.open()
        for line in self.file.read():
            if line.startswith('>'):
                seqID = line.lstrip('>').strip()
                s     = Sequence(sequence_id = seqID)
                self._sequences[self._sequenceID[seqID]] = s
            elif len(line.strip()) > 0:
                self._sequences[self._sequenceID[seqID]].append(line.strip())
        self.file.close()
        self._loaded = True

    def live_show(self):
        '''
        Yields the different sequences in the file without actually storing
        them to memory.

        @yields: {Sequence}
        '''
        if self.is_loaded:
            for s in self.sequences:
                yield s

        else:
            n, s = 0, ''
            self.file.open()
            for line in self.file.read():
                if line.startswith('>'):
                    if n > 0:
                        yield s
                    n += 1
                    s = Sequence(sequence_id = line.lstrip('>').strip())
                elif len(line.strip()) > 0:
                    s.append(line.strip())
            self.file.close()
            yield s

    def sequence_lengths(self):
        '''
        Returns a dictionary with the length of each sequence

        @return: {Dictionary}
        '''
        seq = {}
        for s in self.live_show:
            seq[s.id] = len(s)
        return seq

    def retrieve(self, sequence_ids, all_but = False, prefix_size = None):
        '''
        Get specific sequences from the FASTA file.

        @param:    sequence_ids
        @pdef:     sequence identifier(s)
        @ptype:    {String}, {List} or {Set}

        @param:    all_but
        @pdef:     Flag. Instead of retrieving the given ids, we retrieve all
                   except the given ids.
        @pdefault: _False_
        @ptype:    {Boolean}

        @param:    prefix_size
        @pdef:     maximum characters for the prefix. If _None_, all the
                   characters are included.
        @pdefault: _None_
        @ptype:    {Integer}

        @raises: {AttributeError} if sequence_ids is not a valid type.
        @return: {List} of {Sequence}
        '''
        info = 'Skipping sequence {0}' if all_but else 'Getting sequence {0}'

        if not isinstance(sequence_ids, (list, set)):
            SBIg.alert('debug', self, info.format(sequence_ids))
        else:
            SBIg.alert('debug', self, [info.format(x) for x in sequence_ids])

        if isinstance(sequence_ids, basestring):
            sequence_ids = set([sequence_ids])
        if isinstance(sequence_ids, list):
            sequence_ids = set(sequence_ids)
        if isinstance(sequence_ids, set):
            sequences = []
            for s in self.live_show():
                seq_id = s.id if prefix_size is None else s.id[:prefix_size]
                if seq_id in sequence_ids and not all_but:
                    sequences.append(s)
                if seq_id not in sequence_ids and all_but:
                    sequences.append(s)
            return sequences
        else:
            raise AttributeError('sequence_ids must be a string, list or set.')

    def subset(self, sequence_ids, new_fasta_file, all_but = False,
               prefix_size = None, index = False, force = None):
        '''
        Creates a new {Fasta} with the requested subset of sequences.

        @param:    sequence_ids
        @pdef:     sequence identifier(s)
        @ptype:    {String}, {List} or {Set}

        @param:    new_fasta_file
        @pdef:     name of the new fasta file
        @ptype:    {String}

        @param:    all_but
        @pdef:     Flag. Instead of retrieving the given ids, we retrieve all
                   except the given ids.
        @pdefault: _False_
        @ptype:    {Boolean}

        @param:    prefix_size
        @pdef:     maximum characters for the prefix. If _None_, all the
                   characters are included.
        @pdefault: _None_
        @ptype:    {Integer}

        @param:    index
        @pdef:     create the index file also, in case it does exist
        @pdefault: _False_
        @ptype:    {Boolean}

        @param:    force
        @pdef:     overwrite previous files with the same name
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @raises: {AttributeError} if sequence_ids is not a valid type.
        @return: {Fasta}
        '''
        sequences  = self.retrieve(sequence_ids, all_but, prefix_size)
        fasta_file = Fasta.build_multifasta(new_fasta_file, sequences, force)
        if self.has_index and index:
            idxfile = File(self.index_file)
            newidx  = File(fasta_file.file.full + '.idx', 'w')
            seqids  = set(fasta_file.sequence_identifiers)
            for idx in idxfile.read():
                if idx.split()[0].strip('>') in seqids:
                    newidx.write(idx)
            idxfile.close()
            newidx.close()
            fasta_file.index_file = newidx.full
        return fasta_file

    def remove(self, other, new_fasta_file, force = None):
        '''
        Two {Fasta} objects can be subtracted. In that case, sequences from the
        first {Fasta} who appear in the second {Fasta} are removed. Sequences are
        identified only by sequence identifier.

        @param:    other
        @pdef:     {Fasta} containing the sequences to remove from the {Fasta}
        @ptype:    {Fasta}

        @param:    new_fasta_file
        @pdef:     name of the new fasta file
        @ptype:    {String}

        @param:    force
        @pdef:     overwrite previous files with the same name
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @return: {Fasta}
        '''
        remove_seq = set(other.sequence_identifiers)
        sequences  = []

        for seq in self.live_show():
            if seq.id not in remove_seq:
                SBIg.alert('verbose', self, 'ACCEPT sequence {0}'.format(seq.id))
                sequences.append(seq)
            else:
                SBIg.alert('verbose', self, 'REJECT sequence {0}'.format(seq.id))
        return Fasta.build_multifasta(new_fasta_file, sequences, force)

    def reduce(self, new_fasta_file, list_file, force = None):
        '''
        Reduces the {Fasta} by removing identical sequences.

        @param:    new_fasta_file
        @pdef:     name of the new fasta file
        @ptype:    {String}

        @param:    list_file
        @pdef:     name of the repetition list file
        @ptype:    {String}

        @param:    force
        @pdef:     overwrite previous files with the same name
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @return: {Fasta} and {File} with the list of identical sequences.
        '''
        seq_md5   = {}
        sequences = []
        for seq in self.live_show():
            md5 = seq.md5
            if not md5 in seq_md5:
                sequences.append(seq)
                seq_md5.setdefault(md5, [])
            else:
                SBIg.alert('debug', self, '{0} repeats of {1}'.format(seq.id, seq_md5[md5][0]))
            seq_md5[md5].append(seq.id)
        fasta = Fasta.build_multifasta(new_fasta_file, sequences, force)
        listfile = File(list_file, 'w')
        for md5 in seq_md5:
            listfile.write('\t'.join(seq_md5[md5]) + '\n')
        listfile.close()

        return fasta, listfile

    @staticmethod
    def build(file_name, sequence_id, sequence, force = None):
        '''
        Creates a Fasta object and a FASTA file from a sequence.

        @param:    file_name
        @pdef:     name of the fasta file (with path, if necessary)
        @ptype:    {String}

        @param:    sequence_id
        @pdef:     name of the sequence
        @ptype:    {String}

        @param:    sequence
        @pdef:     sequence
        @ptype:    {String} or {List}

        @param:    force
        @pdef:     overwrite previous files with the same name
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @return: {Fasta}
        '''
        newFasta = File(file_name, 'w', overwrite = force)
        newSeq   = Sequence(sequence_id = sequence_id, sequence = sequence)
        newFasta.write(newSeq.format('FASTA'))
        newFasta.close()
        return Fasta(fasta_file = newFasta.full, auto_load = 0)

    @staticmethod
    def build_multifasta(file_name, sequence_list, force = None):
        '''
        Creates a Fasta object and a FASTA file. For multiple sequences.

        @param:    file_name
        @pdef:     name of the fasta file (with path, if necessary)
        @ptype:    {String}

        @param:    sequence_list
        @pdef:     list of sequences to create the FASTA from.
        @ptype:    {List} or {Set} of {Sequence}

        @param:    force
        @pdef:     overwrite previous files with the same name
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @return: {Fasta}
        '''
        newFasta = File(file_name, 'w', overwrite = force)
        for sequence in sequence_list:
            newFasta.write(sequence.format('FASTA') + '\n')
        newFasta.close()
        fasta_file = Fasta(fasta_file = newFasta.full, auto_load = 0)
        return fasta_file

    ###################
    # PRIVATE METHODS #
    ###################
    def _check_multifasta(self):
        '''
        Counts the total number of sequences in the FASTA file.
        Adds all the identifiers to the object.

        '''
        counter = 0
        self.file.open()
        for line in self.file.read():
            if line.startswith('>'):
                self._sequenceID[line.lstrip('>').strip()] = counter
                counter += 1
        self.file.close()
        self._total_sequences = len(self._sequenceID)

        if self._total_sequences <= self._auto_load:
            self.load()

    def _check_index(self):
        '''
        Checks if a idx file associated to the {Fasta} exists.
        Index files are named as the {Fasta} plus the '.idx' termination
        '''
        if os.path.isfile(self.file.full + '.idx'):
            self._index_file = self.file.full + '.idx'

    #################
    # CLASS METHODS #
    #################
    def __len__(self):
        return self._total_sequences

    def __getitem__(self, key):
        try:
            int(key)
        except:
            raise TypeError
        return self._sequences[int(key)]

    def __iter__(self):
        for s in self._sequences:
            yield s

    def __repr__(self):
        return '<{0.__class__.__name__}: [{0.file}]>'.format(self)

    def __str__(self):
        text = []
        for seq in self.live_show():
            text.append(seq.format('FASTA'))
        return "\n".join(text)
