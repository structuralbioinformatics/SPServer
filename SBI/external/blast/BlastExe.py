'''
@file: BlastExe.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2013

@ [oliva's lab](http://sbi.imim.es)

@class: BlastExe
@class: BlastParser
@class: BlastError
'''
import os
import time
import re
import ConfigParser
from abc import ABCMeta

from bs4 import BeautifulSoup

from SBI.external import ExternalExe
from SBI.beans    import Executable
from SBI.beans    import File
from SBI.beans    import Path
from SBI.sequence import Fasta
from SBI          import SBIglobals as SBIg
from BlastResult  import BlastResult
from BlastResult  import BlastHeader
from BlastHit     import BlastHit


class BlastExe(ExternalExe):
    '''
    Executes [BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi).
    See the [user guide](http://www.ncbi.nlm.nih.gov/books/NBK1763/) to learn
    which parameters to add.

    In case a multi fasta is given as database but it is not formatted,
    it will automatically call the formatting of the file.
    '''
    _DBFORMATER = None
    _EXEC_TYPES = {'protein':    ['.phr', '.pin', '.psq'],
                   'nucleotide': ['.nhr', '.nin', '.nsq']}
    _HIT_ID_FRM = frozenset(['single', 'double', 'all'])
    _BLOCKED    = frozenset(['-db', '-outfmt'])

    def __init__(self, database, search_type = 'protein', self_hit = False,
                 hitid_format = 'single', overwrite = None, clean = True):
        '''
        @param:    database
        @pdef:     database to blast upon.
        @ptype:    {String}

        @param:    search_type
        @pdef:     type of database search.
        @pdefault: 'protein'
        @ptype:    {String}

        @@param:   self_hit
        @pdef:     when _True_ if the query is found in the database, it is
                   retrieved.
        @pdefault: _False_
        @ptype:    {Boolean}

        @param:    hitid_format
        @pdef:     format of the name of the hit. If given a wrong option,
                   it defaults to 'single'
        @pdefault: 'single'
        @poptions: 'single' -> first word of the name,
                   'double' -> first two words of the hit name,
                   'all'    -> all the text in the hit name
        @ptype:    {String}

        @param:    overwrite
        @pdef:     For writing actions. Decides whether it can overwrite an
                   existing file.
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @param:    clean
        @pdef:     remove the temporary files after the data is read.
        @pdefault: _True_
        @pclash:   if _SBIglobals.debug_ is _True_, clean is _False_
        @ptype:    {Boolean}

        @raises: {BlastError}
        '''
        self._error = BlastError()

        self._search_type = self._check_execution_type(search_type)

        # Blast executable configuration
        if self._EXE is None:
            self._set_default_executable('blast')
            BlastExe._DBFORMATER = BlastExe._CONFIG.get('blast', 'dbformatexe')

        # Database Configuration
        self._idx      = None
        self._database = self._check_database(os.path.abspath(database))

        # Optional execution parameters
        self._parameters = []

        msg  = ['New Blast Executable created.', ]
        msg.append('Blast executable at {0}'.format(self._EXE.full_executable))
        SBIg.alert('debug', self, msg)

        self._self_hit     = self_hit
        self._hitid_format = self._check_hitid_format(hitid_format)
        self._clean_files  = clean

        # Local overwrite takes precedence over Global overwrite
        self._overwrite = SBIg.decide_overwrite(overwrite)

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def database(self):
        '''
        Database to blast against.

        @returns: {String}
        '''
        return self._database

    @property
    def database_index(self):
        '''
        Index to correct database positions (if exists)

        @return: {String}
        '''
        return self._idx

    @property
    def self_hit(self):
        '''
        Capture self hit in the database.

        @returns: {Boolean}
        '''
        return self._self_hit

    @property
    def hitid_format(self):
        '''
        Amount of hit name to store

        @returns: {String}
        '''
        return self._hitid_format

    @property
    def overwrite(self):
        '''
        Overwrite old files.

        @returns: {Boolean}
        '''
        return self._overwrite

    @property
    def clean_files(self):
        '''
        Remove temporary files after process is finished.

        @returns: {Boolean}
        '''
        return self._clean_files

    ###########
    # METHODS #
    ###########
    @staticmethod
    def dynamic(executable, path, dbformater):
        '''
        Manually set the values for blast path, executable and db format.

        @param:    executable
        @pdef:     name of the executable file
        @ptype:    {String}

        @param:    path
        @pdef:     path to the executable file
        @ptype:    {String}

        @param:    dbformater
        @pdef:     name of the database format executable
        @ptype:    {String}
        '''
        BlastExe._set_dynamic_executable(executable, path)
        BlastExe._DBFORMATER = dbformater

    def add_attribute(self, attribute_value, attribute_id):
        '''
        Adds specific attributes to the blast execution.

        @param:    attribute_value
        @pdef:     value of the attribute to add.
        @ptype:    {String}. Transformed otherwise.

        @param:    attribute_id
        @pdef:     label of the attribute to add.
        @pclash:   labels 'db' and 'outfmt' are blocked, as they need to be
                   specified when creating the instance.
        @ptype:    {String}
        '''
        if attribute_id in BlastExe._BLOCKED:
            raise self._error.blocked_parameter(BlastExe._BLOCKED)

        self._parameters.append([str(attribute_value), attribute_id])

    def execute_query_seq(self, sequenceID = None, sequence          = None,
                          blast_input_file = None, blast_output_file = None,
                          work_directory   = os.getcwd()):
        '''
        Execute BLAST given a query sequence.

        @param:    sequenceID
        @pdef:     name of the query sequence.
        @pdefault: 'QuerySequence'
        @pclash:   If sequence is not provided, it assumes that the sequenceID
                   belongs to a protein in the database and, thus, it searches
                   for it. Either sequenceID or sequence needs to be provided.
        @ptype:    {String}

        @param:    sequence
        @pdef:     query sequence.
        @pdefault: _None_
        @pclash:   Either sequenceID or sequence needs to be provided.
        @ptype:    {String}

        @param:    blast_input_file
        @pdef:     name of the temporary fasta file to use as BLAST input.
        @pdefault: job.pid + clock + .tmp.fa
        @ptype:    {String}

        @param:    blast_output_file
        @pdef:     name of the temporary BLAST output file.
        @pdefault: job.pid + clock + .blast.xml.out
        @ptype:    {String}

        @param:    work_directory
        @pdef:     Directory to which the temporary files will be created.
        @pdefault: Current working directory.
        @ptype:    {String}

        @raises: {AttributeError} if neither sequenceID nor sequence are
                  provided or if sequenceID is a list of sequence names.
        @raises: {BlastError} in BLAST execution or output parsing errors.

        @returns: {BlastResult}
        '''
        if sequenceID is None and sequence is None:
            msg = 'Either a sequence or sequenceID is needed to perform the blast.'
            raise AttributeError(msg)

        if isinstance(sequenceID, (list, set, tuple)):
            msg = 'Blasts can only be executed one at a time due to XML output restrictions.'
            raise AttributeError(msg)

        sequenceID = 'QuerySequence' if sequenceID is None else sequenceID

        # Given only a code implies that the protein of interest is in the
        # database itself
        if sequence is None:
            grabbedSequence = self._database.retrieve(sequenceID)
            sequenceID      = grabbedSequence[0].id
            sequence        = grabbedSequence[0].sequence

        # All the sequence is unknown, it will crash blast
        if len(re.sub(r'[Xx]', '', sequence)) == 0:
            SBIg.warn(self, 'Created an empty BlastResult.')
            return BlastResult(query_name     = sequenceID,
                               query_sequence = sequence)

        Path.mkdir(work_directory)
        file_prefixes = ".".join([str(os.getpid()), str(int(time.clock()*100000))])
        file_prefixes = os.path.join(work_directory, file_prefixes)
        tmp_input     = file_prefixes + ".tmp.fa"
        tmp_output    = file_prefixes + ".blast.xml.out"

        tmp_input  = tmp_input  if blast_input_file  is None else os.path.join(work_directory, blast_input_file)
        tmp_output = tmp_output if blast_output_file is None else os.path.join(work_directory, blast_output_file)

        QueryFasta = Fasta.build(file_name = tmp_input, sequence_id = sequenceID,
                                 sequence  = sequence,  force       = True)

        self._execute(input_file = QueryFasta, output_file = tmp_output)

        blast_result = self._parse_blast(sequence, tmp_output)

        self._clean([tmp_input, tmp_output])

        return blast_result

    def execute_query(self, query_file, blast_output_file = None,
                      work_directory = os.getcwd()):
        '''
        Execute BLAST given a query sequence.

        @param:    query_file
        @pdef:     Fasta file with the query sequence.
        @pdefault: 'QuerySequence'
        @ptype:    {String} or {File} or {Fasta}

        @param:    blast_output_file
        @pdef:     name of the temporary BLAST output file.
        @pdefault: query_file.prefix + job.pid + .blast.xml.out
        @ptype:    {String}

        @param:    work_directory
        @pdef:     Directory to which the temporary files will be created.
        @pdefault: Current working directory.
        @ptype:    {String}

        @raises: {AttributeError} if query_file is multi-fasta.
        @raises: {BlastError} in BLAST execution or output parsing errors.

        @returns: {BlastResult}
        '''
        if isinstance(query_file, basestring) or isinstance(query_file, File):
            newFasta = Fasta(fasta_file = query_file)
        elif isinstance(query_file, Fasta):
            newFasta = query_file

        if newFasta.is_multifasta:
            msg = 'Blasts can only be executed one at a time due to XML output restrictions.'
            raise AttributeError(msg)

        # All the sequence is unknown, it will crash blast
        newFasta.load()
        query_sequence = newFasta.sequence
        if len(re.sub(r'[Xx]', '', query_sequence.sequence)) == 0:
            SBIg.warn(self, 'Created an empty BlastResult.')
            return BlastResult(query_name     = query_sequence.id,
                               query_sequence = query_sequence.sequence)

        Path.mkdir(work_directory)
        file_prefixes = ".".join([newFasta.file.prefix, str(os.getpid())])
        file_prefixes = os.path.join(work_directory, file_prefixes)
        tmp_output    = file_prefixes + ".blast.xml.out"

        tmp_output = tmp_output if blast_output_file is None else os.path.join(work_directory, blast_output_file)

        self._execute(input_file = newFasta, output_file = tmp_output)

        blast_result = self._parse_blast(newFasta.sequence.sequence, tmp_output)

        self._clean([tmp_output, ])

        return blast_result

    def clean_optional_parameters(self):
        '''
        Deletes all optional parameters added to the {BlastExe}.
        '''
        self._parameters = []

    ###################
    # PRIVATE METHODS #
    ###################
    def _execute(self, input_file, output_file):
        '''
        Executes BLAST.

        @param:    input_file
        @pdef:     file with the sequence to blast. (SINGLE FASTA)
        @ptype:    {Fasta}

        @param:    output_file
        @pdef:     name of the blast output file.
        @ptype:    {String}

        @raises: {BlastError}
        '''
        if not os.path.isfile(output_file) or self.overwrite:
            self._EXE.backup_command()

            # Adding fixed blast parameters
            self._EXE.add_attribute(self._database.file.full, '-db')
            self._EXE.add_attribute('5', '-outfmt')
            self._EXE.add_parameter('-lcase_masking')

            # Adding optional parameters
            for parameter in self._parameters:
                self._EXE.add_attribute(parameter[0], parameter[1])

            # Adding input-output
            self._EXE.add_attribute(input_file.file.full, '-query')
            self._EXE.add_attribute(output_file,          '-out')

            try:
                self._EXE.execute()
                self._EXE.restore_command()
            except SystemError as e:
                default_warning = bool(re.search(self._error.default_blast_warning(),  str(e)))
                seleno_warning  = bool(re.search(self._error.selenocysteine_warning(), str(e)))
                repeat_warning  = bool(re.search(self._error.repetitive_sequence_warning(), str(e)))
                if not default_warning and not seleno_warning and not repeat_warning:
                    raise self._error.blast_execution_failed(e)

    def _clean(self, files):
        '''
        If clean_files is _True_ and not in debug mode, it removes the
        temporary files.

        @param:    files
        @pdef:     list of files to remove
        @ptype:    {List}
        '''
        if self.clean_files and not SBIg.debug:
            for temp_file in files:
                os.unlink(temp_file)

    def _parse_blast(self, query_sequence, blast_output_file):
        '''
        Calls the parsing of the blast output.

        @param:    query_sequence
        @pdef:     sequence that has been searched against the database.
        @ptype:    {String}

        @param:    blast_output_file
        @pdef:     output file of blast.
        @ptype:    {String}
        '''
        return BlastParser.parse(query_sequence, blast_output_file,
                                 self._self_hit, self._hitid_format)

    def _check_database(self, database):
        '''
        Ensures that the given database to blast upon exists and that it is
        formated for blast.
        I also sets the index file for the database if it exists.

        @param:    database
        @pdef:     database to blast upon.
        @ptype:    {String}

        @returns: {Fasta} object pointed to the database fasta file.
        '''
        # Database file does not exist
        if not os.path.isfile(database):
            return self._error.database_does_not_exist(database)

        # Database is not formated (if dbformatexe is added in the
        # configuration path it will be auto-formated)
        formatdb_files = []
        formatdb_sufix = BlastExe._EXEC_TYPES[self._search_type]

        for sufix in formatdb_sufix:
            if not os.path.isfile(database + sufix):
                formatdb_files.append(database + sufix)

        if len(formatdb_files) > 0:
            try:
                self._format_database(database)
            except ConfigParser.NoOptionError as e:
                raise self._error.no_blast_format_exe(e)
            except SystemError as e:
                raise self._error.wrong_db_format(database, e)

        idx       = os.path.abspath(database) + ".idx"
        self._idx = idx if os.path.isfile(idx) else None

        return Fasta(fasta_file = database)

    def _format_database(self, database):
        '''
        Executes the blast script to format the database.

        @param:    database
        @pdef:     database to blast upon.
        @ptype:    {String}
        '''
        SBIg.warn(self, 'Formating {0} for blast.'.format(database))
        dbexe  = Executable(executable = BlastExe._DBFORMATER,
                            path       = self._EXE.path)

        dbexe.add_attribute(database,               '-in')
        dbexe.add_attribute(self._search_type[0:4], '-dbtype')

        SBIg.alert('debug', self, 'Executing command {0}\n'.format(dbexe))

        dbexe.execute()

    def _check_hitid_format(self, hitid):
        '''
        Ensures that the hitid_format requested is amongst the available.

        @param:    hitid
        @pdef:     format of the name of the hit. If given a wrong option,
                   it defaults to 'single'
        @pdefault: 'single'
        @poptions: 'single' -> first word of the name,
                   'double' -> first two words of the hit name,
                   'all'    -> all the text in the hit name
        @ptype:    {String}

        @returns: {String}
        '''
        return 'single' if hitid not in BlastExe._HIT_ID_FRM else hitid

    def _check_execution_type(self, search_type):
        '''
        Ensures that the blast is performed in a known search mode.

        @param:    search_type
        @pdef:     type of blast search
        @poptions: 'protein', 'nucleotide'
        @ptype:    {String}

        @return: {String}
        '''
        if search_type not in BlastExe._EXEC_TYPES:
            raise self._error.incorrect_search_type(search_type)
        return search_type


class BlastParser(object):
    '''
    Processes a blast xml formated output into a {BlastResult} object.
    '''
    __metaclass__ = ABCMeta

    @staticmethod
    def parse(query_sequence, blast_output_file, self_hit, hitid_format):
        '''
        Processes a blast xml formated output into a {BlastResult} object.

        @param:    query_sequence
        @pdef:     sequence of the query protein/nucleotide.
        @ptype:    {String}

        @param:    blast_output_file
        @pdef:     output file from BLAST.
        @ptype:    {String}

        @param:   self_hit
        @pdef:     when _True_ if the query is found in the database, it is
                   retrieved.
        @pdefault: _False_
        @ptype:    {Boolean}

        @param:    hitid_format
        @pdef:     format of the name of the hit. If given a wrong option,
                   it defaults to 'single'
        @pdefault: 'single'
        @poptions: 'single' -> first word of the name,
                   'double' -> first two words of the hit name,
                   'all'    -> all the text in the hit name
        @ptype:    {String}

        @raises: {BlastError} if there are problems while parsing the XML file.
        @returns: {BlastResult}
        '''
        f = File(blast_output_file)
        s = BeautifulSoup(f.read())

        h = BlastHeader(version    = str(s.find('blastoutput_version').string),
                        matrix     = str(s.find('parameters_matrix').string),
                        gap_open   = int(s.find('parameters_gap-open').string),
                        gap_extend = int(s.find('parameters_gap-extend').string),
                        database   = str(s.find('blastoutput_db').string),
                        self_hit   = self_hit)
        b = BlastResult(query_name     = str(s.find('blastoutput_query-def').string),
                        query_sequence = query_sequence,
                        header         = h)

        SBIg.alert('debug', BlastParser(), b.str_blast_header())

        error_bool = False
        error_str  = []
        for iteration in s.find_all('iteration'):
            iternum = int(iteration.find('iteration_iter-num').string)
            for hit in iteration.find_all('hit'):
                hit_name  = BlastParser.hit_name(str(hit.find('hit_def').string), hitid_format)
                hit_lenth = int(hit.find("hit_len").string)
                for subhit in hit.find_all("hsp"):
                    data = BlastParser.parse_subhit(subhit)
                    r = BlastHit(hit            = [hit_name, hit_lenth],
                                 sequences      = [data['qs'], data['hs'], data['sc']],
                                 sequence_inits = [data['qp'], data['hp']],
                                 iteration      = iternum,
                                 stats          = [data['hi'], data['h+'],
                                                   data['hg'], data['ev']])
                    if not BlastParser.same_query_hit_names(b.query, hit_name, self_hit):
                        dbug_info = 'Added hit {0} in iteration {1}'
                        SBIg.alert('debug', BlastParser(), dbug_info.format(hit_name, r.iteration))
                        b.add_hit(r)
                        if not r.are_segments_ok:
                            error_bool = True
                            error_str.append("Check the alignment's fragmentation")
                            error_str.append("for the query %s with %s\n".format(b.query, hit_name))
                            error_str.append("{0}\n".format(r))
        b.set_last_iteration()
        if error_bool:
            SBIg.warn(BlastParser(), error_str)
            be = BlastError()
            raise be.parse_error()
        return b

    @staticmethod
    def parse_subhit(subhit):
        '''
        Extracts the relevant data from a subhit representation.

        @param:    subhit
        @pdef:     XML section that represents a subhit
        @ptype:    {bs4.element.Tag}

        @returns: {Dictionary}
        '''
        data = {}

        gaps       = subhit.find("hsp_gaps")
        data['hg'] = int(gaps.string) if gaps is not None else 0
        idnt       = subhit.find("hsp_identity")
        data['hi'] = int(idnt.string) if idnt is not None else 0
        posi       = subhit.find("hsp_positive")
        data['h+'] = int(posi.string) if posi is not None else 0

        data['ev'] = float(subhit.find("hsp_evalue").string)
        data['al'] = int(subhit.find("hsp_align-len").string)
        data['qs'] = str(subhit.find("hsp_qseq").string).strip()
        data['hs'] = str(subhit.find("hsp_hseq").string).strip()
        data['qp'] = int(subhit.find("hsp_query-from").string)
        data['hp'] = int(subhit.find("hsp_hit-from").string)
        data['sc'] = str(subhit.find("hsp_midline").string).strip()

        return data

    @staticmethod
    def same_query_hit_names(query, hit, self_hit):
        '''
        Compares the query and hit names to check if they are the same.
        Comparison is based on the first section of the name.
        If self_hit is _True_ it always assumes that they are NOT the same.

        @param:    query
        @pdef:     query name.
        @ptype:    {String}

        @param:    hit
        @pdef:     hit name.
        @ptype:    {String}

        @param:    self_hit
        @pdef:     tag to decide whether to compare the names.
        @ptype:    {Boolean}
        '''
        if self_hit:
            return False
        return query.split()[0].strip() == hit.split()[0].strip()

    @staticmethod
    def hit_name(full_name, hitid_format):
        '''
        Processes the name of the hit depending on the selected format.

        @param:    full_name
        @pdef:     original name of the hit
        @ptype:    {String}

        @param:    hitid_format
        @pdef:     format expected for the hit name
        @ptype:    {String}

        '''
        if hitid_format   == 'single':
            return full_name.split()[0].strip()
        elif hitid_format == 'double':
            return " ".join(full_name.split()[:2]).strip()
        else:
            return full_name


class BlastError(Exception):
    '''
    Manages different errors that can occur during the blast execution or
    parsing.

    '''
    _MSG     = ''
    _EXE     = False
    _NOBLAST = 'Blast has NOT been executed.'

    def __init__(self):
        pass

    def incorrect_search_type(self, qsearch):
        BlastError._EXE = False
        BlastError._MSG = '{0} is not a valid blast search method.'.format(qsearch)
        return self

    def database_does_not_exist(self, database):
        BlastError._EXE = False
        BlastError._MSG = '{0} is not a file.'.format(database)
        return self

    def wrong_db_format(self, database, e):
        BlastError._EXE = False
        BlastError._MSG = '{0} could not be formated for blast.'.format(database)
        BlastError._MSG += '\n' + str(e)
        return self

    def no_blast_format_exe(self, e):
        BlastError._EXE = False
        BlastError._MSG = 'Executable to format blast is not specified.'
        BlastError._MSG += '\n' + str(e)
        return self

    def blast_execution_failed(self, e):
        BlastError._EXE = False
        BlastError._MSG = 'Blast execution failed.'
        BlastError._MSG += '\n' + str(e)
        return self

    def blocked_parameter(self, param):
        BlastError._EXE = False
        BlastError._MSG = '{0} are blocked.'.format(param)
        return self

    def parse_error(self):
        BlastError._EXE = True
        BlastError._MSG = 'An error occurred while parsing the BLAST output.'
        return self

    def corrected_hits(self):
        BlastError._EXE = True
        BlastError._MSG = 'Hits have already been corrected.'
        return self

    def default_blast_warning(self):
        msg  = 'Warning: Composition-based score adjustment conditioned on '
        msg += 'sequence properties and unconditional composition-based score '
        msg += 'adjustment is not supported with PSSMs, resetting to default '
        msg += 'value of standard composition-based statistics'
        return msg

    def selenocysteine_warning(self):
        return 'Selenocysteine \(U\) at position'

    def repetitive_sequence_warning(self):
        return 'bad probabilities from'

    def __str__(self):
        d  = BlastError._MSG + '\n'
        d += BlastError._NOBLAST if not BlastError._EXE else ''
        return d
