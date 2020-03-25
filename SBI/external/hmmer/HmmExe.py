'''
@file: HmmExe.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2014

@ [oliva's lab](http://sbi.imim.es)

@class: HmmExe
@class: HmmParser
@class: HmmError
'''
import os
import re
import time
import subprocess
import ConfigParser
from abc import ABCMeta

from SBI.external import ExternalExe
from SBI.beans    import Executable
from SBI.beans    import File
from SBI.sequence import Fasta
from HmmResult    import HmmResult
from HmmHit       import HmmHit
from SBI          import SBIglobals as SBIg


class HmmExe(ExternalExe):
    '''
    Executes [HMMER](http://hmmer.janelia.org/).
    See the [user guide](ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/Userguide.pdf)
    to learn which parameters to add.

    In case a multi fasta is given as database but it is not formatted,
    it will automatically call the formatting of the file.
    '''
    _DBFORMATER = None
    _FORMATSDB  = ['.h3f', '.h3i', '.h3m', '.h3p']
    _BLOCKED    = frozenset(['--domtblout'])

    def __init__(self, database, overwrite = None, clean = True):
        '''
        @param:    database
        @pdef:     database to blast upon.
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

        @raises: {HmmError}
        '''
        self._error    = HmmError()

        #hmmer executable configuration
        if HmmExe._EXE is None:
            self._set_default_executable('hmmer')
            HmmExe._DBFORMATER = HmmExe._CONFIG.get('hmmer', 'dbformatexe')

        # Local overwrite takes precedence over Global overwrite
        self._overwrite   = SBIg.decide_overwrite(overwrite)

        self._database    = self._check_database(os.path.abspath(database))
        self._clean_files = clean

        #Optional execution parameters
        self._parameters  = {'attr': [], 'flag': []}

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def database(self):
        '''
        Database to hmmer against.

        @returns: {String}
        '''
        return self._database

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
        Manually set the values for hmmer path, executable and db format.

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
        HmmExe._set_dynamic_executable(executable, path)
        HmmExe._DBFORMATER = dbformater

    def add_attribute(self, attribute_value, attribute_id):
        '''
        Adds specific attributes to the hmmer execution.

        @param:    attribute_value
        @pdef:     value of the attribute to add.
        @ptype:    {String}. Transformed otherwise.

        @param:    attribute_id
        @pdef:     label of the attribute to add.
        @pclash:   label 'domtblout' is blocked, as they need to be
                   specified when creating the instance.
        @ptype:    {String}
        '''
        if attribute_id in HmmExe._BLOCKED:
            raise self._error.blocked_parameter(HmmExe._BLOCKED)

        self._parameters['attr'].append([str(attribute_value), attribute_id])

    def add_parameter(self, attribute_id):
        '''
        Adds specific flag to the hmmer execution.

        @param:    attribute_id
        @pdef:     label of the attribute to add.
        @pclash:   label 'domtblout' is blocked, as they need to be
                   specified when creating the instance.
        @ptype:    {String}
        '''
        if attribute_id in HmmExe._BLOCKED:
            raise self._error.blocked_parameter(HmmExe._BLOCKED)

        self._parameters['flag'].append(attribute_id)

    def execute_query_seq(self, sequenceID = None, sequence          = None,
                          hmmer_input_file = None, hmmer_output_file = None,
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

        @param:    hmmer_input_file
        @pdef:     name of the temporary fasta file to use as HMMER input.
        @pdefault: job.pid + clock + .tmp.fa
        @ptype:    {String}

        @param:    blast_output_file
        @pdef:     name of the temporary HMMER output file.
        @pdefault: job.pid + clock + .hmmer.out
        @ptype:    {String}

        @param:    work_directory
        @pdef:     Directory to which the temporary files will be created.
        @pdefault: Current working directory.
        @ptype:    {String}

        @raises: {AttributeError} if neither sequenceID nor sequence are
                  provided or if sequenceID is a list of sequence names.
        @raises: {HmmError} in HMMER execution or output parsing errors.

        @returns: {HmmResult}
        '''
        if sequenceID is None or sequence is None:
            msg = 'Both a sequence and sequenceID are needed to perform the blast.'
            raise AttributeError(msg)

        if isinstance(sequenceID, (list, set, tuple)):
            msg = 'Blasts can only be executed one at a time due to parse restrictions.'
            raise AttributeError(msg)

        # All the sequence is unknown, it will crash blast
        if len(re.sub(r'[Xx]', '', sequence)) == 0:
            SBIg.warn(self, 'Created an empty BlastResult.')
            return HmmResult(query_name     = sequenceID,
                             query_sequence = sequence,
                             database       = self.database)

        file_prefixes = ".".join([str(os.getpid()), str(int(time.clock()*100000))])
        file_prefixes = os.path.join(work_directory, file_prefixes)
        tmp_input     = file_prefixes + ".tmp.fa"
        tmp_output    = file_prefixes + ".blast.xml.out"

        tmp_input  = tmp_input  if hmmer_input_file  is None else hmmer_input_file
        tmp_output = tmp_output if hmmer_output_file is None else hmmer_output_file

        QueryFasta = Fasta.build(file_name = tmp_input, sequence_id = sequenceID,
                                 sequence  = sequence,  force       = True)

        self._execute(input_file = QueryFasta, output_file = tmp_output)

        hmmer_result = self._parse_hmmer(sequenceID, sequence, tmp_output)

        self._clean([tmp_input, tmp_output])

        return hmmer_result

    def execute_query(self, query_file, hmmer_output_file = None,
                      work_directory = os.getcwd()):
        '''
        Execute HMMER given a query sequence.

        @param:    query_file
        @pdef:     Fasta file with the query sequence.
        @pdefault: 'QuerySequence'
        @ptype:    {String} or {File} or {Fasta}

        @param:    hmmer_output_file
        @pdef:     name of the temporary HMMER output file.
        @pdefault: query_file.prefix + job.pid + .hmmer.out
        @ptype:    {String}

        @param:    work_directory
        @pdef:     Directory to which the temporary files will be created.
        @pdefault: Current working directory.
        @ptype:    {String}

        @raises: {AttributeError} if query_file is multi-fasta.
        @raises: {HmmError} in HMMER execution or output parsing errors.

        @returns: {HmmResult}
        '''
        if isinstance(query_file, basestring) or isinstance(query_file, File):
            newFasta = Fasta(fasta_file = query_file)
        elif isinstance(query_file, Fasta):
            newFasta = query_file

        if newFasta.is_multifasta:
            msg = 'Hmmer can only be executed one at a time due to parsing restrictions.'
            raise AttributeError(msg)

        #All the sequence is unknown, it will crash blast
        newFasta.load()
        query_sequence = newFasta.sequence
        if len(re.sub(r'[Xx]', '', query_sequence.sequence)) == 0:
            SBIg.warn(self, 'Created an empty BlastResult.')
            return HmmResult(query_name     = query_sequence.id,
                             query_sequence = query_sequence.sequence,
                             database       = self.database)

        file_prefixes = ".".join([newFasta.file.prefix, str(os.getpid())])
        file_prefixes = os.path.join(work_directory, file_prefixes)
        tmp_output    = file_prefixes + ".hmmer.out"

        tmp_output = tmp_output if hmmer_output_file is None else hmmer_output_file

        self._execute(input_file = newFasta, output_file = tmp_output)

        hmm_result = self._parse_hmmer(newFasta.sequence.id,
                                       newFasta.sequence.sequence, tmp_output)

        self._clean([tmp_output, ])

        return hmm_result

    def clean_optional_parameters(self):
        '''
        Deletes all optional parameters added to the {HmmExe}.
        '''
        self._parameters = {'attr': [], 'flag': []}

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

        @raises: {HmmError}
        '''
        if not os.path.isfile(output_file) or self.overwrite:
            self._EXE.backup_command()

            #Adding fixed blast parameters
            self._EXE.add_attribute(output_file, '--domtblout')

            # Adding optional parameters
            for parameter in self._parameters['attr']:
                self._EXE.add_attribute(parameter[0], parameter[1])
            for flag in self._parameters['flag']:
                self._EXE.add_parameter(flag)

            # Adding input-output
            self._EXE.add_parameter(self.database)
            self._EXE.add_parameter(input_file.file.full)

            try:
                self._EXE.execute()
                self._EXE.restore_command()
            except SystemError as e:
                raise self._error.hmmer_execution_failed(e)

    def _check_database(self, database):
        '''
        Ensures that the given database to hmmer upon exists and that it is
        formated for hmmer.
        I also sets the index file for the database if it exists.

        @param:    database
        @pdef:     database to hmmer upon.
        @ptype:    {String}

        @returns: {String}, name of the database
        '''
        # Database file does not exist
        if not os.path.isfile(database):
            return self._error.database_does_not_exist(database)

        # Database is not formated (if dbformatexe is added in the
        # configuration path it will be auto-formated)
        formatdb_files = []

        for sufix in HmmExe._FORMATSDB:
            if not os.path.isfile(database + sufix):
                formatdb_files.append(database + sufix)

        if len(formatdb_files) > 0:
            try:
                self._format_database(database)
            except ConfigParser.NoOptionError as e:
                raise self._error.no_hmmer_format_exe(e)
            except SystemError as e:
                raise self._error.wrong_db_format(database, e)

        return database

    def _parse_hmmer(self, query_name, query_sequence, hmmer_output_file):
        '''
        Calls the parsing of the blast output.

        @param:    query_name
        @pdef:     name of the query protein
        @ptype:    {String}

        @param:    query_sequence
        @pdef:     sequence that has been searched against the database.
        @ptype:    {String}

        @param:    hmmer_output_file
        @pdef:     output file of blast.
        @ptype:    {String}
        '''
        return HmmParser.parse(query_name, query_sequence, self.database,
                               hmmer_output_file)

    def _format_database(self, database):
        '''
        Executes the hmmer script to format the database.

        @param:    database
        @pdef:     database to blast upon.
        @ptype:    {String}
        '''
        SBIg.warn(self, 'Formating {0} for hmmer.'.format(database))

        dbexe  = Executable(executable = HmmExe._DBFORMATER,
                            path       = self._EXE.path)

        if self.overwrite:
            dbexe.add_parameter('-f')

        dbexe.add_parameter(database)

        SBIg.alert('debug', self, 'Executing command {0}\n'.format(dbexe))

        dbexe.execute()

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


class HmmParser(object):
    '''
    Processes a cd-hit output into a {HmmResult} object.
    '''
    __metaclass__ = ABCMeta

    @staticmethod
    def parse(query_name, query_sequence, database, hmmer_output_file):
        '''
        Processes a cd-hit output into a {HmmResult} object.

        @param:    query_name
        @pdef:     name of the query protein.
        @ptype:    {String}

        @param:    query_sequence
        @pdef:     sequence of the query protein.
        @ptype:    {String}

        @param:    hmmer_output_file
        @pdef:     output file from BLAST.
        @ptype:    {String}

        @raises: {HmmError} if there are problems while parsing the XML file.
        @returns: {HmmResult}
        '''
        e = HmmError()
        if not os.path.isfile(hmmer_output_file):
            raise e.no_output_file()

        hmm = HmmResult(query_name, query_sequence, database)
        readfile = ["sort", "-gk7", "-gk13", hmmer_output_file]
        p = subprocess.Popen(readfile, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err.strip() == '':
            for line in out.split('\n'):
                if not line.startswith('#') and line.strip() != "":
                    line_data = line.strip().split()
                    new_domain = HmmHit(name      = line_data[3],  accession = line_data[4],
                                        length    = line_data[5],
                                        query_lim = [line_data[17], line_data[18]],
                                        hit_lim   = [line_data[15], line_data[16]],
                                        values    = [line_data[6],  line_data[12]],
                                        repeats   = line_data[10])
                    hmm.add_hit(new_domain)
        else:
            raise e.parse_error()

        return hmm


class HmmError(Exception):
    '''
    Manages different errors that can occur during the blast execution or
    parsing.

    '''
    _MSG     = ''
    _EXE     = False
    _NOBLAST = 'Hmmer has NOT been executed.'

    def __init__(self):
        pass

    def database_does_not_exist(self, database):
        HmmError._EXE = False
        HmmError._MSG = '{0} is not a file.'.format(database)
        return self

    def wrong_db_format(self, database, e):
        HmmError._EXE = False
        HmmError._MSG = '{0} could not be formated for hmmer.'.format(database)
        HmmError._MSG += '\n' + e.message
        return self

    def no_hmmer_format_exe(self, e):
        HmmError._EXE = False
        HmmError._MSG = 'Executable to format hmmer is not specified.'
        HmmError._MSG += '\n' + e.message
        return self

    def hmmer_execution_failed(self, e):
        HmmError._EXE = False
        HmmError._MSG = 'HMMER execution failed.'
        HmmError._MSG += '\n' + e.message
        return self

    def blocked_parameter(self, param):
        HmmError._EXE = False
        HmmError._MSG = '{0} are blocked.'.format(param)
        return self

    def no_output_file(self):
        HmmError._EXE = False
        HmmError._MSG = 'HMMER output not found.'
        return self

    def parse_error(self):
        HmmError._EXE = True
        HmmError._MSG = 'An error occurred while parsing the HMMER output.'
        return self

    def __str__(self):
        d  = HmmError._MSG + '\n'
        d += HmmError._NOBLAST if not HmmError._EXE else ''
        return d
