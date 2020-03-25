'''
@file: file.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   02/2013

@ [oliva's lab](http://sbi.imim.es)

@class: File
@class: FileError
'''

import os
import math
import gzip
import zipfile
import tarfile
import bz2

from . import Butler
SBIg = Butler()


class File(object):
    '''
    Manages certain aspects of the IO process.

    It allows to seamless control some aspects about working with files.
    Namely, it will allow to:
        * work with regular, tar, zip, gzip and bzip2 files.
        * quickly obtain paths, file name, extension...
        * check the viability reading/writing a given file
        * avoid/ignore overwrite an existing file
        * split a file into several
        * get files contained inside (for tar and zip) and extract them

    '''
    EMPTY_ACTION     = frozenset(['n'])
    WRITE_ACTION     = frozenset(['w', 'a', 'ar', 'wb'])
    READ_ACTION      = frozenset(['r', 'rb', 'r|*', 'r|'])
    AVAILABLE_ACTION = WRITE_ACTION.union(READ_ACTION, EMPTY_ACTION)

    TAR_EXTENSIONS    = frozenset(['.tar',    '.tgz', '.tb2', '.tbz2',
                                   '.tar.gz', '.tar.bz2'])
    GZIP_EXTENSIONS   = frozenset(['.gz',  '.tar.gz', '.tgz'])
    BZIP_EXTENSIONS   = frozenset(['.bz2', 'tar.bz2', '.tb2', '.tbz2'])
    ZIP_EXTENSIONS    = frozenset(['.zip'])

    ACTIVE_EXTENSIONS = TAR_EXTENSIONS.union(GZIP_EXTENSIONS,
                                             BZIP_EXTENSIONS,
                                             ZIP_EXTENSIONS)

    def __init__(self, file_name, action = 'r', overwrite = None):
        '''
        @param:    file_name
        @pdef:     File name.
        @ptype:    {String}

        @param:    action
        @pdef:     Whether we are opening a file to read ('r') or write ('w').
        @pdefault: 'r'
        @ptype:    {String}

        @param:    overwrite
        @pdef:     For writing actions. Decides whether it can overwrite an
                   existing file.
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @raises: {FileError} after several checks.
        '''
        self._file  = file_name
        self._error = FileError(file_name)

        SBIg.alert(0, 'Opening: {0}'.format(self.full))

        self._type  = ''
        self._multi = 1
        self._check_file_type()

        SBIg.alert(2, '\tFile is {0}'.format(self._type))
        SBIg.alert(2, '\tContains {0} file'.format(self._multi))

        # The selected action must be valid according to the active extension
        # or simply to a regular file.
        self._action = None
        self._check_action(action.lower())

        SBIg.alert(1, '\tApply Action {0}'.format(self._action))

        # Local overwrite takes precedence over Global overwrite
        self._overwrite = SBIg.decide_overwrite(overwrite)

        ostring = 'overwrite' if self._overwrite else 'not overwrite'
        if not self._action.startswith('r'):
            SBIg.alert(1, '\tMode set to {0}'.format(ostring))

        # Check that the requested action can be performed over that
        # particular file
        self._check_file()

        self._is_open = False
        self._fd      = None

        # For files generated via split
        self._section = None

        self.open()

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def full(self):
        '''
        Full file path.

        @return:  {String}
        '''
        return os.path.abspath(self._file)

    @property
    def relative(self):
        '''
        Path relative to the current working directory

        @return: {String}
        '''
        return os.path.relpath(self.full)

    @property
    def dir(self):
        '''
        Full file directory path.

        @return:  {String}
        '''
        return os.path.split(self.full)[0]

    @property
    def lastdir(self):
        '''
        Last directory name.

        @return:  {String}
        '''
        return os.path.basename(self.dir)

    @property
    def name(self):
        '''
        File name. Without path.

        @return:  {String}
        '''
        return os.path.basename(self.full)

    @property
    def prefix(self):
        '''
        File name. Without extension.

        @return:  {String}
        '''
        return os.path.splitext(self.name)[0]

    @property
    def first_prefix(self):
        '''
        File name. First section before a '.'.

        @return:  {String}
        '''
        return self.name.split('.')[0]

    @property
    def extension(self):
        '''
        File extension.

        @return:  {String}
        '''
        return os.path.splitext(self.name)[-1]

    @property
    def action(self):
        '''
        File selected action

        @return:  {String}
        '''
        return self._action

    @property
    def descriptor(self):
        msg = 'Use read() instead of descriptor'
        SBIg.warn(msg)
        return self._fd

    @property
    def size(self):
        '''
        File size.

        @return:  {String}
        '''
        return os.path.getsize(self.full)

    @property
    def split_section(self):
        '''
        Value when the file has been generated by spliting another file_name

        @return: {String}
        '''
        return self._section

    @split_section.setter
    def split_section(self, value):
        '''
        @param:    value
        @pdef:     new value for the attribute
        @ptype:    {String}
        '''
        self._section = value

    ############
    # BOOLEANS #
    ############
    @property
    def is_regular_file(self):
        '''
        @return:  {Boolean}
        '''
        return self._type == 'regular'

    @property
    def is_tarfile(self):
        '''
        @return:  {Boolean}
        '''
        return self._type.startswith('tar')

    @property
    def is_compressed(self):
        '''
        @return:  {Boolean}
        '''
        return self.is_gzipped or self.is_zipped or self.is_bzipped

    @property
    def is_gzipped(self):
        '''
        @return:  {Boolean}
        '''
        return self._type.endswith('gzip')

    @property
    def is_zipped(self):
        '''
        @return:  {Boolean}
        '''
        return self._type.endswith('zip') and not self._type.endswith('gzip')

    @property
    def is_bzipped(self):
        '''
        @return:  {Boolean}
        '''
        return self._type.endswith('bzip2')

    @property
    def has_multiple_files(self):
        '''
        @return:  {Boolean}
        '''
        return self._multi > 1

    ###########
    # METHODS #
    ###########
    def open(self):
        '''
        Upon starting a {File} object, it needs to be opened.

        @return: {file descriptor}
        '''
        if self._is_open:
            return

        self._is_open = True
        if self.is_regular_file:
            self._fd = open(self.full, self.action)
        if self.is_gzipped and not self.is_tarfile:
            self._fd = gzip.open(self.full, self.action)
        if self.is_bzipped and not self.is_tarfile:
            self._fd = bz2.BZ2File(self.full, self.action)
        if self.is_tarfile:  # Only to read
            self._fd = tarfile.open(self.full)
        if self.is_zipped:   # Only to read
            self._fd = zipfile.ZipFile(self.full)

    def close(self, clean = False):
        '''
        Closes the file once we are done with it

        @param:    clean
        @pdef:     erase the file on closing if file size is 0.
        @pdefault: _False_
        @ptype:    {Boolean}

        @warns: that files will be deleted if _clean_ is _True_
        '''
        self._fd.close()
        if clean and self.size == 0:   # Delete empty files
            w = '{0} is empty and it\'s going to be deleted'.format(self.full)
            SBIg.warn(w)
            os.unlink(self.full)
        self._is_open = False

    def read(self):
        '''
        Read from file.

        @raises: {FileError} if wrong action or multiple files
        @return: {file descriptor}
        '''
        if not self._action in File.READ_ACTION:
            raise self._error.wrong_action('read')
        if self.has_multiple_files:
            raise self._error.multiple()

        if self.is_regular_file or self.is_gzipped:
            return self._fd
        if self.is_bzipped:
            return self._fd.read()
        if self.is_zipped:
            return self._fd.read(self._fd.namelist()[0]).split('\n')
        if self.is_tarfile:
            return self._fd.extractfile(self._fd.getmembers()[0])

    def skip_lines(self, number):
        '''
        Activate to skip n lines before reading

        @param:    number
        @pdef:     lines to skip
        @ptype:    {Integer}
        '''
        for i in range(number):
            self._fd.readline()

    def read_file(self, file_name):
        '''
        Reads a specific file inside a multi-file container.

        @param: file_name
        @pdef:  name of the internal file to read.
        @ptype: {String}

        @raises: {FileError} if file_name not in container or
                 if not multi-file.
        @return: {file descriptor}
        '''
        if not self.has_multiple_files:
            raise self._error.no_multiple()
        if self.is_zipped:
            if file_name not in self._fd.namelist():
                raise self._error.no_file_in_file(file_name)
            return self._fd.read(file_name).split('\n')
        if self.is_tarfile:
            if file_name not in self._fd.getnames():
                raise self._error.no_file_in_file(file_name)
            return self._fd.extractfile(self._fd.getmember(file_name))

    def has_file(self, file_name):
        '''
        Checks if a specific file exists inside a multi-file container.

        @param: file_name
        @pdef:  name of the internal file to read.
        @ptype: {String}

        @raises: {FileError} if not multi-file.
        @return: {Boolean}
        '''
        if not self.has_multiple_files:
            raise self._error.no_multiple()
        if self.is_zipped:
            return file_name in self._fd.namelist()
        if self.is_tarfile:
            return file_name in self._fd.getnames()

    def write(self, line):
        '''
        Write into the file.

        @param: line
        @pdef:  content to write in the file.
        @ptype: {String}

        @raises: {FileError} if wrong action
        '''
        if not self._action in File.WRITE_ACTION:
            raise self._error.wrong_action('write')
        self._fd.write(line)

    def flush(self):
        '''
        Force to flush the contained data

        @raises: {FileError} if wrong action
        '''
        if not self._action in File.WRITE_ACTION:
            raise self._error.wrong_action('write')
        self._fd.flush()

    def list_files(self):
        '''
        List internal files inside a compressed or tar file.

        @raises: {FileError} if file does not contain multiple files.
        @return: {List} of file names
        '''
        if not self.has_multiple_files:
            raise self._error.no_multiple()
        else:
            if self.is_tarfile:
                return self._fd.getnames()
            if self.is_zipped:
                return self._fd.namelist()

    def split(self, size    = 1073741824, ini_separator        = None,
              end_separator = None, sections = None, overwrite = None,
              target_dir    = os.getcwd(), compress            = False):
        '''
        Splits a given file in multiple files.
        Splitting can be performed by file size, number of expected files
        and/or according to a given separator. Number of expected files takes
        precedence over size.

        A defined separator will avoid splitting data.
        For example, uniprot files separate between proteins with '//',
        then by declaring and end separator we can avoid a protein info from
        being split between two files.
        Similarly, fasta info starts with '>'. Thus, declaring a initial
        separator will ensure not to split a protein sequence.

        @param:    size
        @pdef:     Max size expected for a file. IN BYTES.
        @pdefault: 1073741824. (1GB).
        @ptype:    {Integer}

        @param:    ini_separator
        @pdef:     Separator definition.
        @pdefault: _None_
        @pclash:   end_separator
        @ptype:    {String}

        @param:    end_separator
        @pdef:     Separator definition.
        @pdefault: _None_
        @pclash:   ini_separator
        @ptype:    {String}

        @param:    sections
        @pdef:     Number of expected files.
        @pdefault: _None_
        @ptype:    {Integer}

        @param:    overwrite
        @pdef:     Overwrite files if they already exist.
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @param:    target_dir
        @pdef:     Directory to which the new files will be created.
        @pdefault: Current working directory.
        @ptype:    {String}

        @param:    compress
        @pdef:     New files will be gzipped.
        @pdefault: _False_
        @ptype:    {Boolean}

        @warns:  If self is a multi-file container, it will ignore parameters
                 and call self.extract.
        @raises: {FileError} if can not write in target directory
        @raises: {AttributeError} if both start_separator and end_separator
                 are defined OR if size AND sections are both not defined,
                 OR if size OR sections are not Integers.
        '''
        if self.has_multiple_files:
            SBIg.warn('For multi-file containers use File.extract()')
            self.extract(target_dir)
        if not os.path.isdir(target_dir) or not os.access(target_dir, os.W_OK):
            raise self._error.file_permissions('write', target_dir)
        if (size is None and sections is None):
            errmsg = 'Either size or sections must be defined.\n'
            raise AttributeError(errmsg)
        if (ini_separator is not None and end_separator is not None):
            errmsg = 'Both separators can not be defined simultaneously\n'
            raise AttributeError(errmsg)
        if size is not None:
            try:
                size = int(size)
            except:
                errmsg = 'Attribute size must be numeric.\n'
                raise AttributeError(errmsg)
        if sections is not None:
            try:
                sections = int(sections)
            except:
                errmsg = 'Attribute sections must be numeric.\n'
                raise AttributeError(errmsg)

        counter  = 1
        newfiles = []

        # section TAKES PRECEDENCE over size...
        if sections == 1:
            return [self, ]
        if sections is not None:
            size     = math.ceil(self.size / float(sections))
        else:
            sections = math.ceil(self.size / float(size))

        # preparing partition file names
        sectionid      = '{0:003d}'
        globalfilename = self.prefix + '.' + sectionid + self.extension
        outputfile     = os.path.join(target_dir, globalfilename)
        if compress:
            outputfile += '.gz'

        # Local overwrite takes precedence over Global overwrite
        overwrite = SBIg.decide_overwrite(overwrite)
        if not self._is_open:
            self.open()

        msg = 'Dividing {0.full} into {1} files (approx.)'.format(self, sections)
        SBIg.alert(0, msg)

        newfiles.append(File(file_name = outputfile.format(counter),
                             action    = 'w', overwrite = overwrite))
        newfiles[-1].split_section = sectionid.format(counter)
        for line in self.read():
            if ini_separator is not None and line.startswith(ini_separator):
                if newfiles[-1].size >= size:
                    newfiles[-1].close()
                    counter += 1
                    newfiles.append(File(file_name = outputfile.format(counter),
                                         action    = 'w', overwrite = overwrite))
                    newfiles[-1].split_section = sectionid.format(counter)

            newfiles[-1].write(line)

            if end_separator is not None and line.startswith(end_separator):
                if newfiles[-1].size >= size:
                    newfiles[-1].close()
                    counter += 1
                    newfiles.append(File(file_name = outputfile.format(counter),
                                         action    = 'w', overwrite = overwrite))
                    newfiles[-1].split_section = sectionid.format(counter)

            if end_separator is None and ini_separator is None:
                if newfiles[-1].size >= size:
                    newfiles[-1].close()
                    counter += 1
                    newfiles.append(File(file_name = outputfile.format(counter),
                                         action    = 'w', overwrite = overwrite))
                    newfiles[-1].split_section = sectionid.format(counter)

        newfiles[-1].close()
        self.close()
        self.open()
        return newfiles

    def extract(self, target_dir = os.getcwd()):
        '''
        Extract the files inside a multi-file container

        @param:    target_dir
        @pdef:     Directory to which the new files will be extracted.
        @pdefault: Current working directory.
        @ptype:    {String}

        @raises: {FileError} if can not write in target directory or if not
                 multi-file container.
        '''
        if not self.has_multiple_files:
            raise self._error.no_multiple()
        if not os.path.isdir(target_dir) or not os.access(target_dir, os.W_OK):
            raise self._error.file_permissions('write', target_dir)

        if self.is_tarfile:
            self._fd.extractall(path = target_dir)
        if self.is_zipped:
            self._fd.extractall(path = target_dir)

    ###################
    # PRIVATE METHODS #
    ###################
    def _check_file_type(self):
        '''
        Decide which kind of file we are working with.
        It uses the ACTIVE_EXTENSIONS to decide if some checking is necessary.
        '''
        if self._file.endswith(tuple(File.ACTIVE_EXTENSIONS)):
            if self._file.endswith(tuple(File.ZIP_EXTENSIONS)):
                self._type = 'zip'
                self._check_number_of_files()
            if self._file.endswith(tuple(File.TAR_EXTENSIONS)):
                self._type = 'tar'
                self._check_number_of_files()
            if self._file.endswith(tuple(File.GZIP_EXTENSIONS)):
                self._type = 'gzip' if self._type == '' else self._type + '.gzip'
            if self._file.endswith(tuple(File.BZIP_EXTENSIONS)):
                self._type = 'bzip2' if self._type == '' else self._type + '.bzip2'
        else:
            self._type = 'regular'

    def _check_number_of_files(self):
        '''
        In tar and zip files.
        Check the number of files contained inside.
        '''
        if self.is_tarfile:
            tar         = tarfile.open(self.full)
            self._multi = len(tar.getmembers())
        if self.is_zipped:
            zzip        = zipfile.ZipFile(self.full)
            self._multi = len(zzip.namelist())

    def _check_action(self, action):
        '''
        Sees that the requested action for the given file is valid and adapts
        it in case the file is compressed.

        @param: action
        @pdef:  type of action intended to apply to the File.
        @ptype: {String}

        @raises: {FileError}
        '''
        if not action in File.AVAILABLE_ACTION:  # Raise 'Wrong action' Error
            raise self._error.wrong_file_action(action)

        if self.is_tarfile and action.startswith('r'):
            action = 'r|*'

        if self.is_gzipped:
            if action.startswith('r'):
                action = 'rb'
            elif action.startswith('w'):
                action = 'wb'

        self._action = action

    def _check_file(self):
        '''
        Sees that the requested validated action can be performed over the
        selected file.

        @raises: {FileError}
        '''
        if self._action.startswith('r'):
            if os.path.isdir(self.full):
                raise self._error.is_dir()
            if not os.path.isfile(self.full):
                raise self._error.file_does_not_exist()
            if not os.access(self.full, os.R_OK):
                raise self._error.file_permissions('read')
        if self._action.startswith('w') or self._action.startswith('a'):
            if self.is_tarfile or self.is_zipped:
                raise self._error.unwritable_format()
            if os.path.isfile(self.full)and not self._overwrite:
                raise self._error.file_overwrite()
            if not os.path.isdir(self.dir):
                raise self._error.dir_does_not_exist()
            if not os.access(self.dir, os.W_OK):
                raise self._error.file_permissions('write')

    #################
    # CLASS METHODS #
    #################
    def __repr__(self):
        header = '[ {0.__class__.__name__} ][{1}]: '
        data   = '{0}{1.full}\n'.format(header.format(self,   'NAME  '), self)
        data  += '{0}{1._type}\n'.format(header.format(self,  'TYPE  '), self)
        data  += '{0}{1._multi}\n'.format(header.format(self, 'FILES '), self)
        data  += '{0}{1.action}\n'.format(header.format(self, 'ACTION'), self)
        return data

    def __str__(self):
        return self.full


class FileError(Exception):
    '''
    Manages different error produced by the {File} object when trying to
    read or write.

    '''
    def __init__(self, file_name):
        '''
        @param: file_name
        @pdef:  name of the file to which the error is related.
        @ptype: {String}
        '''
        self._file = file_name
        self._msg  = ''

    ###########
    # METHODS #
    ###########
    def wrong_file_action(self, action):
        msg  = '{0} is not an available action for file {1}.\n'
        msg += 'Accepted actions are {0}'.format(repr(File.AVAILABLE_ACTION))

        self._msg = msg.format(action, self._file)
        return self

    def file_does_not_exist(self):
        self._msg = 'File {0} does not exist.\n'.format(self._file)
        return self

    def file_permissions(self, what, where = None):
        where     = self._file if where is None else where
        self._msg = 'Permission denied to {0} {1}'.format(what, where)
        return self

    def file_overwrite(self):
        self._msg = 'Can not overwrite {0}'.format(self._file)
        return self

    def dir_does_not_exist(self):
        self._msg = 'Directory for {0} does not exist'.format(self._file)
        return self

    def is_dir(self):
        self._msg = 'Requested file {0} is a directory'.format(self._file)
        return self

    def is_file(self):
        self._msg = 'Requested directory {0} is a file'.format(self._file)
        return self

    def unwritable_format(self):
        self._msg = 'This format does not support write'.format(self._file)
        return self

    def wrong_action(self, what):
        msg  = 'Unable to {0}.\n'.format(what)
        msg += 'File {1} not in {0} mode'.format(what, self._file)

        self._msg = msg
        return self

    def no_multiple(self):
        self._msg = '{0} does not contain multiple files'.format(self._file)
        return self

    def multiple(self):
        self._msg = '{0} contain multiple files'.format(self._file)
        return self

    def no_file_in_file(self, search):
        self._msg = 'No {0} in {1}'.format(search, self._file)
        return self

    #################
    # CLASS METHODS #
    #################
    def __str__(self):
        return self._msg
