'''
@file: Executable.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   03/2013

@ [oliva's lab](http://sbi.imim.es)

@class: Executable
@class: ExecutableError
'''

import os
import subprocess
import copy

from .. import SBIglobals as SBIg


class Executable(object):
    '''
    Manages executing external programs through the SBI library.
    '''
    def __init__(self, executable, path=None):
        '''
        @param:    executable
        @pdef:     name of the executable program.
        @ptype:    {String}

        @param:    path
        @pdef:     path to the executable. Not needed if the executable is in
                   the $PATH environment variable.
        @pdefault: _None_
        @ptype:    {String}
        '''
        self._exec = executable
        self._path = path

        self._error = ExecutableError(self.executable)

        if path is not None:
            self._path = os.path.abspath(path)
        else:
            found = self._load_executable_path()
            if not found:
                raise self._error.not_in_path()

        self._error = ExecutableError(self.full_executable)

        self._check_executable()

        self._command = []
        self._command.append(self.full_executable)

        self._backup_command = []

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def executable(self):
        '''
        Name of the executable.

        @return: {String}
        '''
        return self._exec

    @property
    def path(self):
        '''
        Path to the executable.

        @return: {String}
        '''
        return self._path

    @property
    def command(self):
        '''
        Full command to execute.

        @return {List}
        '''
        return self._command

    @property
    def full_executable(self):
        '''
        Full executable with path.

        @return: {String}
        '''
        return os.path.join(self._path, self._exec)

    ###########
    # METHODS #
    ###########
    def add_attribute(self, attribute_value, attribute_id=None):
        '''
        Adds a new parameter to the command.
        Specifically for parameters with 'tags' like '-i'.

        @param:    attribute_value
        @pdef:     value of the attribute to add.
        @ptype:    {String}. Transformed otherwise.

        @param:    attribute_id
        @pdef:     label of the attribute to add.
        @pdefault: _None_
        @ptype:    {String}
        '''
        if attribute_id is None:
            self.add_parameter(attribute_value)
        else:
            self._command.append(attribute_id)
            self._command.append(str(attribute_value))

    def add_parameter(self, parameter):
        '''
        Adds a new stand alone parameter to the command.

        @param:    parameter
        @pdef:     value of the parameter to add.
        @ptype:    {String}. Transformed otherwise.
        '''
        self._command.append(str(parameter))

    def clean_command(self):
        '''
        Removes all attributes and parameters added to the command.
        '''
        self._command = []
        self._command.append(self.full_executable)

    def backup_command(self):
        '''
        Store a copy of the command up to that point to retrieve import
        afterwards.
        '''
        self._backup_command = copy.deepcopy(self._command)

    def restore_command(self):
        '''
        Retrieve the backup command into the working command.
        '''
        self._command        = self._backup_command
        self._backup_command = []

    def execute(self, silent=False):
        '''
        Executes the commands.

        @param:    silent
        @pdef:     External program STDERR is shown through STDERR if _True_
        @pdefault: _False_
        @ptype:    {Boolean}

        @raises: {SystemError} if an error occurs in the external program.
        '''
        outPIPE = subprocess.PIPE if SBIg.debug else open('/dev/null', 'w')
        errPIPE = open('/dev/null', 'w') if silent else subprocess.PIPE

        command = " ".join(self.command)
        SBIg.alert('debug', self, '\tExecuting:\t{0}'.format(command))

        try:
            p = subprocess.Popen(self.command, stdout = outPIPE,
                                               stderr = errPIPE)

            if SBIg.debug:
                for out in iter(p.stdout.readline, b''):
                    SBIg.alert('debug', self, out.strip())
            p.communicate()
        except:
            raise SystemError()

    ###################
    # PRIVATE METHODS #
    ###################
    def _load_executable_path(self):
        '''
        Retrieves the executable path in case self._path is None.

        @return: {Boolean}. _True_ if the path is found, _False_ otherwise.
        '''
        search = ["which", self.executable]
        p = subprocess.Popen(search, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()

        if out != '':
            self._path = os.path.split(out.strip())[0]
            return True
        else:
            return False

    def _check_executable(self):
        '''
        Checks that the final executable can be executed.

        @raises: {ExecutableError} otherwise.
        '''
        if not os.path.isfile(self.full_executable):
            raise self._error.no_exists()
        if not os.access(self.full_executable, os.X_OK):
            raise self._error.no_permission()

    #################
    # CLASS METHODS #
    #################
    def __repr__(self):
        return " ".join(self._command)

    def __str__(self):
        return repr(self)


class ExecutableError(Exception):
    '''
    Manages different error produced by the {Executable} object.

    '''
    def __init__(self, executable):
        '''
        @param: executable
        @pdef:  full path of the executable.
        @ptype: {String}
        '''
        self._exec = executable
        self._msg  = ''

    ###########
    # METHODS #
    ###########
    def not_in_path(self):
        self._msg = '{0} can not be found. Please specify.'.format(self._exec)
        return self

    def no_exists(self):
        self._msg = '{0} does not exist.'.format(self._exec)
        return self

    def no_permission(self):
        self._msg = 'No permission to execute {0}.'.format(self._exec)
        return self

    #################
    # CLASS METHODS #
    #################
    def __str__(self):
        return self._msg
