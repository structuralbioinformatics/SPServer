'''
@file: __init__.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2013

@ [oliva's lab](http://sbi.imim.es)

@class: Parameters
'''
import sys
import os
import traceback
import datetime
import time
import logging

__version__ = '0.0.1'


class Parameters(object):
    '''
    Designed to work through all the {SBI} library.

    It contains several parameters that will control (a) the amount of data
    shown to the user during the execution of {SBI} subroutines and (b) the
    file overwrite settings in a execution. See {SBI.beans.File} for more info
    in that topic.

    '''
    _LOGNAME = 'SBILOG'

    def __init__(self):
        self._verbose  = False       # > Minimum level of progress info.
                                     # Mainly, to inform of advance.
        self._debug    = False       # > Medium level of progress info.
                                     # Use to better understand the process
                                     # behind a particular function.
        self._ddebug   = False       # > Maximum level of progress info.
                                     # To check on processes that are going
                                     # to be highly repeated.
        self._overwrte = False       # > This will set the general overwrite
                                     # setting for all the execution.
                                     # It is OVER-RULED by local overwrite
                                     # parameters.
        self._warning  = True        # > Prompt warning info.
        self._error    = True        # > Prompt error info.

        self._stdout   = False       # STDOUT status

        self._fd       = logging.getLogger()  # > All messages can be redirected to
                                              # some file-handle different than STDERR
        self._fd.setLevel(logging.NOTSET)

    ##############
    # ATTRIBUTES #
    ##############

    @property
    def verbose(self):
        '''
        Minimum level of progress info.
        Mainly, to inform of advance.

        @return:  {Boolean}
        @default: _False_
        '''
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        '''
        @param: {Boolean}
        '''
        self._verbose = value

    @property
    def debug(self):
        '''
        Medium level of progress info.
        Use to better understand the process behind a particular function.

        @return:  {Boolean}
        @default: _False_
        '''
        return self._debug

    @debug.setter
    def debug(self, value):
        '''
        @param: {Boolean}
        '''
        if value:
            self._verbose, self._debug = value, value
        else:
            self._debug                = value

    @property
    def deepdebug(self):
        '''
        Maximum level of progress info.
        To check on processes that are going to be highly repeated.

        @return:  {Boolean}
        @default: _False_
        '''
        return self._ddebug

    @deepdebug.setter
    def deepdebug(self, value):
        '''
        If deepdebug is set to _True_, it also sets _debug_ and _verbose_
        to _True_.

        @param: {Boolean}
        '''
        if value:
            self._ddebug = value
            self.debug   = value
        else:
            self._ddebug = value

    @property
    def warnings(self):
        '''
        Prompts warning info in places where something is not completely wrong
        but there is no need to stop the execution.

        @return:  {Boolean}
        @default: _True_
        '''
        return self._warning

    @warnings.setter
    def warnings(self, value):
        '''
        @param: {Boolean}
        '''
        self._warning = value

    @property
    def errors(self):
        '''
        Prompts error info in places where something goes wrong.

        @return:  {Boolean}
        @default: _True_
        '''
        return self._error

    @errors.setter
    def errors(self, value):
        '''
        @param: {Boolean}
        '''
        self._error = value

    @property
    def overwrite(self):
        '''
        General overwrite setting for all the execution.
        It is OVER-RULED by local overwrite

        @return:  {Boolean}
        @default: _False_
        '''
        return self._overwrte

    @overwrite.setter
    def overwrite(self, value):
        '''
        @param: {Boolean}
        '''
        self._overwrte = value

    @property
    def stdout(self):
        '''
        STDOUT output status

        @return:  {Boolean}
        @default: _False_
        '''
        return self._overwrte

    @stdout.setter
    def stdout(self, value):
        '''
        @param: {Boolean}
        '''
        self._stdout = value
        if value:
            self.log_file()

    ###########
    # METHODS #
    ###########
    def log_file(self, log_file = None):
        '''
        Manually select a file to print the output information.

        The resultant log-file is independent from the overwrite parameter.
        It opens as any regular python file (that is, it does overwrite).

        @param:    log_file
        @pdef:     Name of the log file
        @pdefault: _None_, defaults to STDERR
        @ptype:    {String}
        '''
        frmt = logging.Formatter('[{0}] '.format(self._LOGNAME) +
                                 '%(asctime)s - %(levelname)-7.7s - %(message)s',
                                 '%Y-%m-%d %H:%M')
        if log_file is None:
            handler = logging.StreamHandler()
        else:
            handler = logging.FileHandler(log_file)
        handler.setFormatter(frmt)
        self._fd.addHandler(handler)
        if log_file is not None:
            self._fd.info('[SBIglobals]: Logfile {0} created.'.format(log_file))

    def set_file(self, work_dir = os.getcwd()):
        '''
        Automatically create a file to print the output information.

        The log-file consists in the name of the current executing process,
        the process-identifier (pid) number and the .log extension.

        The resultant log-file is independent from the overwrite parameter.
        It opens as any regular python file (that is, it does overwrite).

        @param:    work_dir
        @pdef:     Name of the directory into which create the file.
        @pdefault: Current working directory.
        @ptype:    {String}
        '''
        script_name = os.path.split(os.path.splitext(sys.argv[0])[0])[1]
        log_file    = ".".join([script_name, str(os.getpid()), 'log'])
        self.log_file(os.path.join(work_dir, log_file))

    def decide_overwrite(self, local_overwrite):
        '''
        When local file overwrite has not been set (is None), the *Parameters*
        _overwrite_ attribute is called. Otherwise local overwrite
        specification has priority.

        @param: local_overwrite
        @pdef:  Selected status to control the overwrite of an already existing
                file.
        @ptype: {Boolean}

        @return: {Boolean}
        '''
        return self.overwrite if local_overwrite is None else local_overwrite

    def alert(self, level = 'verbose', source_object = None, message = None):
        '''
        Throws a message for the user.

        @param:    level
        @pdef:     specifies the minimum active level to actually show the
                   message.
        @pdefault: verbose
        @poptions: 'verbose', 'debug', 'deepdebug'
        @ptype:    {String}

        @param:    source_object
        @pdef:     object that is currently throwing a message.
        @pdefault: {None}
        @ptype:    {object}

        param:     message
        @pdef:     specific message sent by the user.
        @pdefault: {None}
        @ptype:    {String}
        '''
        if not self._active_level(level):
            return
        if source_object is None and message is None:
            return

        name = self._source_name(source_object)
        if message is None:
            message = ''
        if isinstance(message, list):
            for line in message:
                self.alert(level, source_object, line)
        else:
            if level == 'verbose':
                self._fd.info('{0}{1}'.format(name, message))
            else:
                self._fd.debug('{0}{1}'.format(name, message))

    def warn(self, source_object = None, message = None):
        '''
        Throw a warning message. Call when something should be said to the user
        about what's going on.

        @param:    source_object
        @pdef:     object that is currently throwing a warning.
        @pdefault: {None}
        @ptype:    {object}

        @param:    message
        @pdef:     specific message sent by the user.
        @pdefault: {None}
        @ptype:    {String}

        '''
        if self.warnings is False:
            return

        name = self._source_name(source_object)
        if isinstance(message, list):
            for line in message:
                self.warn(source_object, line)
        else:
            self._fd.warning('{0}{1}'.format(name, message))

    def throw(self, source_object = None, message = None,
                    error         = None, killit  = True):
        '''
        Throw an error. Call it when an error occur.
        If the _error_ attribute is _False_, it skips the effect, unless
        _killit_ is _True_, which will imply that we can not proceed given the
        error

        @param:    source_object
        @pdef:     object that is currently reporting an error.
        @pdefault: {None}
        @ptype:    {object}

        @param:    message
        @pdef:     specific message sent by the user.
        @pdefault: {None}
        @ptype:    {String}

        @param:    error
        @pdef:     {Exception} that has been captured with the error.
        @pdefault: {None}
        @ptype:    {Exception}

        @param:    killit
        @pdef:     order to stop the execution after the error has been
                   reported.
        @pdefault: _True_
        @ptype:    {Boolean}
        '''
        if self.errors is False and killit is False:
            return

        name = self._source_name(source_object)

        self._fd.error('{0}{1}'.format(name, message))

        if error is not None:
            self._fd.error(error.message)
            self._fd.error(traceback.format_exc())

        if killit:
            sys.exit()

    def success(self):
        '''
        Tells that the program has successfully finished
        '''
        self._fd.info('[SUCCESS!!]: -- Program ended as expected.')
        sys.exit(0)

    def countdown(self, max_time):
        '''
        Prints a countdown in place.
        Put it in a loop if you are waiting for something.

        @param:    max_time
        @pdef:     time to wait. in seconds.
        @ptype:    {integer}
        '''
        t  = str(datetime.timedelta(seconds=max_time))
        n  = time.localtime()
        s1 = 'Waiting for: {0} hours'.format(t)
        s2 = 'Wait started at {0}'.format(time.strftime('%X', n))
        s3 = 'on {0}'.format(time.strftime('%Y-%m-%d', n))
        sys.stderr.write('{0}\t{1} {2}\n\n'.format(s1, s2, s3))
        while max_time > 0:
            t = str(datetime.timedelta(seconds=max_time))
            sys.stderr.write('Remaining: {0} hours'.format(t))
            time.sleep(1)
            max_time -= 1
            if bool(max_time):
                sys.stderr.write('\r')
            else:
                sys.stderr.write('\r')
                t = str(datetime.timedelta(seconds=max_time))
                sys.stderr.write('Remaining: {0} hours'.format(t))
                sys.stderr.write('\n')

    ###################
    # PRIVATE METHODS #
    ###################
    def _active_level(self, level):
        '''
        Decide if a given level is active or not.

        @param:    level
        @pdef:     level to check if it is active.
        @pdefault: {None}
        @ptype:    {String}

        @return: {Boolean}
        '''
        return ((level == 'verbose'   and self.verbose) or
               (level  == 'debug'     and self.debug)   or
               (level  == 'deepdebug' and self.deepdebug))

    def _source_name(self, source_object):
        '''
        Format the name of the object calling SBIglobals

        @param:    source_object
        @pdef:     object
        @ptype:    {object}
        '''
        if isinstance(source_object, basestring):
            return '[' + source_object.upper() + ']: '
        elif source_object is not None:
            return '[' + source_object.__class__.__name__.upper() + ']: '
        else:
            return ''

SBIglobals = Parameters()
'''
The {SBIglobals} variable is used through all the library, being the one that
actually sets the {Parameters} around the different classes and functions
of {SBI}.
'''
