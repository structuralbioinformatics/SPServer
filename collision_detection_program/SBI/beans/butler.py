'''
\usepackage{hyperref}

@file: butler.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2014

@ [oliva's lab](http://sbi.imim.es)

@class: Singleton
@class: Butler

\section{Description}
The \href{}{Singleton} object ensures that only one instance of a given object
is created. This is required as only one \href{}{Buttler} can go on around the
library consistently.

The \href{}{Buttler} class will create an instance that will move around the
library managing different aspects of it, such as overwriting permissions on
files or the logging.
'''
import sys
import os
import traceback
import datetime
import time
import logging
import inspect

from .singleton import Singleton


class Butler(object):
    '''
    Designed to work through all the \textbf{SBI library}.

    It contains several parameters that will control:
        \begin{itemize}
            \item the amount of data shown to the user during the execution of
                  \textbf{SBI} subroutines
            \item the file overwrite settings in a execution as long as it is
                  managed through a \href{}{File} object
        \end{itemize}

    Regarding the looging, it provides functions similar to those of a logger.
    '''
    __metaclass__   = Singleton

    _LOGNAME        = 'SBILOG'
    _GENERAL_FORMAT = '[{0}] ' + '%(asctime)s - %(levelname)-7.7s - %(message)s'
    _TIME_FORMAT    = '%Y-%m-%d %H:%M'

    def __init__(self):
        self._report = [
            False,  # Verbose
            False,  # Debug
            False,  # Deep Debug
            True,   # Warning
            True    # Error
        ]

        self._ios = [
            False,  # Overwrite
            False   # STDOUT
        ]

        self._fd = logging.getLogger(self._LOGNAME)
        self._fd.setLevel(logging.NOTSET)
        self._fd.addHandler(logging.NullHandler())

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
        return self._report[0]

    @verbose.setter
    def verbose(self, value):
        '''
        @param: {Boolean}
        '''
        self._report[0] = bool(value)

    @property
    def debug(self):
        '''
        Medium level of progress info.
        Use to better understand the process behind a particular function.

        @return:  {Boolean}
        @default: _False_
        '''
        return self._report[1]

    @debug.setter
    def debug(self, value):
        '''
        @param: {Boolean}
        '''
        value = bool(value)
        if value:
            self._report[0:2] = [value]*2
        else:
            self._debug       = value

    @property
    def deepdebug(self):
        '''
        Maximum level of progress info.
        To check on processes that are going to be highly repeated.

        @return:  {Boolean}
        @default: _False_
        '''
        return self._report[2]

    @deepdebug.setter
    def deepdebug(self, value):
        '''
        If deepdebug is set to _True_, it also sets _debug_ and _verbose_
        to _True_.

        @param: {Boolean}
        '''
        value = bool(value)
        if value:
            self._report[0:3] = [value]*3
        else:
            self._report[2]   = value

    @property
    def warnings(self):
        '''
        Prompts warning info in places where something is not completely wrong
        but there is no need to stop the execution.

        @return:  {Boolean}
        @default: _True_
        '''
        return self._report[3]

    @warnings.setter
    def warnings(self, value):
        '''
        @param: {Boolean}
        '''
        self._report[3] = bool(value)

    @property
    def errors(self):
        '''
        Prompts error info in places where something goes wrong.

        @return:  {Boolean}
        @default: _True_
        '''
        return self._report[4]

    @errors.setter
    def errors(self, value):
        '''
        @param: {Boolean}
        '''
        self._report[4] = bool(value)

    @property
    def overwrite(self):
        '''
        General overwrite setting for all the execution.
        It is OVER-RULED by local overwrite

        @return:  {Boolean}
        @default: _False_
        '''
        return self._ios[0]

    @overwrite.setter
    def overwrite(self, value):
        '''
        @param: {Boolean}
        '''
        self._ios[0] = bool(value)

    @property
    def stdout(self):
        '''
        STDOUT output status

        @return:  {Boolean}
        @default: _False_
        '''
        return self._ios[1]

    @stdout.setter
    def stdout(self, value):
        '''
        @param: {Boolean}
        '''
        self._ios[1] = bool(value)
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
        frmt = logging.Formatter(self._GENERAL_FORMAT.format(self._LOGNAME),
                                 self._TIME_FORMAT)
        if log_file is None:
            handler = logging.StreamHandler()
        else:
            handler = logging.FileHandler(log_file)
        handler.setFormatter(frmt)
        self._fd.addHandler(handler)
        if log_file is not None:
            self._fd.info('Logfile {0} created.'.format(log_file))

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

    def alert(self, level = 'verbose', message = None):
        '''
        Throws a message for the user.

        @param:    level
        @pdef:     specifies the minimum active level to actually show the
                   message.
        @pdefault: verbose
        @poptions: 'verbose', 'debug', 'deepdebug' or, equivalently,
                   0, 1, 2
        @ptype:    {String} or {Integer}

        param:     message
        @pdef:     specific message sent by the user.
        @pdefault: {None}
        @ptype:    {String}
        '''
        level, stat = self._active_level(level)
        if not stat:
            return

        callerID = inspect.getmodule(inspect.stack()[1][0]).__name__
        if callerID is '__main__':
            callerID = inspect.getmodule(inspect.stack()[1][0]).__file__
            callerID = os.path.splitext(callerID)[-1]
        callerID = '[' + callerID.upper() + '] '

        if message is None:
            message = ''

        if isinstance(message, list):
            for line in message:
                self.alert(level, line)
        else:
            if level == 'verbose':
                self._fd.info('{0}{1}'.format(callerID, message))
            else:
                self._fd.debug('{0}{1}'.format(callerID, message))

    def warn(self, message = None):
        '''
        Throw a warning message. Call when something should be said to the user
        about what's going on.

        @param:    message
        @pdef:     specific message sent by the user.
        @pdefault: {None}
        @ptype:    {String}

        '''
        if self.warnings is False:
            return

        callerID = inspect.getmodule(inspect.stack()[1][0]).__name__
        if callerID is '__main__':
            callerID = inspect.getmodule(inspect.stack()[1][0]).__file__
            callerID = os.path.splitext(callerID)[-1]
        callerID = '[' + callerID.upper() + '] '

        if isinstance(message, list):
            for line in message:
                self.warn(line)
        else:
            self._fd.warning('{0}{1}'.format(callerID, message))

    def throw(self, message = None, error = None, killit  = True):
        '''
        Throw an error. Call it when an error occur.
        If the _error_ attribute is _False_, it skips the effect, unless
        _killit_ is _True_, which will imply that we can not proceed given the
        error

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

        callerID = inspect.getmodule(inspect.stack()[1][0]).__name__
        if callerID is '__main__':
            callerID = inspect.getmodule(inspect.stack()[1][0]).__file__
            callerID = os.path.splitext(callerID)[-1]
        callerID = '[' + callerID.upper() + '] '

        self._fd.error('{0}{1}'.format(callerID, message))

        if error is not None:
            self._fd.error(error.message)
            self._fd.error(traceback.format_exc())

        if killit:
            sys.exit(-9)

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
        try:
            level = int(level)
            stat  = ((level == 0 and self.verbose) or
                    (level  == 1 and self.debug)   or
                    (level  == 2 and self.deepdebug))
            return ['verbose', 'debug', 'deepdebug'][level], stat
        except Exception:
            stat = ((level == 'verbose'   and self.verbose) or
                   (level  == 'debug'     and self.debug)   or
                   (level  == 'deepdebug' and self.deepdebug))
            return level, stat
