'''
@file: ExternalExe.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2014

@ [oliva's lab](http://sbi.imim.es)

@class: ExternalExe
'''
import os
import ConfigParser
from os.path import join     as j
from os.path import normpath as n
from os.path import dirname  as d
from abc     import ABCMeta

from SBI.beans import Executable


class ExternalExe(object):
    '''
    Main class to derive others that will control the execution of external
    programs.

    Program configuration is defined in the configSBI.txt file or can be
    defined in a file linked to a environment variable called SBI_CONFIG_FILE.

    '''
    __metaclass__ = ABCMeta

    DEFAULT_CONFIG_FILE = j(n(d(__file__)), 'configSBI.txt')

    # Executable configuration
    _CONFIG     = ConfigParser.RawConfigParser(allow_no_value=True)
    _EXE        = None

    def __new__(cls, *args, **kwargs):
        cls._CONFIG.read(os.getenv('SBI_CONFIG_FILE', cls.DEFAULT_CONFIG_FILE))
        return super(ExternalExe, cls).__new__(cls, *args, **kwargs)

    ###################
    # PRIVATE METHODS #
    ###################
    def _set_default_executable(self, external_id):
        '''
        Configure the {Executable} of the class according to the default
        configuration as defined in the default configSBI.txt or the
        SBI_CONFIG_FILE

        @param:    external_id
        @pdef:     name of the external program as referred in the
                   configSBI.txt file.
        @ptype:    {String}
        '''
        i, e, p   = external_id, 'executable', 'path'
        self._EXE = Executable(executable= self._CONFIG.get(i, e),
                               path      = self._CONFIG.get(i, p))

    @staticmethod
    def _set_dynamic_executable(executable, path):
        '''
        Manual configuration of the {Executable}.

        @param:    executable
        @pdef:     name of the executable file
        @ptype:    {String}

        @param:    path
        @pdef:     path to the executable file
        @ptype:    {String}
        '''
        ExternalExe._EXE = Executable(executable = executable, path = path)
