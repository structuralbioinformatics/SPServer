'''
@file: StorableObject.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   01/2013

@ [oliva's lab](http://sbi.imim.es)

@class: StorableObject
'''
from abc import ABCMeta
try:
    import cPickle as pickle
except:
    import pickle

from .  import File
from .. import SBIglobals as SBIg


class StorableObject(object):
    '''
    An abstract "dumping" class.

    This means that it is basically useful for those who would like to extend
    this library.
    Basically, it gives the object the ability to be "dumped" on disk and be
    recovered afterwards.

    '''
    __metaclass__ = ABCMeta

    def dump(self, object_file, overwrite = None):
        '''
        Stores the object into a file.

        @param:    object_file
        @pdef:     Name for the output file
        @ptype:    {String}

        @param:    overwrite
        @pdef:     write over a existing file
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}
        '''
        SBIg.alert('verbose', self, 'Writting object to file {0}'.format(object_file))
        dumpFile = File(file_name = object_file, action = 'w', overwrite = overwrite)
        pickle.dump(self, dumpFile._fd)
        dumpFile.close()

    @staticmethod
    def load(object_file):
        '''
        Retrieves the object from a python object file

        @param:    object_file
        @pdef:     Name of the file containing the object
        @ptype:    {String}

        @return: child of {StorableObject}
        '''
        SBIg.alert('verbose', StorableObject(), 'Preparing to load object from file {0}'.format(object_file))
        Object   = None
        loadFile = File(file_name = object_file, action='r')
        Object   = pickle.load(loadFile._fd)
        loadFile.close()
        return Object
