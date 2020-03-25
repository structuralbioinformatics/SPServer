"""SCOPlink

author: jbonet
date:   02/2014

@oliva's lab
"""

"""
Import Standard Libraries
"""
import os, re
import urllib

"""
Dependences in SBI library
"""
from SBI.databases import SCOPftp
from SBI.beans     import Path
from SBI.beans     import File

class SCOPlink(object):
    def __init__(self, local = None):
        self._local     = os.path.abspath(local)
        self._desc      = None
        self._rel       = None
        if local is not None:
            self.local = local

    """ATTRIBUTES"""
    @property
    def local(self):         return self._local
    @local.setter
    def local(self, value):
        self._local    = os.path.abspath(value)
        self._desc     = os.path.join(self._local, os.path.split(SCOPftp['desc'])[-1])
        self._rel      = os.path.join(self._local, os.path.split(SCOPftp['rel'])[-1])

    @property
    def source(self):
        return SCOPftp['show']

    @property
    def descriptions(self):
        dscFile = File(self._desc, 'r')
        for dsc_line in dscFile.descriptor:
            if not dsc_line.startswith('#'):
                yield dsc_line

    @property
    def relations(self):
        relFile = File(self._rel, 'r')
        for rel_line in relFile.descriptor:
            if not rel_line.startswith('#'):
                yield rel_line

    """BOOLEANS"""
    @property
    def has_local(self):    return self._local is not None

    """METHODS"""
    def download(self):
        if not self.has_local:
            raise NameError('A local SCOP database directory must be defined.')

        Path.mkdir(self.local)
        urllib.urlretrieve(SCOPftp['desc'], self._desc)
        urllib.urlretrieve(SCOPftp['rel'],  self._rel)

        return True
