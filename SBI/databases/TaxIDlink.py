"""TaxIDlink

author: jbonet
date:   10/2013

@oliva's lab
"""

"""
Import Standard Libraries
"""
import sys, os, re
import subprocess
import warnings
import urllib

"""
Dependences in SBI library
"""
from SBI               import SBIglobals
from SBI.databases import taxIDftp
from SBI.beans     import Path
from SBI.beans     import File

class TaxIDlink(object):
    """The TaxIDlink class controls the download and parsing of TaxID database

    """
    def __init__(self, local = None):
        self._local   = os.path.abspath(local)
        self.__name__ = 'databases.TaxIDlink'    # This must be included in every class for the SBIglobals.alert()

        self._nodes  = 'nodes.dmp'
        self._names  = 'names.dmp'
        self._delet  = 'delnodes.dmp'
        self._merged = 'merged.dmp'
        self._taxid  = 'taxid.gz'

        if local is not None:
            self.local = local

    """ATTRIBUTES"""
    @property
    def local(self):         return self._local
    @local.setter
    def local(self, value):
        self._local  = os.path.abspath(value)
        self._nodes  = os.path.join(self._local, self._nodes)
        self._names  = os.path.join(self._local, self._names)
        self._delet  = os.path.join(self._local, self._delet)
        self._merged = os.path.join(self._local, self._merged)
        self._taxid  = os.path.join(self._local, self._taxid)

    @property
    def localTaxIDs(self):
        taxFile = File(self._taxid, 'r')
        for tax_line in taxFile.descriptor:
            yield tax_line
        taxFile.close()

    @property
    def source(self):
        return taxIDftp['show']

    """BOOLEANS"""
    @property
    def has_local(self):    return self._local is not None

    """METHODS"""
    def download(self):
        if not self.has_local:
            raise NameError('A local TaxID database directory must be defined.')

        Path.mkdir(self.local)
        destination = os.path.join(self.local, 'taxdmp.zip')
        urllib.urlretrieve(taxIDftp['global'], destination)
        command = ['unzip', '-o', destination, '-d', self.local]
        p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = p.communicate()

        self._process()

        return True

    def get_TaxID(self, TaxID):
        if self.has_local:
            for tax_line in self.localTaxIDs:
                if tax_line.split('\t')[0] == TaxID:
                    return tax_line

        else:
            raise NameError('A local TaxID database directory must be defined.')

    def get_TaxIDs(self, TAXset):
        if isintance(TAXset, str):
            warnings.warn('For single TaxID search the get_TaxID function is recomended.')
            yield self.get_TaxID(TAXset)

        if self.has_local:
            for tax_line in self.localTaxIDs:
                if tax_line.split('\t')[0] in TAXset:
                    yield tax_linej
        else:
            raise NameError('A local TaxID database directory must be defined.')

    """PRIVATE METHODS"""
    def _process(self):
        inh = {}
        nodefile = File(file_name = self._nodes, action = 'r')
        for line in nodefile.descriptor:
            line = re.sub('\'', '\\\'', line)
            line_data = line.split('|')
            inh[line_data[0].strip()] = TaxID(line_data[0].strip())
            inh[line_data[0].strip()].parent = line_data[1].strip()
            inh[line_data[0].strip()].rank   = line_data[2].strip()
        nodefile.close()

        namefile = File(file_name = self._names, action = 'r')
        for line in namefile.descriptor:
            line = re.sub('\'', '\\\'', line)
            line_data = line.split('|')
            if line_data[3].strip() == 'scientific name':
                inh[line_data[0].strip()].name = line_data[1].strip()
        namefile.close()

        delefile = File(file_name = self._delet, action = 'r')
        for line in delefile.descriptor:
            data = line.split('|')
            inh[data[0].strip()]     = TaxID(data[0].strip())
            inh[data[0].strip()].old = True
        delefile.close()

        mrgefile = File(file_name = self._merged, action = 'r')
        for line in mrgefile.descriptor:
            data = line.split('|')
            inh[data[0].strip()]     = TaxID(data[0].strip())
            inh[data[0].strip()].old = True
            inh[data[0].strip()].new = data[1].strip()
        mrgefile.close()

        taxFile = File(self._taxid, 'w', True)
        for taxid in inh:
            taxFile.write(str(inh[taxid]) + "\n")
        taxFile.close()

class TaxID(object):
    def __init__(self, taxid = None, inline = None):
        if inline is not None:
            inline      = inline.strip()
        self.taxid  = taxid if inline is None else inline.split('\t')[0]
        self.name   = None  if inline is None else inline.split('\t')[1] if inline.split('\t')[1] != 'None' else None
        self.rank   = None  if inline is None else inline.split('\t')[2] if inline.split('\t')[2] != 'None' else None
        self.parent = None  if inline is None else inline.split('\t')[3] if inline.split('\t')[3] != 'None' else None
        self.old    = False if inline is None else eval(inline.split('\t')[4])
        self.new    = None  if inline is None else inline.split('\t')[5] if inline.split('\t')[5] != 'None' else None

    """BOOLEANS"""
    @property
    def has_old(self): return self.old
    @property
    def has_new(self): return False if self.new is None else True

    """OVERWRITE PARENT METHODS"""
    def __str__(self):
        return "{0.taxid}\t{0.name}\t{0.rank}\t{0.parent}\t{0.old}\t{0.new}".format(self)
