"""PDBeChemlink

author: jbonet
date:   10/2013

@oliva's lab
"""

"""
Import Standard Libraries
"""
import sys
import os
import re
import subprocess
import warnings
import urllib

"""
Dependences in SBI library
"""
from SBI           import SBIglobals
from SBI.databases import PDBeChemftp
from SBI.beans     import Path
from SBI.beans     import File
from SBI.data      import element_dic


class PDBeChemlink(object):

    """The PDBeChemlink class controls the download and parsing of PDBeChem database

    """

    def __init__(self, local = None):
        self._local   = os.path.abspath(local)
        self.__name__ = 'databases.PDBeChemlink'    # This must be included in every class for the SBIglobals.alert()

    """ATTRIBUTES"""
    @property
    def local(self): return self._local

    @local.setter
    def local(self, value): self._local = os.path.abspath(value)

    @property
    def localPDBeChems(self):
        for chem_file in Path.list_files(root = self.local, pattern = '*.cif'):
            yield chem_file

    @property
    def source(self):
        return PDBeChemftp['show']

    """BOOLEANS"""
    @property
    def has_local(self): return self._local is not None
    """METHODS"""

    def download(self):
        if not self.has_local:
            raise NameError('A local PDBeChem database directory must be defined.')

        Path.mkdir(self.local)
        destination = os.path.join(self.local, 'mmcif.tar.gz')
        try:
            urllib.urlretrieve(PDBeChemftp['global'], destination)
        except:
            return False
        command = ['tar', 'zxvf', destination, '-C', self.local]
        p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = p.communicate()

        return True

    def get_PDBeChem(self, chemID):
        if self.has_local:
            for chem_file in self.localPDBeChems:
                newfile = File(file_name = chem_file, action = 'r')
                if newfile.prefix.upper() == chemID.upper():
                    return chem_file

        # If we do not find it in local (or we do not have a local) we search it on the FTP
        chem_file = chemID.upper() + '.cif'
        source = PDBeChemftp['single'] + chem_file
        try:
            urllib.urlretrieve(source, chem_file)
        except:
            return False
        return os.path.abspath(chem_file)

    def get_PDBeChems(self, chemIDset):
        if isintance(chemIDset, str):
            warnings.warn('For single PDBeChem search the get_PDBeChem function is recomended.')
            yield self.get_PDBeChem(chemIDset)
        else:
            chemIDset = set([x.upper() for x in chemIDset])

        if self.has_local:
            for chem_file in self.localPDBeChems:
                newfile = File(file_name = chem_file, action = 'r')
                if newfile.prefix.lstrip('pdb').upper() in chemIDset:
                    yield chem_file
        else:
            for chemID in chemIDset:
                yield self.get_PDBeChem(chemID)


class PDBeChem(object):

    """
    """

    def __init__(self, cif_file):
        self._file    = File(file_name = cif_file, action = 'r')
        self.__name__ = 'databases.PDBeChem'    # This must be included in every class for the SBIglobals.alert()

        self._id        = None
        self._name      = None
        self._type      = None
        self._formula   = None
        self._parent    = None
        self._weight    = None
        self._fcharge   = None
        self._code1l    = None
        self._flformula = {}

        self._parse()
        self._decompose_formula()

    """ATTRIBUTES"""
    @property
    def id(self): return self._id

    @property
    def name(self): return self._name

    @property
    def type(self): return self._type

    @property
    def formula(self): return self._formula

    @property
    def full_formula(self): return self._flformula

    @property
    def parent(self): return self._parent

    @property
    def weight(self): return self._weight

    @property
    def formal_charge(self): return self._fcharge

    @property
    def code1(self): return self._code1l

    @property
    def code3(self): return self._id
    """PRIVATE METHODS"""

    def _parse(self):
        for line in self._file.descriptor:
            if line.startswith('_chem_comp.'):
                line  = line.replace('_chem_comp.', '')
                value = line[35:].strip().strip('"')
                value = value.replace(' (NON-PREFERRED NAME)', '')
                value = value if value != '?' else None
                if line.startswith('id'): self._id      = value
                if line.startswith('pdbx_type'): self._type    = value
                if line.startswith('formula '): self._formula = value
                if line.startswith('formula_weight'): self._weight  = value
                if line.startswith('pdbx_formal_charge'): self._fcharge = value
                if line.startswith('one_letter_code'): self._code1l  = value
                if line.startswith('name'): self._name    = value.upper()
                if line.startswith('mon_nstd_parent_comp_id'):
                    self._parent  = set([x.strip() for x in value.split(',')]) if value is not None else None
            if line.startswith(';') and self._name == '': self._name  += line.strip().lstrip(';').upper()
        self._file.close()

    def _decompose_formula(self):
        if self.formula is not None:
            data    = self.formula.split()
            atregex = re.compile('(\D+)(\d*)')
            for atom in data:
                m = atregex.search(atom)
                if m.group(1) in element_dic:
                    self._flformula[m.group(1)] = m.group(2) if m.group(2) != '' else 1

    """OVERWRITE INHERITED FUNCTIONS"""

    def __str__(self):
        if self.code1 is not None and self.parent is not None:
            return "[{0.id} - {0.code1} from {0.parent}: {0.weight} - {0.formula} - {0.formal_charge}] {0.name} - {0.type}".format(self)
        elif self.code1 is not None:
            return "[{0.id} - {0.code1}: {0.weight} - {0.formula} - {0.formal_charge}] {0.name} - {0.type}".format(self)
        elif self.parent is not None:
            return "[{0.id} from {0.parent}: {0.weight} - {0.formula} - {0.formal_charge}] {0.name} - {0.type}".format(self)
        else:
            return "[{0.id}: {0.weight} - {0.formula} - {0.formal_charge}] {0.name} - {0.type}".format(self)
