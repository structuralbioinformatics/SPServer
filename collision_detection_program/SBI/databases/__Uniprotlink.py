"""Uniprotlink

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
from SBI.databases import Uniprotftp
from SBI.beans.Path     import Path
from SBI.beans.File     import File
from SBI.beans.StorableObject     import StorableObject


class Uniprotlink(object):
    """The Uniprotlink class controls the download and parsing of Uniprot database

    """
    def __init__(self, local = None):
        self._local   = os.path.abspath(local)
        self.__name__ = 'databases.Uniprotlink'    # This must be included in every class for the SBIglobals.alert()
        self._swsfile = 'swisprot.gz'
        self._swsfast = 'swisprot.fasta'
        self._trbfile = 'trembl.gz'
        self._trbfast = 'trembl.fasta'

        if local is not None:
            self.local = local

    """ATTRIBUTES"""
    @property
    def local(self):         return self._local
    @local.setter
    def local(self, value):
        self._local   = os.path.abspath(value)
        self._swsfile = os.path.join(self._local, self._swsfile)
        self._swsfast = os.path.join(self._local, self._swsfast)
        self._trbfile = os.path.join(self._local, self._trbfile)
        self._trbfast = os.path.join(self._local, self._trbfast)

    @property
    def localUniprots(self):
        for uni_line in self.localSwissprots:
            yield uni_line
        for uni_line in self.localTrembls:
            yield uni_line

    @property
    def localSwissprots(self):
        swsFile = File(self._swsfile, 'r')
        for uni_line in swsFile.descriptor:
            yield uni_line

    @property
    def localTrembls(self):
        tblFile = File(self._trbfile, 'r')
        for uni_line in tblFile.descriptor:
            yield uni_line

    @property
    def source(self):
        return Uniprotftp['show']

    """BOOLEANS"""
    @property
    def has_local(self):    return self._local is not None

    """METHODS"""
    def download(self):
        if not self.has_local:
            raise NameError('A local Uniprot database directory must be defined.')

        Path.mkdir(self.local)
        destination = os.path.join(self.local, 'uniprot_sprot.dat.gz')
        urllib.urlretrieve(Uniprotftp['swissprot'], destination)
        destination = os.path.join(self.local, 'uniprot_trembl.dat.gz')
        urllib.urlretrieve(Uniprotftp['trembl'], destination)

        self._process()

        return True

    def get_Uniprot(self, Uniprot):
        if self.has_local:
            for uni_line in self.localUniprots:
                if uni_line.split('\t')[0] == Uniprot:
                    return uni_obj
        else:
            raise NameError('A local Uniprot database directory must be defined.')

    def get_Uniprots(self, UNIset):
        if isintance(UNIset, str):
            warnings.warn('For single Uniprot search the get_Uniprot function is recomended.')
            yield self.get_Uniprot(UNIset)

        if self.has_local:
            for uni_line in self.localUniprots:
                if uni_line.split('\t')[0] in UNIset:
                    yield uni_line
        else:
            raise NameError('A local Uniprot database directory must be defined.')

    """PRIVATE METHODS"""
    def _process(self):
        self._parse_uniprot_file(source = os.path.join(self.local, 'uniprot_sprot.dat.gz'),
                                 destination = self._swsfile, fasta = self._swsfast, code = 'swissprot')
        self._parse_uniprot_file(source = os.path.join(self.local, 'uniprot_trembl.dat.gz'),
                                 destination = self._trbfile, fasta = self._trbfast, code = 'trembl')

    def _parse_uniprot_file(self, source, destination, fasta, code):
        sourceFile = File(source, 'r')
        destinFile = File(destination, 'w', True)
        fastaFile  = File(fasta, 'w', True)
        protein = None
        for line in sourceFile.descriptor:
            if line.startswith('ID'):
                protein = Uniprot(line.split()[1].strip(), code)
            if line.startswith('AC'):
                protein.accession = line.split()[1:]
            if line.startswith('OX'):
                protein.taxid     = line.split()[1]
            if line.startswith('OH'):
                protein.hosts     = line.split()[1]
            if line.startswith('DR'):
                protein.databases = line.split()[1:3]
            if line.startswith('  '):
                protein.sequence  = line.strip().replace(' ','')
            if line.startswith('//'):
                destinFile.write(str(protein) + "\n")
                fastaFile.write(repr(protein) + "\n")
        sourceFile.close()
        destinFile.close()


class Uniprot(object):
    def __init__(self, entry = None, source = None, inline = None):
        if inline is not None:
            inline      = inline.strip()
        self.entry      = entry  if inline is None else inline.split('\t')[0]
        self.source     = source if inline is None else inline.split('\t')[1]
        self._accession = []     if inline is None else eval(inline.split('\t')[2])
        self._ox        = None   if inline is None else inline.split('\t')[3] if inline.split('\t')[3] != 'None' else None
        self._oc        = []     if inline is None else eval(inline.split('\t')[4])
        self._db        = []     if inline is None else eval(inline.split('\t')[5])
        self._seq       = None   if inline is None else inline.split('\t')[6] if inline.split('\t')[6] != 'None' else None

    @property
    def accession(self): return self._accession
    @accession.setter
    def accession(self, value):
        self._accession.extend([x.strip(';') for x in value])

    @property
    def taxid(self): return self._ox
    @taxid.setter
    def taxid(self, value):
        if self._ox is None:
            self._ox = self._extract_taxid(value)
        else:
            raise('More than one taxid! {0}'.format(value))

    @property
    def hosts(self): return self._oc
    @hosts.setter
    def hosts(self, value):
        self._oc.append(self._extract_taxid(value))

    @property
    def sequence(self): return self._seq
    @sequence.setter
    def sequence(self, value):
        if self._seq is None:
            self._seq = value.strip()
        else:
            self._seq += value.strip()

    @property
    def databases(self): return self._db
    @databases.setter
    def databases(self, value):
        self._db.append((value[0].strip(';'),value[1].strip(';')))

    def _extract_taxid(self, value):
        import re
        taxid_id = re.compile('NCBI_TaxID\=(\d+);')
        m        = taxid_id.search(value)
        return m.group(1)

    """OVERWRITE PARENT METHODS"""
    def __repr__(self):
        return ">{0.entry}\n{0.sequence}".format(self)

    def __str__(self):
        return "{0.entry}\t{0.source}\t{0.accession}\t{0.taxid}\t{0.hosts}\t{0.databases}\t{0.sequence}".format(self)
