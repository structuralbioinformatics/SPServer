"""PDBTMlink

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
from SBI.databases import PDBTMftp
from SBI.beans    import Path
from SBI.beans     import File

class PDBTMlink(object):

    def __init__(self, local = None):
        self._local     = os.path.abspath(local)
        self._main      = None
        self._pdbtmfile = None
        if local is not None:
            self.local = local

    """ATTRIBUTES"""
    @property
    def local(self):         return self._local
    @local.setter
    def local(self, value):
        self._local     = os.path.abspath(value)
        self._pdbtmfile = os.path.join(self._local, 'pdbtm.gz')

    @property
    def source(self):
        return PDBTMftp['show']

    """BOOLEANS"""
    @property
    def has_local(self):    return self._local is not None

    """METHODS"""
    def download(self):
        if not self.has_local:
            raise NameError('A local drugBank database directory must be defined.')

        Path.mkdir(self.local)
        here = os.getcwd()
        os.chdir(self.local)
        os.system("svn export {0}".format(PDBTMftp['svn']))

        self._process()

        return True

    @property
    def localTM(self):
        tmoFile = File(self._pdbtmfile, 'r')
        for tm_line in tmoFile.descriptor:
            yield tm_line

    """PRIVATE METHODS"""
    def _process(self):
        tmoFile = File(self._pdbtmfile,'w', True)
        for xmlfile in Path.list_files(os.path.join(self._local,'pdbtm/database/'), '*.xml'):
            xmldata = TM(pdb = os.path.splitext(os.path.split(xmlfile)[1])[0].upper())
            skip_chains = set()
            read = False
            fdxml = open(xmlfile)
            for line in fdxml:
                if line.startswith('    <TMRES>'):     xmldata.tmres  = line
                elif line.startswith('    <TMTYPE'):   xmldata.tmtype = line
                elif line.startswith('    <PDBKWRES'): xmldata.kwres  = line
                elif line.startswith('  <SIDEDEFINITION'):
                    m = re.search('Side1="(\S+)"', line)
                    xmldata.side = m.group(1)
                elif line.startswith('      <APPLY_TO_CHAIN'):
                    m = re.search('NEW_CHAINID=\"(\S{1})\"', line)
                    if m: skip_chains.add(m.group(1))
                elif line.startswith('  <CHAIN '):
                    m = re.search('CHAINID=\"(\S{1})\" NUM_TM=\"(\d{1})\" TYPE=\"(\S+)\"', line)
                    if m:
                        chain, num, tmtype = m.group(1), m.group(2), m.group(3)
                        if not chain in skip_chains:
                            cdata = tuple([chain, num, tmtype])
                            xmldata.set_chain(cdata)
                            read = True
                elif line.startswith('    <REGION ') and read:
                    m = re.search('pdb_beg=\"(\-*\d+\w*)\"[\s\S]+pdb_end=\"(\-*\d+\w*)\"\s+type=\"(\w{1})\"', line)
                    ini, end, tmtype = m.group(1), m.group(2), m.group(3)
                    xmldata.set_chain(cdata, tuple([ini, end, tmtype]))
                elif line.startswith('  </CHAIN>'): read = False
            fdxml.close()
            if len(xmldata.chains) > 0:
                tmoFile.write(str(xmldata)+"\n")
        tmoFile.close()


class TM(object):
    section_types = {'1':'Side1', '2':'Side2', 'i':'Inside', 'o':'Outside',
                     'B':'Beta-strand', 'H':'alpha-helix', 'C':'coil',
                     'I':'membrane-inside', 'L':'membrane-loop', 'F':'interfacial helix', 'U':'unknown localizations'}

    def __init__(self, pdb= None, inline = None):
        if inline is not None: inline = inline.strip().split('\t')

        self._pdb    = pdb  if inline is None else eval(inline[0])
        self._tmres  = ''   if inline is None else eval(inline[1])
        self._tmtype = ''   if inline is None else eval(inline[2])
        self._kwres  = ''   if inline is None else eval(inline[3])
        self._side   = ''   if inline is None else eval(inline[4])
        self._chains = {}   if inline is None else eval(inline[5])

    @property
    def pdb(self): return self._pdb
    @property
    def tmres(self): return self._tmres
    @tmres.setter
    def tmres(self, value): self._tmres = TM._clean_xml(value)
    @property
    def tmtype(self): return self._tmtype
    @tmtype.setter
    def tmtype(self, value): self._tmtype = TM._clean_xml(value)
    @property
    def kwres(self): return self._kwres
    @kwres.setter
    def kwres(self, value): self._kwres = TM._clean_xml(value)
    @property
    def side(self): return self._side
    @side.setter
    def side(self, value): self._side = value
    @property
    def chains(self): return self._chains

    def set_chain(self, key, value = None):
        if value is None:
            self._chains.setdefault(key,[])
        else:
            self._chains[key].append(value)

    @staticmethod
    def _clean_xml(xml):
        return re.search(">([\S\s]+)</", xml).group(1)

    def __str__(self):
        data = []
        data.extend([self._pdb,self._tmres,self._tmtype,self._kwres, self._side])
        data.append(self._chains)
        return "\t".join([repr(x) for x in data])


