"""drugBanklink

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
from SBI.databases import drugBankftp
from SBI.beans     import Path
from SBI.beans     import File

class DrugBanklink(object):
    def __init__(self, local = None):
        self._local    = os.path.abspath(local)
        self._tagets   = None
        self._main     = None
        self._drugfile = None
        if local is not None:
            self.local = local

    """ATTRIBUTES"""
    @property
    def local(self):         return self._local
    @local.setter
    def local(self, value):
        self._local    = os.path.abspath(value)
        self._target   = os.path.join(self._local, os.path.split(drugBankftp['targets'])[-1])
        self._main     = os.path.join(self._local, os.path.split(drugBankftp['main'])[-1])
        self._drugfile = os.path.join(self._local, 'drugbank.gz')

    @property
    def source(self):
        return drugBankftp['show']

    """BOOLEANS"""
    @property
    def has_local(self):    return self._local is not None

    """METHODS"""
    def download(self):
        if not self.has_local:
            raise NameError('A local drugBank database directory must be defined.')

        Path.mkdir(self.local)
        urllib.urlretrieve(drugBankftp['targets'], self._target)
        urllib.urlretrieve(drugBankftp['main'],    self._main)

        self._process()

        return True

    @property
    def localDrugs(self):
        drgFile = File(self._drugfile, 'r')
        for drg_line in drgFile.descriptor:
            yield drg_line

    """PRIVATE METHODS"""
    def _process(self):

        targets = self._process_targets()
        drugs   = self._process_drugs(targets)

        drugFile = File(self._drugfile,'w', True)
        for d in drugs:
            drugFile.write(repr(d)+"\n")
        drugFile.close()

    def _process_targets(self):
        import csv, StringIO, zipfile

        targets = {}
        z       = zipfile.ZipFile(self._target, 'r')
        data    = StringIO.StringIO(z.read(z.namelist()[0]))
        reader  = csv.reader(data)
        for row in reader:
            if row[5] != '':
                targets[row[0]] = (row[5], row[6])
        z.close()
        return targets

    def _process_drugs(self, targets):
        import zipfile

        drugs = []
        z = zipfile.ZipFile(self._main, 'r')
        for line in z.open(z.namelist()[0]):
            if line.startswith('  <drug'):
                drugs.append(Drug())
                drug        = drugs[-1]
                drug.type   = line
            elif line.startswith('    <drugbank-id'):                   drug.dbid       = line
            elif line.startswith('    <name'):                          drug.name       = line
            elif line.startswith('    <description>'):                  drug.desc       = line
            elif line.startswith('    <cas-number'):                    drug.cas        = line
            elif line.startswith('      <secondary-accession-number'):  drug.secondary_accession_numbers = line
            elif line.startswith('      <group'):                       drug.groups     = line
            elif line.startswith('        <kind'):                      propkind        = line
            elif line.startswith('        <value'):                     drug.dproperty  = (propkind, line)
            elif line.startswith('      <category'):                    drug.categories = line
            elif line.startswith('        <resource>'):                 propkind        = line
            elif line.startswith('        <identifier>'):               drug.external   = (propkind, line)
            elif line.startswith('      <synonym>'):                    drug.synonym    = line
            elif 'partner=' in line:
                m = re.search('<(\w+) partner="(\d+)"', line)
                if m.group(2) in targets:
                    drug.partners = [m.group(1), targets[m.group(2)]]
                else:
                    drug.partners = [m.group(1), ('','')]
                partner = drug.partners[-1]
            elif line.startswith('          <action>'):
                partner.append(re.search(">([\s\S]+)</", line).group(1))
            elif line.startswith('        <known-action>'):
                partner.append(re.search(">(\S+)</", line).group(1).lower())
            # elif line.startswith('  </drug'):
            #     # print line.strip()
            #     print drug
        z.close()
        return drugs

class Drug(object):

    properties_of_interest = ["Hydrophobicity", "Isoelectric Point", "Molecular Weight", "Molecular Formula"]

    def __init__(self, inline = None):
        if inline is not None: inline = inline.strip().split('\t')

        self._type = None if inline is None else inline[0]
        self._dbid = None if inline is None else inline[1]
        self._name = None if inline is None else inline[2]
        self._cas  = None if inline is None else inline[3]
        self._desc = None if inline is None else inline[4]
        self._san  = []   if inline is None else eval(inline[5])
        self._grps = []   if inline is None else eval(inline[6])
        self._prop = {}   if inline is None else eval(inline[7])
        self._cate = []   if inline is None else eval(inline[8])
        self._extr = []   if inline is None else eval(inline[9])
        self._synm = []   if inline is None else eval(inline[10])
        self._part = []   if inline is None else eval(inline[11])

    @property
    def type(self): return self._type
    @type.setter
    def type(self, value): self._type = re.search('type\=\"(\w+\s*\w+)\"', value).group(1)

    @property
    def dbid(self): return self._dbid
    @dbid.setter
    def dbid(self, value): self._dbid = re.search(">(DB\d+)</", value).group(1)

    @property
    def name(self): return self._name
    @name.setter
    def name(self, value): self._name = re.search(">([\S\s]+)</", value).group(1)

    @property
    def cas(self): return self._cas
    @cas.setter
    def cas(self, value):
        m = re.search(">([\d\-]+)</", value)
        self._cas = m.group(1) if m else ''

    @property
    def desc(self): return self._desc
    @desc.setter
    def desc(self, value):
        m = re.search(">([\S\s]+)\s*[<|&]", value)
        self._desc = m.group(1).rstrip('< ') if m else ''

    @property
    def secondary_accession_numbers(self): return self._san
    @secondary_accession_numbers.setter
    def secondary_accession_numbers(self, value): self._san.append(re.search(">([\S\s]+)</", value).group(1))

    @property
    def groups(self): return self._grps
    @groups.setter
    def groups(self, value): self._grps.append(re.search(">([\S\s]+)</", value).group(1))

    @property
    def dproperty(self): return self._prop
    @dproperty.setter
    def dproperty(self, value):
        if re.search(">([\S\s]+)</", value[0]).group(1) in self.properties_of_interest:
            self._prop[re.search(">([\S\s]+)</", value[0]).group(1)] = re.search(">([\S\s]+)</", value[1]).group(1).split()[0]

    @property
    def categories(self): return self._cate
    @categories.setter
    def categories(self, value): self._cate.append(re.search(">([\S\s]+)</", value).group(1))

    @property
    def external(self): return self._extr
    @external.setter
    def external(self, value):
        self._extr.append(tuple([re.search(">([\S\s]+)</", value[0]).group(1), re.search(">([\S\s]+)</", value[1]).group(1)]))

    @property
    def synonym(self): return self._synm
    @synonym.setter
    def synonym(self, value): self._synm.append(re.search(">([\S\s]+)</", value).group(1))

    @property
    def partners(self): return self._part
    @partners.setter
    def partners(self, value): self._part.append(value)

    def __repr__(self):
        data = []
        data.append(self._type)
        data.append(self._dbid)
        data.append(self._name)
        data.append(self._cas)
        data.append(self._desc)
        data.append(repr(self._san))
        data.append(repr(self._grps))
        data.append(repr(self._prop))
        data.append(repr(self._cate))
        data.append(repr(self._extr))
        data.append(repr(self._synm))
        data.append(repr(self._part))

        return "\t".join(data)
