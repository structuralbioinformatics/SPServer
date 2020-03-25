"""GOlink

author: jbonet
date:   10/2013

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
from SBI.databases import GOftp, GOnamespace
from SBI.beans     import Path
from SBI.beans     import File

class GOlink(object):
    """The GOlink class controls the download and parsing of GO database

    """
    def __init__(self, local = None):
        self._local   = os.path.abspath(local)
        self._gfile   = 'go-simple.obo'
        self.__name__ = 'databases.GOlink'    # This must be included in every class for the SBIglobals.alert()
        self._gofile  = 'go.gz'
        if local is not None:
            self.local = local

    """ATTRIBUTES"""
    @property
    def local(self):         return self._local
    @local.setter
    def local(self, value):
        self._local  = os.path.abspath(value)
        self._gofile = os.path.join(self._local, self._gofile)

    @property
    def localGOs(self):
        goFile = File(self._gofile, 'r')
        for go_line in goFile.descriptor:
            yield go_line

    @property
    def source(self):
        return GOftp['show']

    """BOOLEANS"""
    @property
    def has_local(self):    return self._local is not None

    """METHODS"""
    def download(self):
        if not self.has_local:
            raise NameError('A local GO database directory must be defined.')

        Path.mkdir(self.local)
        destination = os.path.join(self.local, self._gfile)

        urllib.urlretrieve(GOftp['source'], destination)
        self._process()

        return True

    def get_GO(self, GO):
        if self.has_local:
            for go_line in self.localGOs:
                if go_line.split('\t')[0] == GO.upper():
                    return go_line
        else:
            raise NameError('A local GO database directory must be defined.')

    def get_GOs(self, GOset):
        if isintance(GOset, str):
            warnings.warn('For single GO search the get_GO function is recomended.')
            yield self.get_GO(GOset)
        else:
            GOset = set([x.upper() for x in GOset])

        if self.has_local:
            for go_line in self.localGOs:
                if go_line.split('\t')[0] in GOset:
                    yield go_line
        else:
            raise NameError('A local GO database directory must be defined.')

    """PRIVATE METHODS"""
    def _process(self):
        go_dic = {}
        parseFile = File(os.path.join(self.local, self._gfile), 'r')
        go = None
        for line in parseFile.descriptor:
            line = re.sub('\'', '\\\'', line)
            if line.startswith('[Term]'):
                if go is not None:
                    go_dic[go.id] = go
            if line.startswith('id:'):
                go = GOterm(id = line.split()[1].strip())
                continue
            if line.startswith('name:'):
                go.name = " ".join(line.split()[1:]).strip()
                continue
            if line.startswith('namespace:'):
                go.namespace = line.split()[1].strip()
                continue
            if line.startswith('alt_id:'):
                go.alt_id.append(line.split()[1].strip())
                continue
            if line.startswith('is_obsolete:'):
                go.obsolete = True
                continue
            if line.startswith('is_a:'):
                go.parents.add(line.split()[1].strip())
                continue
            if line.startswith('relationship:'):
                go.relations.append((line.split()[1].strip(),line.split()[2].strip()))
                continue
            if line.startswith('[Typedef]'):
                go_dic[go.id] = go
                break
        parseFile.close()

        for go in go_dic:
            go_dic[go].parents = self._search_parents(go_dic, go)

        goFile = File(self._gofile, 'w', True)
        for go in go_dic:
            go_dic[go].parents.add(go)
            goFile.write(str(go_dic[go]) + "\n")
        goFile.close()

    def _search_parents(self, godict, go):
        parents = set()
        for parent in godict[go].parents:
            parents.add(parent)
            parents.update(self._search_parents(godict, parent))
        return parents

class GOterm(object):

    def __init__(self, id = None, inline = None):
        if inline is not None:
            inline      = inline.strip()
        self.id         = id    if inline is None else inline.split('\t')[0]
        self.name       = None  if inline is None else inline.split('\t')[1]
        self._namespace = None  if inline is None else inline.split('\t')[2]
        self.obsolete   = False if inline is None else False if int(inline.split('\t')[3]) == 0 else True
        self.alt_id     = []    if inline is None else eval(inline.split('\t')[4])
        self.parents    = set() if inline is None else eval(inline.split('\t')[5])
        self.relations  = []    if inline is None else eval(inline.split('\t')[6])

    """ATTRIBUTES"""
    @property
    def namespace(self):        return self._namespace
    @namespace.setter
    def namespace(self, value): self._namespace = value
    @property
    def namespace_mini(self):   return GOnamespace[self._namespace]
    @property
    def obsolete_binary(self):  return 1 if self.obsolete else 0

    """OVERWRITE PARENT METHODS"""
    def __str__(self):
        return "{0.id}\t{0.name}\t{0.namespace}\t{0.obsolete_binary}\t{0.alt_id}\t{0.parents}\t{0.relations}".format(self)

