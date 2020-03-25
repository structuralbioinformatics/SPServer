"""Enzymelink

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
from SBI.databases import Enzymeftp
from SBI.beans     import Path
from SBI.beans     import File

class Enzymelink(object):
    """The Enzymelink class controls the download and parsing of Enzyme database

    """
    def __init__(self, local = None):
        self._local   = os.path.abspath(local)
        self._dfile   = 'enzyme.dat'
        self._cfile   = 'enzclass.txt'
        self.__name__ = 'databases.Enzymelink'    # This must be included in every class for the SBIglobals.alert()
        self._enzfile = 'enzyme.gz'
        if local is not None:
            self.local = local

    """ATTRIBUTES"""
    @property
    def local(self):         return self._local
    @local.setter
    def local(self, value):
        self._local   = os.path.abspath(value)
        self._dfile   = os.path.join(self._local, self._dfile)
        self._cfile   = os.path.join(self._local, self._cfile)
        self._enzfile = os.path.join(self._local, self._enzfile)

    @property
    def localEnzymes(self):
        enzFile = File(self._enzfile, 'r')
        for enz_line in enzFile.descriptor:
            yield enz_line

    @property
    def source(self):
        return Enzymeftp['show']

    """BOOLEANS"""
    @property
    def has_local(self):    return self._local is not None

    """METHODS"""
    def download(self):
        if not self.has_local:
            raise NameError('A local Enzyme database directory must be defined.')

        Path.mkdir(self.local)
        urllib.urlretrieve(Enzymeftp['dat'], self._dfile)
        urllib.urlretrieve(Enzymeftp['cls'], self._cfile)

        self._process()

        return True

    def get_Enzyme(self, enzyme):
        if self.has_local:
            for enz_line in self.localEnzymes:
                if enz_line.split('\t')[0] == enzyme:
                    return enz_line
        else:
            raise NameError('A local Enzyme database directory must be defined.')

    def get_Enzymes(self, ENZset):
        if isintance(ENZset, str):
            warnings.warn('For single enzyme search the get_Enzyme function is recomended.')
            yield self.get_Enzyme(ENZset)

        if self.has_local:
            for enz_line in self.localEnzymes:
                if enz_line.split('\t')[0] in ENZset:
                    yield enz_line
        else:
            raise NameError('A local Enzyme database directory must be defined.')

    """PRIVATE METHODS"""
    def _process(self):
        enzymes = self._parse_enzclass() + self._parse_enzymedat()
        enzymes.sort()

        enzFile = File(self._enzfile,'w', True)
        for e in enzymes:
            enzFile.write(repr(e)+"\n")
        enzFile.close()

    def _parse_enzclass(self):
        enzymes = []
        fd = open(self._cfile, 'r')
        for line in fd:
            if re.match(r'^\d', line):
                line = re.sub(r'\. ','.', line)
                line = re.sub('\'', '\\\'', line)
                data = line.split()
                enzymes.append(Enzyme(ec = data[0], description = ' '.join(data[1:])))
        fd.close()
        return enzymes

    def _parse_enzymedat(self):
        code          = None
        name          = None
        cofr          = ''
        reac          = ''
        reactions_dic = {}
        enzymes       = []

        cofactor_separator = re.compile(r"; | or |;| and ")
        multiple_reac_sepr = re.compile(r"\.\(\d+\)\s")
        noisy_words_equal  = re.compile(r"^A |^An |secondary | a | an ")
        unspaced_summator  = re.compile(r" \+|\+ ")
        unspaced_equalate  = re.compile(r" \=|\= ")
        summator_separator = re.compile(r" \+ | or ")
        equalate_separator = re.compile(r" \= ")
        cardinal_numerical = re.compile(r"^(\d+)\s+")

        fd = open(self._dfile)
        for line in fd:
            line = re.sub('\'', '\\\'', line)
            line = re.sub('<element>','', line)
            line = re.sub('</element>','', line)
            if line.startswith('//'):
                if code is not None:
                    enzymes.append(Enzyme(code, name))
                    if cofr is not '':
                        cofr = cofr[:-1]
                        for fct in cofactor_separator.split(cofr):
                            enzymes[-1].cofactors = fct

                    if reac is not '':
                        if re.search(r"\=", reac):
                            for rct in multiple_reac_sepr.split(reac):
                                repeated = False
                                if rct.startswith('(1)'): rct = rct[4:]
                                if rct.endswith('.'):     rct = rct[:-1]
                                rct = re.sub(noisy_words_equal, '', rct)
                                rct = re.sub(unspaced_summator, ' + ', rct)
                                rct = re.sub(unspaced_equalate, ' = ', rct)
                                rct = re.sub('\s+',' ', rct)
                                rct = re.sub('\> \=', '\>\=', rct)
                                rct = re.sub('\< \=', '\<\=', rct)
                                rct = re.sub('\(n \= ', '(n=', rct)

                                if not rct in reactions_dic:
                                    rid = len(reactions_dic) + 1
                                    reactions_dic[rct] = rid
                                else:
                                    rid = reactions_dic[rct]
                                    repeated = True

                                enzymes[-1].reactions.append(Reaction(description = rct, repeated = repeated))
                                if repeated: continue
                                try:
                                    (substrate, product) = equalate_separator.split(rct)
                                except:
                                    (substrate, product) = unspaced_equalate.split(rct)
                                cardinality = 1
                                for subs in summator_separator.split(substrate):
                                    cardinal = cardinal_numerical.search(subs)
                                    if cardinal:
                                        cardinality = cardinal.group(1).strip()
                                        subs = re.sub(cardinal_numerical, '', subs)
                                        subs.strip()
                                    enzymes[-1].reactions[-1].substrates.append((subs,cardinality))
                                    cardinality = 1
                                for prod in summator_separator.split(product):
                                    cardinal = cardinal_numerical.search(prod)
                                    if cardinal:
                                        cardinality = cardinal.group(1).strip()
                                        prod = re.sub(cardinal_numerical, '', prod)
                                        prod.strip()
                                    enzymes[-1].reactions[-1].products.append((prod,cardinality))
                                    cardinality = 1

                        else:
                            if not reac in reactions_dic:
                                rid = len(reactions_dic) + 1
                                reactions_dic[reac] = rid
                            else:
                                rid = reactions_dic[reac]
                            enzymes[-1].reactions.append(Reaction(description = reac))
                    if len(unip) > 0:
                        enzymes[-1].uniprots = unip

                code = None
                name = None
                cofr = ''
                reac = ''
                unip = []
                continue
            if line.startswith('ID'):
                code = line.split()[-1].strip()
                name = ""
                continue
            if line.startswith('DE'):
                if name != "":
                    name += " "
                name += " ".join(line.split()[1:]).strip()
                # name = name[:-1]
                continue
            if line.startswith('CF'):
                cofactors = " ".join(line.split()[1:]).strip()
                cofr = cofr + cofactors
                continue
            if line.startswith('CA'):
                reactions = " ".join(line.split()[1:]).strip()
                reac = reac + reactions
                continue
            if line.startswith('DR'):
                unipcodes = [x.replace(' ','').split(',')[1] for x in line[2:].strip().split(';')[:-1]]
                unip.extend(unipcodes)
                continue
        fd.close()
        return enzymes

class Enzyme(object):

    def __init__(self, ec = None, description = None, inline = None):
        if inline is not None:
            inline = inline.strip().split('\t')
        self.ec          = ec.replace(' ','')     if inline is None else inline[0]
        self.description = description.strip()    if inline is None else inline[1]
        self.deleted     = self.description == 'Deleted entry.'
        self.level       = 4 - self.ec.count('-') if inline is None else int(inline[2])
        self.parents     = []                     if inline is None else eval(inline[3])
        self.direct_p    = 'NULL'                 if inline is None else inline[4]
        self._cofactors  = []                     if inline is None else eval(inline[5])
        self.transfers   = []                     if inline is None else eval(inline[6])
        self.uniprots    = []                     if inline is None else eval(inline[7])
        self.reactions   = []

        if inline is None:
            self._parents()
            self._transfers()
        else:
            self._reactions(eval(inline[8]))

    @property
    def cofactors(self):
        return self._cofactors
    @cofactors.setter
    def cofactors(self, value):
        if re.search('\(',value) and not re.search('\)',value): value = re.sub('\(','',value)
        if re.search('\)',value) and not re.search('\(',value): value = re.sub('\)','',value)
        self._cofactors.append(value.upper())

    @property
    def compounds(self):
        for reaction in self.reactions:
            for subs in reaction.substrates:
                yield subs[0]
            for prod in reaction.products:
                yield prod[0]

    @property
    def has_transfers(self):
        return len(self.transfers) != 0
    @property
    def has_direct_parent(self):
        return self.direct_p != 'NULL'
    @property
    def is_deleted(self):
        return self.deleted

    def _reactions(self, reactarray):
        for react in reactarray:
            self.reactions.append(Reaction(inline = react))

    def _parents(self):
        self.parents.append(self.ec)
        pieces = self.ec.split('.')
        for i in range(3,0,-1):
            if pieces[i] != '-':
                pieces[i] = '-'
                parent = '.'.join(pieces)
                self.parents.append(parent)
                pieces = parent.split('.')
        try:
            self.direct_p = self.parents[1]
        except:
            self.direct_p = 'NULL'
        self.parents.reverse()

    def _transfers(self):
        if self.description.startswith('Transferred entry'):
            d = self.description
            d.strip('.')
            aggregator_separator = re.compile(r",|and")
            self.transfers = [x.strip().rstrip('.') for x in aggregator_separator.split(d.split(':')[1])]

    def __lt__(self, other):
        me = [int(x) for x in self.ec.replace('-','0').replace('n','1000').split('.')]
        he = [int(x) for x in other.ec.replace('-','0').replace('n','1000').split('.')]

        for i in range(4):
            if me[i] > he[i]:
                return False
            elif me[i] < he[i]:
                return True
        return False

    def __gt__(self, other):
        me = [int(x) for x in self.ec.replace('-','0').replace('n','1000').split('.')]
        he = [int(x) for x in other.ec.replace('-','0').replace('n','1000').split('.')]

        for i in range(4):
            if me[i] > he[i]:
                return True
            elif me[i] < he[i]:
                return False
        return False

    def __repr__(self):
        data = []
        data.append('{0.ec}\t{0.description}\t{0.level}\t{0.parents}\t{0.direct_p}'.format(self))
        data.append('{0.cofactors}\t{0.transfers}\t{0.uniprots}'.format(self))
        data.append('{0.reactions}'.format(self))

        return "\t".join(data)
        # return '{0.ec}\t{0.description}\t{0.level}\t'.format(self) + repr(self.parents) + repr(self.reactions)

class Reaction(object):

    def __init__(self, description = None, repeated = False, inline = None):
        self.description = description.rstrip('\\') if inline is None else inline[0]
        self.substrates  = []                       if inline is None else inline[1]
        self.products    = []                       if inline is None else inline[2]
        self.repeated    = repeated                 if inline is None else inline[3]

    def __repr__(self):
        data = []
        data.append(self.description)
        data.append(self.substrates)
        data.append(self.products)
        data.append(self.repeated)
        return '{0}'.format(data)
        # return '{0.description},{0.substrates},{0.products},{0.repeated}'.format(self)

