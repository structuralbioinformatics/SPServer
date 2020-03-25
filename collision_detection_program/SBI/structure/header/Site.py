import re
from .MiniRes import MiniResidue
from SBI.beans.JSONer import JSONer

class Site(JSONer):

    def __init__(self, identifier):
        self._id    = identifier
        self._desc  = ''
        self._evid  = ''
        self._spec  = []
        self._bind  = []

    #
    # ATTRIBUTES
    #
    @property
    def identifier(self):
        return self._id

    @property
    def description(self):
        return self._desc

    @property
    def evidence(self):
        return self._evid

    @property
    def residues(self):
        return self._spec

    @property
    def binding(self):
        return self._bind

    @property
    def chains(self):
        data = []
        if len(self.residues) != 0:
            data.extend(self.residues)
        if len(self.binding) != 0:
            data.extend(self.binding)
        if len(data) > 0:
            data = list(set([x.chain for x in data]))
        return data

    #
    # READ FUNCTIONS
    #
    def add_remark(self, line):
        if line.startswith('EVIDENCE_CODE:'):
            self._evid = line.split(':')[1].strip()
        elif line.startswith('SITE ') or line == 'SITE':
            pass
        elif line.startswith('SITE_DESCRIPTION:'):
            self._add_description(line.split(':')[1].strip())
        elif len(line) > 0:
            self._add_description(line.strip())

    def add_residues(self, residues):
        for r in residues:
            if r[0] != 'HOH':  # SKIP WATERS
                self._spec.append(MiniResidue(r[0], r[1]))

    def _add_description(self, line):
        self._desc += line

        m = re.search('BINDING SITE FOR RESIDUE\s+(\w+)', self._desc)
        if m:
            res = m.group(1)
            pos = re.sub('BINDING SITE FOR RESIDUE\s+\w+', '', line)
            self._bind.append(MiniResidue(res, pos.strip().strip(',')))

    #
    # OUT FUNCTION
    #
    def as_dict(self):
        nobj = {}
        nobj['id']          = self.identifier
        nobj['description'] = self.description
        nobj['evidence']    = self.evidence

        if len(self.binding) > 0:
            nobj['binds'] = [x.as_dict() for x in self.binding]
        else:
            nobj['binds'] = self.binding

        if len(self.residues) > 0:
            nobj['residues'] = [x.as_dict() for x in self.residues]
        else:
            nobj['residues'] = []

        return nobj

    def __repr__(self):
        return self.json()
