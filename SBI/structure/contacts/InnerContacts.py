"""
InnerContacts

author: jbonet
date:   02/2014

@oliva's lab
"""
from .   import PPInnerContact, PHInnerContact
from SBI import SBIglobals


class InnerContacts(object):
    def __init__(self, pdb, AA = False, AA_type = "cb",  AA_distance = 12,
                            NC = False, NC_type = "min", NC_distance = 8,
                            HT = True,  HT_type = "min", HT_distance = 6):
        self._pdb          = pdb

        self._AAdo        = AA
        self._AAcontacts  = {}
        self._AA_type     = AA_type
        self._AA_distance = AA_distance

        self._NCdo        = NC
        self._NCcontacts  = {}
        self._NC_type     = NC_type
        self._NC_distance = NC_distance

        self._HTdo        = HT
        self._HTcontacts  = {}
        self._HT_type     = HT_type
        self._HT_distance = HT_distance

        self._build()

    #
    # ATTRIBUTES
    #
    @property
    def pdb(self):        return self._pdb

    @property
    def AAcontacts(self): return self._AAcontacts.values()

    @property
    def NCcontacts(self): return self._NCcontacts.values()

    @property
    def HTcontacts(self): return self._HTcontacts.values()

    #
    # PRIVATE FUNCTIONS
    #
    def _build(self):
        count  = 1

        SBIglobals.alert('debug', self, 'Analyzing Inner Contacts of {0:03} chains'.format(len(self.pdb)))

        for chain in self.pdb.chains:
            SBIglobals.alert('debug', self, 'Analyzing Chain {0:03} out of {1:03}'.format(count, len(self.pdb)))
            self._add_AAc(chain)
            self._add_NCc(chain)
            self._add_HTc(chain)
            count += 1

    def _add_AAc(self, chain):
        if self._AAdo and chain.chaintype == 'P':
            raise NotImplementedError
            ppi_id = PPInnerContact.test_identifier(chain)
            SBIglobals.alert('debug', self, 'Analyzing Protein Inner Contacts {0} for {1.chain}'.format(ppi_id, chain))

            ppi = PPInnerContact(chain, self._AA_type, self._AA_distance)

            SBIglobals.alert('debug', self, '\tAdding new Inner Contacts')
            self._AAcontacts[ppi_id] = ppi

    def _add_NCc(self, chain):
        if self._NCdo and chain.chaintype == 'N':
            raise NotImplementedError

    def _add_HTc(self, chain):
        if self._HTdo and chain.chaintype == 'P':
            phi_id = PHInnerContact.test_identifier(chain)
            SBIglobals.alert('debug', self, 'Analyzing Protein-Heteroatom Inner Contacts {0} for {1.chain}'.format(phi_id, chain))

            phi = PHInnerContact(chain, self._HT_type, self._HT_distance)

            SBIglobals.alert('debug', self, '\tAdding new Inner Contacts')
            self._HTcontacts[phi_id] = phi

    #
    # OVERRIDE DEFAULT METHODS
    #
    def __str__(self):
        data = []
        data.append("Inner Contacts for {0}".format(self.pdb.id))
        for ppi in self._AAcontacts:
            if len(self._AAcontacts[ppi]) > 0:
                data.append(self._AAcontacts[ppi].toString(True))
        for pni in self._NCcontacts:
            if len(self._NCcontacts[pni]) > 0:
                data.append(self._NCcontacts[pni].toString(True))
        for phi in self._HTcontacts:
            if len(self._HTcontacts[phi]) > 0:
                data.append(self._HTcontacts[phi].toString(True))
        return "\n".join(data)
