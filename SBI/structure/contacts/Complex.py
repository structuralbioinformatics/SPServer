"""
Complex

author: jbonet
date:   02/2014

@oliva's lab
"""
from .   import PPInterface, PNInterface, PHInterface
from SBI import SBIglobals


class Complex(object):
    def __init__(self, pdb, biomolecule = False, PPI = True, PPI_type = "cb",  PPI_distance = 12,
                                                 PNI = True, PNI_type = "min", PNI_distance = 8,
                                                 PHI = True, PHI_type = "min", PHI_distance = 6):
        self._pdb          = pdb
        self._biomolecule  = biomolecule

        self._PPIdo        = PPI
        self._PPInterface  = {}
        self._PPIblocked   = {}
        self._PPI_type     = PPI_type
        self._PPI_distance = PPI_distance

        self._PNIdo        = PNI
        self._PNInterface  = {}
        self._PNIblocked   = {}
        self._PNI_type     = PNI_type
        self._PNI_distance = PNI_distance

        self._PHIdo        = PHI
        self._PHInterface  = {}
        self._PHIblocked   = {}
        self._PHI_type     = PHI_type
        self._PHI_distance = PHI_distance

        self._build()

    #
    # ATTRIBUTES
    #
    @property
    def pdb(self):          return self._pdb

    @property
    def biomolecule(self):  return self._biomolecule

    @property
    def PPInterfaces(self): return self._PPInterface.values()

    @property
    def PNInterfaces(self): return self._PNInterface.values()

    @property
    def PHInterfaces(self): return self._PHInterface.values()

    #
    # PRIVATE FUNCTIONS
    #
    def _build(self):
        structures = []
        count      = 1

        if not self._biomolecule:
            structures.append(self._pdb)
        else:
            SBIglobals.alert('debug', self, 'Building biomolecules')
            structures.extend(self.pdb.apply_biomolecule_matrices(keepchains = True, water = False))

        SBIglobals.alert('debug', self, 'Analyzing Interfaces of {0:03} biomolecules'.format(len(structures)))

        for biom in structures:
            SBIglobals.alert('debug', self, 'Analyzing Biomolecule {0:03} out of {1:03}'.format(count, len(structures)))

            protein_chains = []
            protein_pgeoms = []
            protein_hgeoms = []
            for p in biom.proteins:
                protein_chains.append(p)
                protein_pgeoms.append(p.geometric_center(structure = True,  hetero = False, water = False, by_residue = True))
                protein_hgeoms.append(p.geometric_center(structure = False, hetero = True,  water = False, by_residue = True))
            nucleotide_chains = []
            nucleotide_ngeoms = []
            nucleotide_hgeoms = []
            for n in biom.nucleotides:
                nucleotide_chains.append(n)
                nucleotide_ngeoms.append(n.geometric_center(structure = True,  hetero = False, water = False, by_residue = True))
                nucleotide_hgeoms.append(n.geometric_center(structure = False, hetero = True,  water = False, by_residue = True))

            total_chains = 0
            if self._PPIdo: total_chains += len(protein_chains)
            if self._PNIdo: total_chains += len(nucleotide_chains)

            SBIglobals.alert('debug', self, '\tBiomolecule has {0:03} chains -> {1:03} max. Interfaces'.format(total_chains, (total_chains*(total_chains-1))/2))

            for i in xrange(len(protein_chains)):
                for j in xrange(i+1, len(protein_chains)):
                    self._add_PPI(protein_chains[i], protein_chains[j], protein_pgeoms[i], protein_pgeoms[j])
                for j in xrange(len(protein_chains)):
                    if i != j:
                        self._add_PHI(protein_chains[i], protein_chains[j], protein_pgeoms[i], protein_hgeoms[j])
                for j in xrange(len(nucleotide_chains)):
                    self._add_PNI(protein_chains[i], nucleotide_chains[j], protein_pgeoms[i], nucleotide_ngeoms[j])
                    self._add_PHI(protein_chains[i], nucleotide_chains[j], protein_pgeoms[i], nucleotide_hgeoms[j])

    def _add_PPI(self, chain1, chain2, geom1, geom2):
        if self._PPIdo:
            ppi_id = PPInterface.test_identifier(chain1, chain2)
            SBIglobals.alert('debug', self, 'Analyzing Protein-Protein Interface {0} for {1.chain} - {2.chain}'.format(ppi_id, chain1, chain2))
            if ppi_id not in self._PPIblocked or self._PPIblocked[ppi_id] < 2:
                SBIglobals.alert('debug', self, '\tInterface is NOT blocked')
                self._PPIblocked.setdefault(ppi_id, 0)
                ppi = PPInterface(chain1, chain2, self._PPI_type, self._PPI_distance, geom1, geom2)
                if ppi_id not in self._PPInterface:
                    SBIglobals.alert('debug', self, '\tAdding new Interface')
                    self._PPInterface[ppi_id] = ppi
                else:
                    SBIglobals.alert('debug', self, '\tUpdating Interface')
                    l = len(self._PPInterface[ppi_id])
                    self._PPInterface[ppi_id] + ppi
                    L = len(self._PPInterface[ppi_id])
                    if l == L:
                        SBIglobals.alert('debug', self, '\t\tInterface does NOT give any new contact')
                        self._PPIblocked[ppi_id] += 1
                    else:
                        SBIglobals.alert('debug', self, '\t\tInterface GIVES new contacts')
            else:
                SBIglobals.alert('debug', self, '\tInterface IS blocked')

    def _add_PNI(self, chain1, chain2, geom1, geom2):
        if self._PNIdo:
            pni_id = PNInterface.test_identifier(chain1, chain2)
            SBIglobals.alert('debug', self, 'Analyzing Protein-Nucleotide Interface {0} for {1.chain} - {2.chain}'.format(pni_id, chain1, chain2))
            if pni_id not in self._PNIblocked or self._PNIblocked[pni_id] < 2:
                SBIglobals.alert('debug', self, '\tInterface is NOT blocked')
                self._PNIblocked.setdefault(pni_id, 0)
                pni = PNInterface(chain1, chain2, self._PNI_type, self._PNI_distance, geom1, geom2)
                if pni_id not in self._PNInterface:
                    SBIglobals.alert('debug', self, '\tAdding new Interface')
                    self._PNInterface[pni_id] = pni
                else:
                    SBIglobals.alert('debug', self, '\tUpdating Interface')
                    l = len(self._PNInterface[pni_id])
                    self._PNInterface[pni_id] + pni
                    L = len(self._PNInterface[pni_id])
                    if l == L:
                        SBIglobals.alert('debug', self, '\t\tInterface does NOT give any new contact')
                        self._PNIblocked[pni_id] += 1
                    else:
                        SBIglobals.alert('debug', self, '\t\tInterface GIVES new contacts')
            else:
                SBIglobals.alert('debug', self, '\tInterface IS blocked')

    def _add_PHI(self, chain1, chain2, geom1, geom2):
        if self._PHIdo:
            phi_id = PHInterface.test_identifier(chain1, chain2)
            SBIglobals.alert('debug', self, 'Analyzing Protein-Heteroatom Interface {0} for {1.chain} - {2.chain}'.format(phi_id, chain1, chain2))
            if phi_id not in self._PHIblocked or self._PHIblocked[phi_id] < 2:
                SBIglobals.alert('debug', self, '\tInterface is NOT blocked')
                self._PHIblocked.setdefault(phi_id, 0)
                phi = PHInterface(chain1, chain2, self._PHI_type, self._PHI_distance, geom1, geom2)
                if phi_id not in self._PHInterface:
                    SBIglobals.alert('debug', self, '\tAdding new Interface')
                    self._PHInterface[phi_id] = phi
                else:
                    SBIglobals.alert('debug', self, '\tUpdating Interface')
                    l = len(self._PHInterface[phi_id])
                    self._PHInterface[phi_id] + phi
                    L = len(self._PHInterface[phi_id])
                    if l == L:
                        SBIglobals.alert('debug', self, '\t\tInterface does NOT give any new contact')
                        self._PHIblocked[phi_id] += 1
                    else:
                        SBIglobals.alert('debug', self, '\t\tInterface GIVES new contacts')
            else:
                SBIglobals.alert('debug', self, '\tInterface IS blocked')

    #
    # OVERRIDE DEFAULT METHODS
    #
    def __str__(self):
        data = []
        data.append("Contacts for {0}".format(self.pdb.id))
        for ppi in self._PPInterface:
            if len(self._PPInterface[ppi]) > 0:
                data.append(self._PPInterface[ppi].toString(True))
        for pni in self._PNInterface:
            if len(self._PNInterface[pni]) > 0:
                data.append(self._PNInterface[pni].toString(True))
        for phi in self._PHInterface:
            if len(self._PHInterface[phi]) > 0:
                data.append(self._PHInterface[phi].toString(True))
        return "\n".join(data)
