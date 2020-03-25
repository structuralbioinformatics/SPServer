import json

from SBI.beans    import JSONer
from SBI.sequence import Sequence


class Uniprot(JSONer):

    def __init__(self, entry_name, status):
        '''
        @param:    entry_name
        @pdef:     entry name of the new uniprot entry
        @ptype:    {String}

        @param:    status
        @pdef:     To distinguish the fully annotated entries in the Swiss-Prot
                   section of the UniProt Knowledgebase from the computer-annotated
                   entries in the TrEMBL section
        @ptype:    {String}
        '''
        self._entry     = entry_name
        self._status    = status
        self._version   = 0
        self._accession = []
        self._taxid     = None
        self._hostid    = []
        self._referals  = {}
        self._sequence  = ''

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def entry_name(self):
        '''
        Uniprot unique identifier

        @return: {String}
        '''
        return self._entry

    @property
    def section(self):
        '''
        To which section the entry belongs.

        @return: {String}
        '''
        return 'swissprot' if self.is_swissprot else 'trembl'

    @property
    def version(self):
        '''
        actual version of the entry

        @return: {Integer}
        '''
        return self._version

    @version.setter
    def version(self, value):
        '''
        @param:    value
        @pdef:     version of the entry
        @ptype:    {Integer}
        '''
        self._version = int(value)

    @property
    def accession(self):
        '''
        List of accession numbers for the entry

        @return: {List}
        '''
        return self._accession

    @accession.setter
    def accession(self, value):
        '''
        @param:    value
        @pdef:     list of accession numbers
        @ptype:    {List}
        '''
        self._accession.extend(value)

    @property
    def taxid(self):
        '''
        Specie to which the protein belongs.

        @return: {String}
        '''
        return self._taxid

    @taxid.setter
    def taxid(self, value):
        '''
        @param:    value
        @pdef:     Specie to which the protein belongs.
        @ptype:    {Integer}
        '''
        self._taxid = int(value)

    @property
    def host(self):
        '''
        List of species to which the protein specie infects.

        @return: {List}
        '''
        return self._hostid

    @host.setter
    def host(self, value):
        '''
        @param:    value
        @pdef:     Specie to which the protein specie infects.
        @ptype:    {Integer}
        '''
        self._hostid.append(int(value))

    @property
    def referals(self):
        '''
        References from other databases.
        Altough being all together, each reference has its specificities.
        Check them in the [Uniprot manual](http://web.expasy.org/docs/userman.html#DR_line)

        @return: {List} of {Dictionary}
        '''
        return self._referals

    @referals.setter
    def referals(self, value):
        '''
        @param:    value
        @pdef:     database reference info
        @ptype:    {List}
        '''
        db = value.pop(0)
        self._referals.setdefault(db, [])
        reference = {}
        reference['id'] = value.pop(0)
        if bool(len(value)):
            opt1 = value.pop(0).strip()
            if db == 'RefSeq':
                reference['nuc_seq_id'] = opt1
            elif db in ['InterPro', 'PANTHER', 'Pfam', 'PIR', 'PRINTS', 'ProDom',
                        'REBASE', 'SMART', 'SUPFAM', 'TIGRFAMs', 'PeroxiBase',
                        'PROSITE']:
                reference['name'] = opt1
            elif db == 'PDB':
                reference['method'] = opt1
            elif db in ['dictyBase', 'EcoGene', 'FlyBase', 'CGD', 'HGNC', 'MGI',
                        'PATRIC', 'RGD', 'SGD', 'Xenbase', 'ZFIN']:
                reference['gene_designation'] = opt1
            elif db == 'GO':
                reference['ontology'] = opt1.split(':')[0]
            elif db in ['HAMAP', 'PIRSF', 'CAZy', 'GeneFarm', 'TCDB']:
                reference['family'] = opt1
            elif db in ['WormBase', 'Ensembl', 'EnsemblBacteria', 'EnsemblFungi',
                        'EnsemblMetazoa', 'EnsemblPlants', 'EnsemblProtists',
                        'EMBL']:
                reference['protein_id'] = opt1
            elif db == 'VectorBase':
                reference['species_of_origin'] = opt1
            elif db in ['IntAct', 'BioGrid']:
                reference['num_interactors'] = opt1
            elif db == 'SMR':
                reference['range'] = opt1
            elif db == 'MIM':
                reference['type'] = opt1
            elif db == 'Allergome':
                reference['allergen'] = opt1
            elif db in ['ArachnoServer', 'ConoServer']:
                reference['toxin'] = opt1
            elif db == 'Orphanet':
                reference['disease'] = opt1
            elif db in ['ChiTaRS', 'UCSC', 'BRENDA']:
                reference['organism'] = opt1
            elif db == 'UniPathway':
                reference['reaction'] = opt1
            elif db == 'DrugBank':
                reference['drug'] = opt1
            if bool(len(value)):
                opt2 = value.pop(0).strip()
                if db in ['Gene3D', 'HAMAP', 'PANTHER', 'Pfam', 'PIRSF', 'ProDom',
                          'SMART', 'SUPFAM', 'TIGRFAMs', 'PROSITE']:
                    reference['hits'] = opt2.rstrip('.')
                elif db in ['WormBase', 'Ensembl', 'EnsemblBacteria', 'EnsemblFungi',
                            'EnsemblMetazoa', 'EnsemblPlants', 'EnsemblProtists']:
                    reference['gene_id'] = opt2
                elif db == 'GO':
                    reference['evidence'] = opt2.split(':')[0]
                elif db == 'PDB':
                    reference['resolution'] = opt2
                elif db == 'EMBL':
                    reference['status'] = opt2
                if bool(len(value)):
                    opt3 = value.pop(0).strip()
                    if db == 'PDB':
                        reference['ranges'] = opt3
                    elif db == 'WormBase':
                        reference['gene_designation'] = opt3
                    elif db == 'EMBL':
                        reference['molecule_type'] = opt3
        self._referals[db].append(reference)

    @property
    def sequence(self):
        '''
        Protein sequence

        @return: {Sequence}
        '''
        return Sequence(self.entry_name, self._sequence)

    @sequence.setter
    def sequence(self, value):
        '''
        @param:    value
        @pdef:     section of the protein sequence
        @ptype:    {String}
        '''
        self._sequence += value

    ############
    # BOOLEANS #
    ############
    @property
    def is_swissprot(self):
        '''
        Checks if the entry belongs to Swiss-Prot

        @return: {Boolean}
        '''
        return self._status == 'Reviewed'

    @property
    def is_trembl(self):
        '''
        Checks if the entry belongs to TrEMBL

        @return: {Boolean}
        '''
        return self._status == 'Unreviewed'

    def has_reference_to(self, external_db):
        '''
        Check if the protein has a reference to a given database.

        @param:    external_db
        @pdef:     identifier of the external database
        @ptype:    {String}

        @return: {Boolean}
        '''
        return external_db in self.referals

    ###########
    # METHODS #
    ###########
    @staticmethod
    def grab(json_line):
        '''
        Retrieve the object data from a json line.

        @param:    json_line
        @pdef:     line of data.
        @ptype:    {String} or {Dictionary}

        @return: {Uniprot}
        '''
        if isinstance(json_line, basestring):
            json_line = json.loads(json_line.strip())

        data = Uniprot('', '')
        data.__dict__ = json_line

        return data

    def as_dict(self):
        return self.__dict__
