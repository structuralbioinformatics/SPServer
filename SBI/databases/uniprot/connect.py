import os
import re

from SBI.databases import DBlink
from uniprot       import Uniprot
from SBI.beans     import File
from SBI           import SBIglobals as SBIg


class Connect(DBlink):
    '''
    Manages the connection to the [Uniprot](http://www.uniprot.org/) database.
    '''
    _FTP = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/'
    _MANDATORY_FILES = ['swisprot.gz', 'swisprot.fasta',
                        'trembl.gz',   'trembl.fasta']
    _ITEM_FILES      = ['swisprot.gz', 'trembl.gz']
    _SOURCES         = ['uniprot_sprot.dat.gz', 'uniprot_trembl.dat.gz']
    _SHOW_LINK       = 'http://www.uniprot.org/'

    _DBOBJECT        = Uniprot

    def __init__(self, local):
        '''
        @param:    local
        @pdef:     local directory to store the database. Create if not exist
        @ptype:    {String}
        '''
        super(Connect, self).__init__(local)
        self._RELEASE['db'] = 'uniprot'

    ###################
    # PRIVATE METHODS #
    ###################
    def _process(self, update = False):
        '''
        Transform the source files into the final local db files.

        @param:    update
        @pdef:     toggles between create and update processing
        @pdefault: _False_
        @ptype:    {Boolean
        '''
        if update:
            old = self._RELEASE['total_items'].copy()
        j = 0
        for i in range(len(self._SOURCES)):
            dfilen  = os.path.join(self.local, self._SOURCES[i])
            ofilen  = os.path.join(self.local, self._MANDATORY_FILES[j])
            ffilen  = os.path.join(self.local, self._MANDATORY_FILES[j + 1])
            if not os.path.isfile(dfilen):
                continue
            SBIg.alert('verbose', self, 'Parsing:       {0}'.format(dfilen))
            SBIg.alert('verbose', self, 'DB file to:    {0}'.format(ofilen))
            SBIg.alert('verbose', self, 'Fasta file to: {0}'.format(ffilen))
            dfile   = File(dfilen)
            ofile   = File(ofilen, 'w', update)
            ffile   = File(ffilen, 'w', update)
            protein = None
            for protein in Connect._parse_uniprot(dfile):
                pname = protein.entry_name
                pvers = protein.version
                SBIg.alert('verbose', self, 'Protein: {0}'.format(pname))
                if not update:
                    self._RELEASE['total_items'][pname] = pvers
                else:
                    if pname not in self._RELEASE['total_items']:
                        self._RELEASE['new_items'][pname] = pvers
                    else:
                        del(old[pname])
                        if self._RELEASE['total_items'][pname] != pvers:
                            self._RELEASE['update_items'][pname] = pvers

                ffile.write(protein.sequence.format('FASTA') + '\n')
                ofile.write(protein.json() + '\n')
            j += 2
            dfile.close()
            ofile.close()
            ffile.close()

        if update:
            self._RELEASE['total_items'].update(self._RELEASE['new_items'])
            self._RELEASE['total_items'].update(self._RELEASE['update_items'])
            self._RELEASE['deleted_items'] = old
            for k in self._RELEASE['deleted_items']:
                del(self._RELEASE['total_items'][k])

    @staticmethod
    def _parse_uniprot(uniprot_file):
        '''
        Reads a uniprot file and generates the different uniprot entries

        @param:    uniprot_file
        @pdef:     uniprot file to parse
        @ptype:    {File}

        @yields: {Uniprot}
        '''
        version = re.compile('DT\s+\d+\-\w+\-\d+\,\sentry\sversion\s(\d+)')
        taxid   = re.compile('NCBI_TaxID\=(\d+);')
        protein = None
        for line in uniprot_file.read():
            line = line.rstrip()
            if line.startswith('ID'):
                protein = Uniprot(line.split()[1].strip(),
                                  line.split()[2].strip(' ;'))
            if line.startswith('DT'):
                m = version.search(line)
                if m:
                    protein.version = m.group(1)
            if line.startswith('AC'):
                protein.accession = [x.strip(';') for x in line.split()[1:]]
            if line.startswith('OX'):
                protein.taxid = taxid.search(line.split()[1]).group(1)
            if line.startswith('OH'):
                protein.host = taxid.search(line.split()[1]).group(1)
            if line.startswith('DR'):
                protein.referals = ''.join(line.split()[1:]).split(';')
            if line.startswith('  '):
                protein.sequence = line.replace(' ', '')
            if line.startswith('//'):
                yield protein
