'''
@file: DSSPExe.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2013

@ [oliva's lab](http://sbi.imim.es)

@class: DSSPExe
'''
import os
import sys

from SBI.external import ExternalExe
from DSSP         import DSSP
from SBI          import SBIglobals as SBIg
from SBI.beans    import File


class DSSPExe(ExternalExe):
    '''
    Manages the execution of DSSP, a secondary structure determination program
    from PDB structures.

    '''
    def __init__(self, pdb_file, dssp_file, cleanpdb = False, cleandssp = False):
        '''
        @param:    pdb_file
        @pdef:     name of the PDB file to calculate SS from
        @ptype:    {String}

        @param:    dssp_file
        @pdef:     name of the DSSP output file
        @ptype:    {String}

        @param:    cleanpdb
        @pdef:     flag. delete pdb file after execution
        @pdefault: _False_
        @ptype:    {Boolean}

        @param:    cleandssp
        @pdef:     flag. delete dssp file after execution
        @pdefault: _False_
        @ptype:    {Boolean}
        '''
        self._pdbfile  = pdb_file
        self._dsspfile = dssp_file
        self._dsspdata = []
        self._gapped   = False

        if DSSPExe._EXE is None:
            self._set_default_executable('dssp')

        self._execute()
        self._parse()
        self._clean(cleanpdb, cleandssp)

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def dsspdata(self):
        '''
        List with the DSSP result.
        One residue-data per position.

        @return: {List}
        '''
        return self._dsspdata

    @property
    def empty_dssp(self):
        '''
        Creates a empty DSSP object to fill the gaps.

        @return: {DSSP}
        '''
        return DSSP(secondary_structure = '-',
                    accessibility       = 1,
                    amino               = 'X')

    ############
    # BOOLEANS #
    ############
    @property
    def gapped(self):
        '''
        Checks if dssp has found gaps in the structure.

        @return: {Boolean}
        '''
        return self._gapped

    ###########
    # METHODS #
    ###########
    @staticmethod
    def dynamic(executable, path):
        '''
        Manually set the values for dssp path and executable .

        @param:    executable
        @pdef:     name of the executable file
        @ptype:    {String}

        @param:    path
        @pdef:     path to the executable file
        @ptype:    {String}
        '''
        DSSPExe._set_dynamic_executable(executable, path)

    ###################
    # PRIVATE METHODS #
    ###################
    def _execute(self):
        '''
        Executes the DSSP call.
        '''
        self._EXE.add_parameter(self._pdbfile)
        self._EXE.add_parameter(self._dsspfile)
        try:
            self._EXE.execute(silent=True)
        except SystemError, e:
            msg = 'Some error occurred while executing dssp\n{0}\n'.format(e)
            SBIg.throw(self, msg, e)

        self._EXE.clean_command()

    def _parse(self):
        file_fd    = File(self._dsspfile)
        read       = False
        continuity = -1000
        readline   = 0
        for line in file_fd.read():
            if line.startswith("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC"):
                read = True
                continue
            if read:
                if line[13:14] != '!':
                    res_num = int(line[6:10].strip())
                    ss      = line[16:17] if line[16:17] != ' ' else '-'
                    buried  = int(line[35:38].strip())
                    aa      = line[13:15].strip()

                    self._dsspdata.append(DSSP(secondary_structure = ss,
                                               accessibility       = buried,
                                               amino               = aa))
                    self._dsspdata[-1].add_hydrogen_links(line[39:50],
                                                          line[50:61],
                                                          line[61:72],
                                                          line[72:84])
                    if readline > 0:
                        if res_num != continuity + 1:
                            self._gapped = True
                        continuity = res_num
                    readline += 1
                else:
                    msg = "truncated chain!{0}\n".format(self._dsspfile)
                    sys.stderr.write(msg)
                    SBIg.warn(self, msg)
                    self._gapped = True
        file_fd.close()

    def _clean(self, cleanpdb, cleandssp):
        '''
        If required, removes the files used in the DSSP execution.

        @param:    cleanpdb
        @pdef:     flag. delete pdb file after execution
        @ptype:    {Boolean}

        @param:    cleandssp
        @pdef:     flag. delete dssp file after execution
        @ptype:    {Boolean}
        '''
        if cleanpdb:
            os.unlink(self._pdbfile)
        if cleandssp:
            os.unlink(self._dsspfile)
