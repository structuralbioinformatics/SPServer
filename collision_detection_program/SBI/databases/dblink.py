'''
@file: dblink.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2014

@ [oliva's lab](http://sbi.imim.es)

@class: DBlink
'''
from abc import ABCMeta, abstractmethod
from datetime import date
# import time
import urllib

import os
import json

from SBI.beans import Path
from SBI.beans import File
from SBI       import SBIglobals as SBIg


class DBlink(object):
    '''
    Manages the connection to databases.
    It creates a local copy of the DB that can be queried a posteriori.
    '''
    __metaclass__ = ABCMeta

    _MANDATORY_FILES = []
    _ITEM_FILES      = []
    _SOURCES         = []
    _SOURCES_SIZES   = []
    _SHOW_LINK       = None

    _DBOBJECT        = None

    _CONTROL_FILE    = 'release.json.gz'
    _TODAY           = date.today().isoformat()
    _RELEASE         = {'db':            None,
                        'date':          _TODAY,
                        'total_items':   {},
                        'new_items':     {},
                        'update_items':  {},
                        'deleted_items': {}}

    def __init__(self, local):
        '''
        @param:    local
        @pdef:     local directory to store the database. Create if not exist
        @ptype:    {String}
        '''
        self._local = local
        Path.mkdir(self._local)

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def local(self):
        '''
        Directory of the local database

        @return: {String}
        '''
        return os.path.abspath(self._local)

    @property
    def release(self):
        '''
        Retrieves release data for the database.
        Not according to the DB release, but to when we downloaded it.

        @returns: {Dictionary}
        '''
        if os.path.isfile(os.path.join(self.local, self._CONTROL_FILE)):
            f    = File(os.path.join(self.local, self._CONTROL_FILE))
            data = json.loads(f.read())
            f.close()
        else:
            data = self._RELEASE
        return data

    @property
    def items(self):
        '''
        Loops through the items of the database

        @yields: Object depending on the database.
        '''
        if not self.has_local:
            SBIg.throw(self, 'A local database needs to be build first', IOError)

        for ifile in self._ITEM_FILES:
            ifile = os.path.join(self.local, ifile)
            f     = File(ifile)
            for line in f.read():
                yield self._DBOBJECT.grab(line.strip())
            f.close()

    ############
    # BOOLEANS #
    ############
    @property
    def has_local(self):
        '''
        Checks whether we already have a local database.

        @return: {Boolean}
        '''
        if not bool(len(self._MANDATORY_FILES)):
            return os.path.isdir(self.local)
        else:
            for mfile in self._MANDATORY_FILES:
                if not os.path.isfile(os.path.join(self.local, mfile)):
                    return False
            return True

    ###########
    # METHODS #
    ###########
    def create(self):
        '''
        Create the local database.
        Returns _True_ on success, _False_ otherwise

        @return: {Boolean}
        '''
        if self.has_local:
            SBIg.warn(self, 'A local copy exist. Executing update.')
            return self.update()

        self._download_sources()
        self._process()
        self._save_release()
        # self._clean_sources()

        return True

    def update(self):
        '''
        Update the local database.
        Returns _True_ on success, _False_ otherwise

        @return: {Boolean}
        '''
        self._RELEASE = self.release
        self._RELEASE['new_items']     = {}
        self._RELEASE['update_items']  = {}
        self._RELEASE['deleted_items'] = {}
        self._download_sources()
        self._process(True)
        self._save_release()
        # self._clean_sources()

        return True

    ###################
    # PRIVATE METHODS #
    ###################
    @abstractmethod
    def _process(self, update = False):
        '''
        Transform the source files into the final local db files.

        @param:    update
        @pdef:     toggles between create and update processing
        @pdefault: _False_
        @ptype:    {Boolean}
        '''
        raise NotImplemented

    def _save_release(self):
        '''
        Store the release data into a file.
        '''
        f    = File(os.path.join(self.local, self._CONTROL_FILE), 'w', True)
        f.write(json.dumps(self._RELEASE))
        f.close()

    def _download_sources(self):
        '''
        Download the source files to local directory.
        '''
        for dfile in self._SOURCES:
            download    = False
            source      = os.path.join(self._FTP, dfile)
            source_size = DBlink._file_size(source)
            # source_date = DBlink._file_date(source)
            DBlink._SOURCES_SIZES.append(source_size)
            destination = os.path.join(self.local, dfile)
            if not os.path.isfile(destination):
                if DBlink._RELEASE['date'] != DBlink._TODAY:
                    # if source_date < DBlink._RELEASE:
                    #     SBIg.alert('verbose', self, 'No new updates in the source side' +
                    #                                 ' for {0}'.format(source) +
                    #                                 ' since the last local update.')
                    # else:
                    download = True
                else:
                    download = True
            else:
                SBIg.alert('verbose', self, 'Looks like {0} has already been downloaded.'.format(source))

            if download:
                SBIg.alert('verbose', self, 'Downloading {0} to {1}'.format(source, destination))
                SBIg.alert('verbose', self, 'Source file size is {0:.3f} MB.'.format(source_size))
                urllib.urlretrieve(source, destination)

    def _clean_sources(self):
        '''
        Remove source files.
        '''
        for dfile in self._SOURCES:
            destination = os.path.join(self.local, dfile)
            if not SBIg.debug:
                os.unlink(destination)

    @staticmethod
    def _file_size(external_file):
        '''
        Determines the size of a remote file. In Mega Bytes.

        @param:    external_file
        @pdef:     remote file
        @ptype:    {String}

        @return: {Float}
        '''
        usock = urllib.urlopen(external_file)
        size  = usock.info().get('Content-Length')
        if size is None:
            return 0
        return float(size) / (1024.0 * 1024.0)

    # @staticmethod
    # def _file_date(external_file):
    #     '''
    #     Determines the date of a remote file.

    #     @param:    external_file
    #     @pdef:     remote file
    #     @ptype:    {String}

    #     @return: {Time}
    #     '''
    #     usock = urllib.urlopen(external_file)
    #     lastm = usock.info().get('last-modified')
    #     return time.strptime(lastm, '%a, %d %b %Y %H:%M:%S %Z')
