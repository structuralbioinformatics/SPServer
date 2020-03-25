'''
@file: Path.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   02/2013

@ [oliva's lab](http://sbi.imim.es)

@class: Path
'''
import os
import fnmatch
import shutil

from .  import FileError
from .. import SBIglobals as SBIg


class Path(object):
    '''
    Collection of static methods (can be called without declaring an instance
    of the class) to help in file and directory management.
    They are defined this way and not as single methods to avoid overwrite
    system functions.
    '''

    @staticmethod
    def list_files(root              = os.getcwd(), pattern  = '*',
                   avoid_empty_files = True,        rootless = False):
        '''
        Returns all file names in a directory tree matching a specific pattern.

        @param:    root
        @pdef:     Root of the directory tree to search.
        @pdefault: Current working directory.
        @ptype:    {String}

        @param:    pattern
        @pdef:     Expression(s) to match (ls-like format)
        @pdefault: '*'
        @ptype:    {String} or {List}

        @param:    avoid_empty_files
        @pdef:     Ignore files with size 0.
        @pdefault: _True_
        @ptype:    {Boolean}

        @param:    rootless
        @pdef:     When False, the name of the files are returned with
                   absolute path. Otherwise the root is removed.
        @pdefault: _False_
        @ptype:    {Boolean}

        @raises: {FileError} if root is a file
        @yields: {String}.
        '''

        if os.path.isfile(root):
            error = FileError(root)
            raise error.is_file()

        search_patterns = []
        if not isinstance(pattern, list):
            search_patterns.append(pattern)
        else:
            search_patterns = pattern

        for pat in search_patterns:
            for path, dirs, files in os.walk(os.path.abspath(root)):
                for filename in fnmatch.filter(files, pat):
                    qfile = os.path.join(path, filename)
                    size = os.path.getsize(qfile)
                    if not avoid_empty_files or size > 0:
                        msg = 'Found file {0}'.format(qfile)
                        SBIg.alert('debug', Path(), msg)
                        if not rootless:
                            yield qfile
                        else:
                            root = os.path.abspath(root) + "/"
                            yield qfile.replace(root, '')

    @staticmethod
    def list_directories(root = os.getcwd(), rootless = False):
        '''
        Returns all directory names in a directory tree.

        @param:    root
        @pdef:     Root of the directory tree to search.
        @pdefault: Current working directory.
        @ptype:    {String}

        @param:    rootless
        @pdef:     When False, the name of the directories are returned with
                   absolute path. Otherwise the root is removed.
        @pdefault: _False_
        @ptype:    {Boolean}


        @raises: {FileError} if root is a file
        @yields: {String}.
        '''

        if os.path.isfile(root):
            error = FileError(root)
            raise error.is_file()

        for path, dirs, files in os.walk(os.path.abspath(root)):
            for onedir in dirs:
                qdir = os.path.join(path, onedir)
                SBIg.alert('debug', Path(), 'Found directory {0}'.format(qdir))
                if not rootless:
                    yield qdir
                else:
                    r = root if root.endswith('/') else root + '/'
                    yield qdir.replace(r, '')

    @staticmethod
    def sync_directories(source_dir, destination_dir, by_dir = True,
                         by_file = False, only_list = False):
        '''
        Synchronizes two directory trees.

        @param:    source_dir
        @pdef:     Name of the original directory.
        @ptype:    {String}

        @param:    destination_dir
        @pdef:     Name of the synchronized directory.
        @ptype:    {String}

        @param:    by_dir
        @pdef:     Sync is performed at directory level.
        @pdefault: _True_
        @pclash:   by_file
        @ptype:    {Boolean}

        @param:    by_file
        @pdef:     Sync is performed at file level.
        @pdefault: _False_
        @pclash:   by_dir
        @ptype:    {Boolean}

        @param:    only_list
        @pdef:     When True, sync is not performed, only listed.
        @pdefault: _False_
        @ptype:    {Boolean}

        @raises: {AttributeError} if by_dir and by_file are both True.
        @raises: {FileError} if either source_dir or destination_dir are files.
        @yields: {String}. Names of the synchronized files/directories
        '''

        if by_file is True and by_dir is True:
            raise AttributeError('Both sync methods can not be active\n')
        if os.path.isfile(source_dir):
            error = FileError(source_dir)
            raise error.is_file()
        if os.path.isfile(destination_dir):
            error = FileError(destination_dir)
            raise error.is_file()

        sdir = os.path.abspath(source_dir)
        ddir = os.path.abspath(destination_dir)
        pdir = '#not_a_dir_at_all'

        source_dirs = set()
        if by_dir:
            for onedir in Path.list_directories(root = sdir, rootless = True):
                source_dirs.add(onedir)
            for onedir in Path.list_directories(root = ddir, rootless = True):
                if onedir in source_dirs:
                    source_dirs.remove(onedir)
        if by_file:
            for onedir in Path.list_files(root = sdir, rootless = True):
                source_dirs.add(onedir)
            for onedir in Path.list_files(root = ddir, rootless = True):
                if onedir in source_dirs:
                    source_dirs.remove(onedir)

        for onedir in sorted(source_dirs):
            onedir = onedir.lstrip('/')
            if not os.path.commonprefix([pdir, onedir]) == pdir:
                msg = '{0} is different from {1} to {2}'.format(onedir,
                                                                sdir, ddir)
                SBIg.alert('verbose', Path(), msg)
                fullori = os.path.join(sdir, onedir)
                fullnew = os.path.join(ddir, onedir)
                if not only_list:
                    shutil.copytree(fullori, fullnew)
                    yield fullori
                else:
                    yield fullori
                pdir = onedir

    @staticmethod
    def do_files_match(root = os.getcwd(), pattern = '*',
                       avoid_empty_files = True):
        '''
        Determines if there is any file in a directory tree matching a
        specific pattern.

        @param:    root
        @pdef:     Root of the directory tree to search.
        @pdefault: Current working directory.
        @ptype:    {String}

        @param:    pattern
        @pdef:     Expression(s) to match (ls-like format)
        @pdefault: '*'
        @ptype:    {String} or {List}

        @param:    avoid_empty_files
        @pdef:     Ignore files with size 0.
        @pdefault: _True_
        @ptype:    {Boolean}

        @raises: {FileError} if root is a file
        @return: {Boolean}
        '''
        if os.path.isfile(root):
            error = FileError(root)
            raise error.is_file()

        search_patterns = []
        if not isinstance(pattern, list):
            search_patterns.append(pattern)
        else:
            search_patterns = pattern

        for pat in search_patterns:
            for path, dirs, files in os.walk(os.path.abspath(root)):
                for filename in fnmatch.filter(files, pat):
                    size = os.path.getsize(os.path.join(path, filename))
                    if not avoid_empty_files or size > 0:
                        return True
        return False

    @staticmethod
    def mkdir(new_dir):
        '''
        Creates a new directory.
        Ignores the command if the directory already exists.
        Creates parent directories as required.

        @recipe: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/82465

        @param:    new_dir
        @pdef:     Name of the new directory.
        @ptype:    {String}

        @raises: {FileError} if the directory can not be created
        '''
        error = FileError(os.path.abspath(new_dir))
        if os.path.isdir(new_dir):
            SBIg.warn(Path(), 'A directory with the name {0} already exists.'.format(new_dir))
        elif os.path.isfile(new_dir):
            raise error.is_file()
        else:
            SBIg.alert('debug', Path(), 'Creating dir {0}.'.format(new_dir))
            head, tail = os.path.split(new_dir)
            if head and not os.path.isdir(head):
                Path.mkdir(head)
            if tail:
                os.mkdir(new_dir)

    @staticmethod
    def copy_file(source, destination, overwrite = False):
        '''
        Copy a file from one place to another.

        @param:    source
        @pdef:     File to copy.
        @ptype:    {String}

        @param:    destination
        @pdef:     Copy file.
        @ptype:    {String}

        @param:    overwrite
        @pdef:     Overwrite files if they already exist.
        @pdefault: _SBIglobals.overwrite_
        @ptype:    {Boolean}

        @raises: {FileError} if source does not exist or destination exists
                 and can not be overwritten
        '''
        # Local overwrite takes precedence over Global overwrite
        overwrite = SBIg.decide_overwrite(overwrite)

        if not os.path.exists(source):
            error = FileError(source)
            raise error.file_does_not_exist()
        if os.path.exists(destination) and not overwrite:
            error = FileError(destination)
            raise error.file_overwrite()

        shutil.copyfile(src = source, dst = destination)
