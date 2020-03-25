'''
@file: IndexedNum.py

@author: Jaume Bonet
@mail:   jaume.bonet@gmail.com
@date:   2013

@ [oliva's lab](http://sbi.imim.es)

@class: IndexedNum
'''
import re


class IndexedNum(object):
    '''
    Created for PDB.
    Some PDB files contain residues in positions that are not only numeric,
    but also have an index.
    {IndexedNum} tries to palliate this by giving the properties of a number
    (as much as possible) to this kind of number-index combination.
    Thus, a 7 becomes a '7 ' but can still be used as a number.
    Comparisons, addition (to integers or fully integer {IndexedNum}) and int(),
    float() casts are implemented. It is hash-able and, thus, can be used as key
    in dictionaries.

    If it is initialized with X, it means that we have a "blank" position.
    To avoid mix-ups, the index is transformed to '_X_'
    '''
    _NUMBER_REGEX = re.compile('(\-*\d+)(\w*)')
    _REVERS_REGEX = re.compile('(\w{1})(\d*)')

    def __init__(self, value):
        '''
        @param:    value
        @pdef:     value to transform to a {IndexedNum}
        @ptype:    {Integer}, {String}, {IndexedNum}
        '''
        if isinstance(value, IndexedNum):
            self._num = value._num
            self._idx = value._idx
        else:
            value = str(value)
            if not value.startswith('X'):
                parts = re.search(self._NUMBER_REGEX, value)
                self._num = int(parts.group(1))
                self._idx = parts.group(2)
                if self._idx == '':
                    self._idx = ' '
            else:
                parts = re.search(self._REVERS_REGEX, value)
                if parts:
                    self._num = int(parts.group(2))
                else:
                    self._num = None
                self._idx = '_X_'

    ##############
    # ATTRIBUTES #
    ##############
    @property
    def number(self):
        '''
        Integer part value.

        @return: {Integer}
        '''
        return self._num

    @property
    def index(self):
        '''
        String part value.

        @return: {String}
        '''
        return self._idx

    ############
    # BOOLEANS #
    ############
    @property
    def is_integer(self):
        '''
        Is the value a pure Integer?

        @return: {Boolean}
        '''
        return self._idx == ' '

    @property
    def is_blank(self):
        '''
        Does it have no value?

        @return: {Boolean}
        '''
        return self._idx == '_X_'

    #################
    # CLASS METHODS #
    #################
    def __lt__(self, other):
        if isinstance(other, int):
            return self.number < other
        elif isinstance(other, IndexedNum):
            if self.number < other.number:
                return True
            elif self.number == other.number:
                return self.index < other.index
            else:
                return False
        return NotImplemented

    def __le__(self, other):
        if isinstance(other, int):
            return self.number <= other
        elif isinstance(other, IndexedNum):
            if self.number < other.number:
                return True
            elif self.number == other.number:
                return self.index <= other.index
            else:
                return False
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, int):
            return self.number == other
        elif isinstance(other, IndexedNum):
            return self.number == other.number and self.index == other.index
        elif isinstance(other, basestring):
            return str(self).strip() == other.strip()
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __gt__(self, other):
        if isinstance(other, int):
            return self.number > other
        elif isinstance(other, IndexedNum):
            if self.number > other.number:
                return True
            elif self.number == other.number:
                return self.index > other.index
            else:
                return False
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, int):
            return self.number >= other
        elif isinstance(other, IndexedNum):
            if self.number > other.number:
                return True
            elif self.number == other.number:
                return self.index >= other.index
            else:
                return False
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, int):
            return IndexedNum(str(self.number + other) + self.index)
        if isinstance(other, IndexedNum) and other.is_integer:
            return IndexedNum(str(self.number + other) + self.index)
        return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, int):
            self.number += other
        if isinstance(other, IndexedNum) and other.is_integer:
            self.number += other
        return NotImplemented

    def __int__(self):
        return int(self.number)

    def __float__(self):
        return float(self.number)

    def __hash__(self):
        return repr(self)

    def __repr__(self):
        if self.is_blank:
            if self._num is None:
                return 'X'
            else:
                return 'X' + str(self._num)
        else:
            return (str(self.number) + self.index).strip()

    def __str__(self):
        return repr(self)
