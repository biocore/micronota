# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from abc import ABC, abstractmethod
from os.path import basename, splitext, join
from os import listdir
import re


class BaseMod(ABC):
    def __init__(self, directory, file_patterns):
        self.result = {}
        self.report = {}
        self.files = {}

        for k in file_patterns:
            matched_f = [f for f in listdir(directory) if re.match(file_patterns[k], f)]
            if len(matched_f) != 1:
                raise ValueError('There are multiple files matching %s' % file_patterns[k])
            self.files[k] = join(directory, matched_f[0])

    @abstractmethod
    def parse(self):
        '''parse result'''
