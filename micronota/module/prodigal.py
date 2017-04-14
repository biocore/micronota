# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re

from skbio.io import read

from . import BaseMod


class Module(BaseMod):
    def __init__(self, directory, file_patterns=None):
        if file_patterns is None:
            file_patterns = {'gff': 'prodigal.gff'}
        super().__init__(directory, file_patterns=file_patterns)

    def parse(self):
        '''Parse the annotation and add it to interval metadata.'''
        self.result = {sid: imd for sid, imd in read(self.files['gff'], format='gff3')}

    def report(self):
        ''''''
