# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio import read

from . import BaseMod


class Module(BaseMod):
    def __init__(self, directory, file_patterns=None):
        if file_patterns is None:
            file_patterns = {'txt': 'tandem_repeats_finder.txt'}
        super().__init__(directory, file_patterns)

    def parse(self):
        '''Parse the annotation and add it to interval metadata. '''
        for seqid, imd in read(self.files['txt'], format='tandem_repeats_finder'):
            self.result[seqid] = imd
