# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

from skbio.io import read

from . import BaseMod


logger = getLogger(__name__)


class Module(BaseMod):
    def __init__(self, directory, name=__file__):
        super().__init__(directory, name=name)
        self.files = {'faa': self.name + '.faa',
                      'gff': self.name + '.gff'}
        self.ok = self.name + '.ok'

    def parse(self):
        '''Parse the annotation and add it to interval metadata.

        Parameters
        ----------
        fn : str
            the file name from prodigal prediction

        Yield
        -----
        tuple of str and IntervalMetadata
            seq_id and interval metadata
        '''
        self.result = {sid: imd for sid, imd in read(self.files['gff'], format='gff3')}

    def report(self):
        ''''''
