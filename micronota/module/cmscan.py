# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

from skbio.metadata import IntervalMetadata

from . import BaseMod


logger = getLogger(__name__)


class Module(BaseMod):
    def __init__(self, directory, name=__file__):
        super().__init__(directory, name=name)
        self.files = {'txt': self.name + '.txt'}
        self.ok = self.name + '.ok'

    def parse(self):
        '''Parse the annotation and add it to interval metadata.

        Parameters
        ----------
        fp : str
            the file path from cmscan run.

        Yield
        -----
        tuple of str and IntervalMetadata
            seq_id and interval metadata
        '''
        for seqid, imd in read(self.files['txt'], format='cmscan'):
            self.result[seqid] = imd
