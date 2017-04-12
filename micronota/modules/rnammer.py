# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

from skbio.metadata import IntervalMetadata

from ..util import split, SplitterID
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
            the file name from RNAmmer prediction

        Yield
        -----
        tuple of str and IntervalMetadata
            seq_id and interval metadata
        '''
        logger.debug('Parsing RNAmmer prediction')
        splitter = split(SplitterID(lambda s: s.split('\t')[0]),
                         construct=lambda s: s.strip(),
                         ignore=lambda s: s.startswith('#'))
        with open(self.files['gff']) as fh:
            for lines in splitter(fh):
                k, v = _parse_record(lines)
                self.result[k] = v


def _parse_record(lines):
    imd = IntervalMetadata(None)
    seq_id = lines[0].split('\t')[0]
    for line in lines:
        bounds, md = _parse_line(line)
        imd.add(bounds, metadata=md)
    return seq_id, imd


def _parse_line(line):
    md = {}
    items = line.split('\t')
    md['source'] = items[1]
    md['type'] = items[2]
    md['score'] = items[5]
    md['strand'] = items[6]
    md['product'] = items[-1]
    start = int(items[3]) - 1
    end = int(items[4])
    return [(start, end)], md
