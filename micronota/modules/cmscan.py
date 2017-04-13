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
from .modules import BaseMod


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
        splitter = split(SplitterID(lambda s: s.split()[2]),
                         ignore=lambda s: s.startswith('#'))
        with open(self.files['txt']) as fh:
            for lines in splitter(fh):
                k, v =  _parse_record(lines)
                self.result[k] = v


def _parse_record(lines):
    '''Return interval metadata'''
    imd = IntervalMetadata(None)
    seq_id = lines[0].split()[2]
    for line in lines:
        bounds, md = _parse_line(line)
        imd.add(bounds, metadata=md)
    return seq_id, imd


def _parse_line(line):
    items = line.split()
    fam_id = items[1]
    md = {'source': 'Rfam', 'ncRNA_class': items[0], 'db_xref': fam_id}
    if fam_id in {'RF00001', 'RF00177', 'RF02541', 'RF01959',
                  'RF02540', 'RF00001', 'RF00002', 'RF01960', 'RF02543'}:
        md['type'] = 'rRNA'
        if fam_id == 'RF00001':
            md['product'] = '5s_rRNA'
        elif fam_id in {'RF00177', 'RF01959'}:
            md['product'] = '16s_rRNA'
        elif fam_id in {'RF02541', 'RF02540'}:
            md['product'] = '23s_rRNA'
    else:
        md['type'] = 'ncRNA'

    strand = items[9]
    md['strand'] = strand
    if strand == '+':
        start = int(items[7]) - 1
        end = int(items[8])
    elif strand == '-':
        start = int(items[8]) - 1
        end = int(items[7])
    else:
        raise ValueError('Unknown strand for the ncRNA: %s' % line)
    return [(start, end)], md
