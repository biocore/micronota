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


logger = getLogger(__name__)


def parse(fp='cmscan.txt'):
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
    logger.debug('Parsing cmscan prediction of Rfam')
    splitter = split(SplitterID(lambda s: s.split()[2]),
                     ignore=lambda s: s.startswith('#'))
    with open(fp) as fh:
        for lines in splitter(fh):
            yield _parse_record(lines)


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
