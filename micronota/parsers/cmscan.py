# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from logging import getLogger

from skbio.io.format._sequence_feature_vocabulary import _yield_section
from skbio.metadata import IntervalMetadata


logger = getLogger(__name__)


def parse(out_dir, fn='cmscan.txt'):
    '''Parse the annotation and add it to interval metadata.

    Parameters
    ----------
    out_dir : str
        the output dir
    fn : str
        the file name from prediction

    Yield
    -----
    tuple of str and IntervalMetadata
        seq_id and interval metadata
    '''
    logger.debug('Parsing cmscan prediction of Rfam')
    fp = join(out_dir, fn)
    return _cmscan_to_interval_metadata(fp)


def _cmscan_to_interval_metadata(fp):
    '''Yield seq_id and its interval metadata.'''
    current = False
    lines = []
    with open(fp) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            items = line.split()
            seq_id = items[2]
            if current == seq_id:
                lines.append(line)
            else:
                if current is not False:
                    yield current, _parse_record(lines)
                current = seq_id
        yield current, _parse_record(lines)


def _parse_record(lines):
    '''Return interval metadata'''
    imd = IntervalMetadata(None)
    for line in lines:
        bounds, md = _parse_line(line)
        imd.add(bounds, metadata=md)
    return imd


def _parse_line(line):
    md = {'type': 'ncRNA'}
    items = line.split()
    md['db_xref'] = items[1]
    md['ncRNA_class'] = items[0]
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
