# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from logging import getLogger
import re

from skbio.io.format._sequence_feature_vocabulary import _yield_section
from skbio.metadata import IntervalMetadata


logger = getLogger(__name__)


def parse(out_dir, fn='aragorn.txt'):
    '''Parse the annotation and add it to interval metadata.

    Parameters
    ----------
    out_dir : str
        the dir where the aragorn output file exist
    fn : str
        the file name from aragorn prediction

    Yield
    -----
    tuple of str and IntervalMetadata
        seq_id and interval metadata
    '''
    logger.debug('Parsing aragorn prediction')
    fp = join(out_dir, fn)
    return _aragorn_to_interval_metadata(fp)


def _aragorn_to_interval_metadata(fp):
    '''Yield seq_id and its interval metadata.'''
    # aragorn output has a final summary line like this:
    # >end    5 sequences 97 tRNA genes 1 tmRNA genes
    # This line should be skipped and not parsed
    p = re.compile(r'>end\s+\d+ sequences \d+ tRNA genes \d+ tmRNA genes')
    splitter = _yield_section(lambda line: line.startswith('>'))
    with open(fp) as fh:
        for lines in splitter(fh):
            headline = lines[0]
            if p.match(headline):
                return
            sid = headline.split(None, 1)[0][1:]
            # the first 2 lines are not actual data lines
            yield sid, _parse_record(lines[2:])


def _parse_record(lines):
    '''Return interval metadata.'''
    imd = IntervalMetadata(None)
    for line in lines:
        bounds, md = _parse_line(line)
        imd.add(bounds, metadata=md)
    return imd


def _parse_line(line):
    _, tRNA, loc, _, _ = line.split()
    md = {'strand': '+', 'product': tRNA, 'type': 'tRNA', 'source': 'Aragorn'}
    if loc[0] == 'c':
        loc = loc[1:]
        md['strand'] = '-'
    start, end = loc.split(',')
    start = int(start[1:]) - 1
    end = int(end[:-1])
    return [(start, end)], md
