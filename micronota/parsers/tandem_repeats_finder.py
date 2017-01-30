# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from logging import getLogger

from skbio.metadata import IntervalMetadata

from ..util import split, split_head


logger = getLogger(__name__)


def parse(fp='tandem_repeats_finder.txt'):
    '''Parse the annotation and add it to interval metadata.

    Parameters
    ----------
    fp : str
        the file path from Tandem Repeat Finder prediction

    Yield
    -----
    tuple of str and IntervalMetadata
        seq_id and interval metadata
    '''
    logger.debug('Parsing tandem repeat prediction')

    splitter = split(split_head, is_head=lambda line: line.startswith('@'))
    with open(fp) as fh:
        for lines in splitter(fh):
            sid = lines[0].split(None, 1)[0][1:]
            yield sid, _parse_record(lines[1:])


def _parse_record(lines):
    '''Return interval metadata.'''
    imd = IntervalMetadata(None)
    for line in lines:
        bounds, md = _parse_line(line)
        imd.add(bounds, metadata=md)
    return imd


def _parse_line(line):
    # the table has following columns:
    # Indices of the repeat relative to the start of the sequence.
    # Period size of the repeat.
    # Number of copies aligned with the consensus pattern.
    # Size of consensus pattern (may differ slightly from the period size).
    # Percent of matches between adjacent copies overall.
    # Percent of indels between adjacent copies overall.
    # Alignment score.
    # Percent composition for each of the four nucleotides.
    # Entropy measure based on percent composition.
    # the repeat
    # the actual sequence has all the repeats
    # left and right flanking sequences
    items = line.split(' ')
    md = {'repeat': items[13], 'type': 'tandem_repeat', 'source': 'Tandem_Repeats_Finder'}
    start = int(items[0]) - 1
    end = int(items[1])
    return [(start, end)], md
