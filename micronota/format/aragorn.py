import re

from skbio.io import create_format
from skbio.metadata import IntervalMetadata

from ..util import split, split_head


aragorn = create_format('aragorn')


@aragorn.reader(None)
def _generator(fh):
    # aragorn output has a final summary line like this:
    # >end    5 sequences 97 tRNA genes 1 tmRNA genes
    # This line should be skipped and not parsed
    p = re.compile(r'>end\s+\d+ sequences \d+ tRNA genes \d+ tmRNA genes')
    splitter = split(split_head)
    for lines in splitter(fh):
        headline = lines[0]
        if p.match(headline):
            return
        sid = headline.split(None, 1)[0][1:]
        yield sid,  _parse_record(lines[2:])


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
    start = int(start[1:])
    end = int(end[:-1])
    if start > end:
        start, end = end, start
    start -= 1
    return [(start, end)], md


