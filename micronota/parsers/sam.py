r'''
SAM Parser
==========

SAM stores input query sequences that have mapped to a reference [#]_.
The SAM format has a header section, where each line is denoted by @.

Each line following the header lines encodes information for a single read.
Each line is structured as follows.

    +--------------------------------+----------------------------------+
    | QNAME - Query template NAME    | (begins each entry; 1 per entry) |
    +--------------------------------+----------------------------------+
    | FLAG  - Bitwise FLAGs          | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | RNAME - Reference sequence     | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | POS - 1-based position         | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | MAPQ - mapping quality         | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | CIGAR - CIGAR string           | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | RNEXT - next alignment name    | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | PNEXT - next alignment position| 1 per entry                      |
    +--------------------------------+----------------------------------+
    | TLEN - signed observed length  | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | SEQ - segment sequence         | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | QUAL - base quality            | 1 per entry                      |
    +--------------------------------+----------------------------------+

The line can be followed by optional fields.


Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+
|Yes   |No    |generator of :mod:`skbio.sequence.Sequence` objects            |
+------+------+---------------------------------------------------------------+


Reference
---------
.. [#] https://samtools.github.io/hts-specs/SAMv1.pdf
'''

from skbio.io import create_format
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.io.format._base import (
    _line_generator, _get_nth_sequence, _too_many_blanks)


sam = create_format('sam')

# Alignment headers
_REQUIRED_FIELDS = [
          'QNAME',   # Query template NAME.
          'FLAG',    # Combination of bitwise FLAGs
          'RNAME',   # Reference sequence NAME of the alignment
          'POS',     # 1-based leftmost mapping position of the first base
          'MAPQ',    # Mapping quality. -10log10(P_err).
          'CIGAR',   # CIGAR string
          'RNEXT',   # Reference sequence name of the primary alignment of NEXT
          'PNEXT',   # Position of the primary alignment of the NEXT read
          'TLEN',    # signed observed template length
          'SEQ',     # segment sequence
          'QUAL',    # ASCII of base quality
]


@sam.sniffer()
def _sam_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}

    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    if line.startswith('@HD'):
        return True, {}
    else:
        return False, {}


def _construct(record, constructor=None, **kwargs):
    seq, md = record
    if constructor is None:
        constructor = Sequence
    if constructor == RNA:
        return DNA(
            seq, metadata=md, **kwargs).transcribe()
    else:
        return constructor(
            seq, metadata=md, **kwargs)


@sam.reader(None)
def _sam_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh):
        yield from _construct(record, constructor, **kwargs)


@sam.reader(Sequence)
def _sam_to_sequence(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, Protein, **kwargs)


@sam.reader(Protein)
def _sam_to_protein(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, Protein, **kwargs)


@sam.reader(DNA)
def _sam_to_DNA(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, DNA, **kwargs)


@sam.reader(RNA)
def _sam_to_RNA(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh), seq_num)
    return _construct(record, RNA, **kwargs)


def _parse_records(fh, constructor=None, **kwargs):
    md = {}
    res = None
    opt_fields = {}
    n = len(_REQUIRED_FIELDS)
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        # parse the header (would be nice to abstract this pattern out)
        if line.startswith('@'):
            key, val = line.split('\t', 1)
            if key == '@CO' and val.startswith('Reporting'):
                _, s = val.split(None, 1)
                opt_fields = _parse_co(s)
            else:
                md[key] = val
        # parse the actual sequences
        else:
            tabs = line.split('\t')
            # zip stops generating after the shorter list of the two
            md = dict(zip(_REQUIRED_FIELDS, tabs))
            seq = md.pop('SEQ')

            opt = (_parse_optional(field) for field in tabs[n:])
            opt = {opt_fields.get(k, k): v for k, v in opt}
            md.update(opt)

            res = seq, md
            yield res
    # this is to fix the skbio's surprising behavior of read into generator
    if res is None:
        return iter([])


def _parse_optional(s):
    field, t, v = s.split(':')
    if t == 'i':
        v = int(v)
    elif t == 'f':
        v = float(v)
    return field, v


def _parse_required(s):
    if s.isdigit():
        return int(s)
    else:
        try:
            return float(s)
        except ValueError:
            return s


def _parse_co(s):
    '''Parse CO string.

    Examples
    --------
    >>> s = 'AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate'
    >>> _parse_co(s)
    {'AS': 'bitScore',
     'ZE': 'expected',
     'ZF': 'frame',
     'ZI': 'percent identity',
     'ZL': 'reference length',
     'ZR': 'rawScore',
     'ZS': 'query start DNA coordinate'}
    '''
    items = s.split(', ')
    fields = (i.split(': ', 1) for i in items)
    return dict(fields)
