import string
import pandas as pd
from skbio.util._misc import merge_dicts
from skbio.io import create_format, FileFormatError
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.io.format._base import (
    _line_generator, _get_nth_sequence, _too_many_blanks)

class SAMFormatError(FileFormatError):
    pass


sam = create_format('sam')
# File headers
_HEADERS=[
          'HD',      # The header line
          'SQ',      # Reference sequence dictionary.
          'RG',      # Read group
          'PG',      # Program
          'CO'      # one-line text comment
          ]

# Alignment headers
_ALIGNMENT_HEADERS=[

          'QNAME',   # Query template NAME.
          'FLAG',    # Combination of bitwise FLAGs
          'RNAME',   # Reference sequence NAME of the alignment
          'POS',     # 1-based leftmost mapping position of the first matching base
          'MAPQ',    # Mapping quality. -10log10(P_err).
          'CIGAR',   # CIGAR string
          'RNEXT',   # Reference sequence name of the primary alignment of the NEXT
          'PNEXT',   # Position of the primary alignment of the NEXT read
          'TLEN',    # signed observed template length
          # 'SEQ',     # segment sequence
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

    try:
        assert line.startswith('@HD')
    except SAMFormatError:
        return False, {}
    return True, {}


def _is_float(input):
  try:
    num = float(input)
  except ValueError:
    return False
  return True

def parse_optional(s):
    seps = s.split(':')
    _field, _type, _val = s.split(':')
    if _type == 'i':
        return {_field: int(_val)}
    elif _type == 'f':
        return {_field: float(_val)}
    else:
        return {_field: _val}

def parse_required(s):
    if s.isdigit():
        return int(s)
    elif _is_float(s):
        return float(s)
    else:
        return s


def _construct(record, constructor=None, **kwargs):
    seq, md = record
    if constructor is None:
        constructor = Sequence
    if constructor == RNA:
        return DNA(
            seq, metadata=md, positional_metadata=pmd, **kwargs).transcribe()
    else:
        return constructor(
            seq, metadata=md, **kwargs)


@sam.reader(None)
def _sam_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh):
        yield _construct(record, constructor, **kwargs)


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
    data_chunks = []
    metadata = {}
    optional_headers = []
    headers = _ALIGNMENT_HEADERS
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        # parse the header (would be nice to abstract this pattern out)
        if line.startswith('@'):
            tabs = line.split('\t')
            key = tabs[0][1:]
            if key == 'CO':
                val = vals[0]
                optional_headers = val.split(',')
                headers = _ALIGNMENT_HEADERS + optional_headers
                headers = headers[:9] + headers[10:]
            else:
                vals = tabs[1:]
                if len(vals) > 1:
                    metadata[key] = vals
                else:
                    metadata[key] = vals[0]

        # parse the actual sequences
        else:
            tabs = line.split('\t')
            # extract sequence
            tabs = tabs[:9] + tabs[10:]
            seq = tabs[9]

            req = list(map(parse_required, tabs[:10]))
            opt = list(map(parse_optional, tabs[10:]))
            req = dict(zip(_ALIGNMENT_HEADERS, req))

            md = merge_dicts(metadata, req, *opt)
            yield seq, md



