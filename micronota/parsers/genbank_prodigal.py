r'''
GenBank variant for Prodigal
============================

.. currentmodule:: micronota.parsers.genbankprodigal

This module parse the Prodigal output of GenBank format.

Format Support
--------------
+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`skbio.sequence.IntervalMetadata`                         |
+------+------+---------------------------------------------------------------+

'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from skbio.io import create_format
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.io.format._base import (
    _line_generator, _get_nth_sequence, _too_many_blanks)


genbankprodigal = create_format('genbankprodigal')

# This list is ordered. From EMBL specification
_HEADERS = ['ID',  # identification            (begins each entry; 1 per entry)
            'AC',  # accession number          (>=1 per entry)
            'PR',  # project identifier        (0 or 1 per entry)
            'DT',  # date                      (2 per entry)
            'DE',  # description               (>=1 per entry)
            'KW',  # keyword                   (>=1 per entry)
            'OS',  # organism species          (>=1 per entry)
            'OC',  # organism classification   (>=1 per entry)
            'OG',  # organelle                 (0 or 1 per entry)
            'RN',  # reference number          (>=1 per entry)
            'RC',  # reference comment         (>=0 per entry)
            'RP',  # reference positions       (>=1 per entry)
            'RX',  # reference cross-reference (>=0 per entry)
            'RG',  # reference group           (>=0 per entry)
            'RA',  # reference author(s)       (>=0 per entry)
            'RT',  # reference title           (>=1 per entry)
            'RL',  # reference location        (>=1 per entry)
            'DR',  # database cross-reference  (>=0 per entry)
            'CC',  # comments or notes         (>=0 per entry)
            'AH',  # assembly header           (0 or 1 per entry)
            'AS',  # assembly information      (0 or >=1 per entry)
            'FH',  # feature table header      (2 per entry)
            'FT',  # feature table data        (>=2 per entry)
            'SQ',  # sequence header           (1 per entry)
            'CO',  # contig/construct line     (0 or >=1 per entry)
            'bb']  # (blanks) sequence data    (>=1 per entry)


@embl.sniffer()
def _embl_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}
    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    try:
        _parse_id([line])
    except:
        return False, {}
    return True, {}


@embl.reader(None)
def _embl_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh, _parse_single_embl):
        yield _construct(record, constructor, **kwargs)


@embl.reader(Sequence)
def _embl_to_sequence(fh, seq_num=1, **kwargs):
    pass


@embl.reader(Protein)
def _embl_to_protein(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh, _parse_single_embl), seq_num)
    return _construct(record, Protein, **kwargs)


@embl.reader(DNA)
def _embl_to_DNA(fh, seq_num=1, **kwargs):
    pass


@embl.reader(RNA)
def _embl_to_RNA(fh, seq_num=1, **kwargs):
    pass


def _construct(record, constructor=None, **kwargs):
    seq, md, pmd = record
    if constructor is None:
        unit = md['ID']['unit'].lower()
        if unit == 'bp':
            constructor = DNA
        elif unit == 'aa':
            constructor = Protein
    if constructor == RNA:
        return DNA(
            seq, metadata=md, positional_metadata=pmd, **kwargs).transcribe()
    else:
        return constructor(
            seq, metadata=md, positional_metadata=pmd, **kwargs)


def _parse_records(fh, parser):
    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True, strip=False):
        if line.startswith('//'):
            yield parser(data_chunks)
            data_chunks = []
        else:
            data_chunks.append(line)


def _parse_single_embl(chunks):
    metadata = {}
    positional_metadata = None
    sequence = ''

    def split(s):
        if s[0].isspace():
            return None, s.strip()
        else:
            return s.split(None, 1)

    for header, data in _yield_section(chunks, split, strip=False):
        if header in ['XX', 'RN', 'RC', 'RP', 'RX', 'RG', 'RA', 'RT', 'RL']:
            continue
        parser = _PARSER_TABLE.get(header)
        if parser:
            data = parser(data)
        if header == 'SQ':
            sequence = data
        else:
            metadata[header] = data

    return sequence, metadata, positional_metadata


def _parse_id(lines, strict=False):
    '''Parse ID line.

    Parameters
    ----------
    strict : bool.
        Whether to parse strictly compliant to EMBL specifications.'''
    # strip the ending dot
    if not strict:
        items = lines[0].rstrip('.').split(';')
        res = dict()
        res['id'] = items[0].split()[0]
        res['size'], res['unit'] = items[-1].split()

        return res


def _parse_ac(lines):
    return lines[0].rstrip(';')


def _parse_oc(lines):
    return ' '.join(lines)


def _parse_ox(lines):
    line = lines[0]
    taxon_id = line.replace('NCBI_TaxID=', '').rstrip(';')
    return taxon_id


def _parse_oh(lines):
    '''Parse the OH (Organism Host) line.

    Optional and appears only in viral entries. It indicates the host
    organism(s) that are susceptible to be infected by a virus.

    Example line:
    'NCBI_TaxID=8295; Ambystoma (mole salamanders).'
    '''
    taxa = dict()
    for line in lines:
        taxon_id, taxon = line.split(';')
        taxon_id = taxon_id.replace('NCBI_TaxID=', '')
        taxon = taxon.strip(' .')
        taxa[taxon_id] = taxon
    return taxa


def _parse_dr(lines):
    '''Parse the DR (Database Cross-reference) line.'''
    res = dict()
    for db_name, db_lines in _yield_section(lines, lambda s: s.split('; ', 1)):
        res[db_name] = db_lines
    return res


def _parse_pe(lines):
    return lines[0].split(':')[0]


def _parse_kw(lines):
    l = ' '.join(lines).rstrip('.')
    return l.split('; ')


def _parse_sq(lines):
    return ''.join(''.join(lines[1:]).split())


def _yield_section(lines, split_header, **kwargs):
    '''Yield the lines with the same header.

    Parameters
    ----------
    split_header : function
        It accepts a string of line and returns the header and
        rest of the needed data as a tuple. If no header exists,
        return ``None`` and the data.

    Notes
    -----
    The following example is a valid section::

    DR   GO GO:000001
    DR   GO GO:000002

    This is also a section::

    SQ   ATGCA ATGCA
         ATGCA ATGCA
    '''
    curr = []
    header, _ = split_header(lines[0])

    for line in _line_generator(lines, **kwargs):
        items = split_header(line)

        # if the header is changed, it is a new section
        if items[0] is not None and items[0] != header:
            if curr:
                yield header, curr
                curr = []
                header = items[0]

        curr.append(items[1].strip())
    # don't forget to return the last section in the file
    if curr:
        yield header, curr


_PARSER_TABLE = {
    'ID': _parse_id,
    'AC': _parse_ac,
    'OC': _parse_oc,
    'OX': _parse_ox,
    'OH': _parse_oh,
    'KW': _parse_kw,
    'DR': _parse_dr,
    'PE': _parse_pe,
    'SQ': _parse_sq}
