# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, basename, splitext
from logging import getLogger
import re

import pandas as pd

from dumpling import (
    check_choice, Dumpling, OptionParam, Parameters)


blast_params = [
    OptionParam('--threads', 'cpus', help='number of cpu threads.'),

    OptionParam('--gapopen', help='Gap open penalty.'),
    OptionParam('--gapextend', help='Gap extension penalty.'),
    OptionParam('--matrix', help='Scoring matrix.'),
    OptionParam('--seg', help='Enable SEG masking.'),
    OptionParam('--max-target-seqs', '-k',
                help='The maximum number of hits per query to keep alignments for.'),
    OptionParam('--top',
                help='Keep alignments within the given percentage range of the top alignment'),
    OptionParam('--evalue', '-e', help='Maximum expected value to keep an alignment.'),
    OptionParam('--min-score',
                help=('Minimum bit score to keep an alignment. Setting this option'
                      'will override the --evalue parameter.')),
    OptionParam('--query-cover',
                help='Report only alignments above the given percentage of query cover.'),
    OptionParam('--salltitles',
                help='Print full length subject titles in output.'),

    OptionParam('--band', help=''),
    OptionParam('--index-chunks', '-c',
                help='The number of chunks for processing the seed index.'),

    OptionParam('--sensitive',
                help=('Trigger the sensitive alignment mode with a 16x9 seed'
                      'shape config.')),
    OptionParam('--tmpdir', '-t',
                help='Directory to be used for temporary storage.'),

    OptionParam('--db', '-d', help='Path to DIAMOND database file'),
    OptionParam('--daa', '-a',
                help='Path to DAA file.'),
    OptionParam('--query', '-q',
                help=('Path to query input file in FASTA or FASTQ format '
                      '(may be gzip compressed).'))]


makedb_params = [
    OptionParam('--threads', '-p', name='cpus'),
    OptionParam('--in', name='fasta',
                help='protein reference database file in FASTA format (may be gzip compressed)'),
    OptionParam('--db', '-d',
                help='DIAMOND database file.'),
    OptionParam('--block-size',
                help='Block size in billions of sequence letters to be processed at a time.')]


view_params = [
    OptionParam('--outfmt', name='fmt', action=check_choice(('tab', 'sam')),
                help=('Format of output file. (tab = BLAST tabular format;'
                      'sam = SAM format)')),
    OptionParam('--daa', '-a',
                help='Path to DAA file.'),
    OptionParam('--out', '-o',
                help='Path to output file.'),
    OptionParam('--compress', action=check_choice((0, 1)),
                help='Compression for output file (0=none, 1=gzip).')]


def run_blast(query, daa, aligner='blastp', **kwargs):
    '''Search query sequences against the database.

    Parameters
    ----------
    query : str
        The file path of the query seq
    daa : str
        The file path of the output daa file
    aligner : str
        The aligner. blastp or blastx
    kwargs : dict
        keyword arguments. Command line parameters for diamond blastp
        or blastx.
    Returns
    -------
    str
        The file path of the blast result.
    '''
    logger = getLogger(__name__)

    blast = Dumpling(['diamond', aligner],
                     params=Parameters(*blast_params))
    blast.update(query=query, daa=daa, **kwargs)
    logger.info('Running {}'.format(blast.command))
    blast()
    return blast


def run_view(daa, out, fmt='sam', **kwargs):
    '''
    Parameters
    ----------
    daa : str
        Input file resulting from diamond blast.
    out : str
        Output file.
    '''
    logger = getLogger(__name__)
    view = Dumpling(['diamond', 'view'],
                    params=Parameters(*view_params))
    view.update(daa=daa, out=out, fmt=fmt, **kwargs)
    logger.info('Running {}'.format(view.command))
    view()
    return view


def run_makedb(fasta, db=None, **kwargs):
    '''Format database from a fasta file.

    This is similar to running ``diamond makedb --in db.faa --db db``.

    Parameters
    ----------
    fasta : str
        Input path for the fasta file.
    db : str or None (default)
        Output path for the formatted database file. It will be named
        after input file in the same directory by default.
    kwargs : dict
        keyword arguments. Other command line parameters for diamond makedb.

    Returns
    -------
    `Dumpling`
    '''
    logger = getLogger(__name__)
    if db is None:
        db = splitext(fasta)[0]
    makedb = Dumpling(['diamond', 'makedb'], params=Parameters(*makedb_params),
                      version='0.7.12', url='https://github.com/bbuchfink/diamond')
    makedb.update(fasta=fasta, db=db, **kwargs)
    logger.info('Running {}'.format(makedb.command))
    makedb()
    return makedb


def run(out_dir, query, db, query_cover=0, fmt='sam', aligner='blastp', **kwargs):
    '''
    Run Diamond search.

    Parameters
    ----------
    out_dir : str
        output dir
    query : str
        file path to query sequence
    db : str
        file path the db
    query_cover : `Numeric`
        report hits above the given percentage of query cover.
    fmt : str ('sam' or 'tab')
        output file format
    aligner : str ('blastp' or 'blastx'
        which aligning mode to use
    kwargs : dict
        keyword arguments passing to `run_blast`

    Returns
    -------
    None
    '''
    logger = getLogger(__name__)
    prefix = splitext(basename(db))[0]
    daa = join(out_dir, '{}.daa'.format(prefix))
    out = join(out_dir, '{0}.{1}'.format(prefix, fmt))
    logger.info('Running Diamond search ...')
    run_blast(daa=daa, db=db, query=query, aligner=aligner, query_cover=query_cover, **kwargs)
    run_view(daa=daa, out=out, fmt=fmt)


def parse_tabular(res, column='bitscore'):
    '''Parse the tabular output of diamond blastp/blastx.

    Parameters
    ----------
    res : str
        file path to Diamond tabular output
    column : str
        The column used to pick the best hits.

    Returns
    -------
    pandas.DataFrame
        The hit records for each query sequence.
    '''
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
               'gapopen', 'qstart', 'qend', 'sstart', 'send',
               'evalue', 'bitscore']
    df = pd.read_table(res, names=columns)
    return df


def parse_sam(res):
    '''Parse the SAM output of diamond blastp/blastx.

    Parameters
    ----------
    res : str
        file path to Diamond SAM output

    Returns
    -------
    pandas.DataFrame
        The hit records for each query sequence.
    '''
    columns = [
        'qseqid',   # Query template NAME.
        'FLAG',    # Combination of bitwise FLAGs
        'sseqid',   # Reference sequence NAME of the alignment
        'POS',     # 1-based leftmost mapping position of the first base
        'MAPQ',    # Mapping quality. -10log10(P_err).
        'CIGAR',   # CIGAR string
        'RNEXT',   # Reference sequence name of the primary alignment of NEXT
        'PNEXT',   # Position of the primary alignment of the NEXT read
        'TLEN',    # signed observed template length
        'SEQ',     # segment sequence
        'QUAL']    # ASCII of base quality
    optional = [
        'bitscore',
        'NM',       # Edit distance to the reference
        'slen',     # subject seq length
        'rawscore',
        'evalue',
        'pident',
        'frame',
        'qstart',   # start position of alignment
        'MD']

    df = pd.read_table(res, names=columns + optional, comment='@')

    for col in optional:
        df[col] = df[col].apply(_convert)

    return df


def filter_best(df, column='evalue'):
    '''Filter out the best hits by their e-value or bitscore.'''
    # pick the rows that have highest bitscore for each qseqid
    # df_max = df.groupby('qseqid').apply(
    #     lambda r: r[r[column] == r[column].max()])
    if column == 'evalue':
        idx = df.groupby('qseqid')[column].idxmin()
    elif column == 'bitscore':
        idx = df.groupby('qseqid')[column].idxmax()
    df_best = df.loc[idx]
    # df_best.set_index('qseqid', drop=True, inplace=True)
    return df_best


def filter_ident_overlap(df, pident=90, overlap=80):
    '''Filter away the hits using the same UniRef clustering standards.

    Parameters
    ----------
    df : `pandas.DataFrame`
        parsed from `parse_sam`
    pident : `Numeric`
        minimal percentage of identity
    overlap : `Numeric`
        minimal percentage of overlap for subject sequences.

    Returns
    -------
    `pandas.DataFrame`
        The data frame only containing hits that pass the thresholds.
    '''
    select_id = df.pident >= pident
    overlap_length = df.CIGAR.apply(_compute_aligned_length)
    select_overlap = overlap_length * 100 / df.slen >= overlap
    # if qlen * 100 / len(row.sequence) >= 80:
    df_filtered = df[select_id & select_overlap]
    # df_filtered.set_index('qseqid', drop=True, inplace=True)
    return df_filtered


def _compute_aligned_length(cigar):
    '''Compute the length of ungapped region in the alignment.

    It includes both matched and mismatched regions.

    Examples
    --------
    >>> [_compute_aligned_length(i) for i in (
    ...     '',
    ...     '45D',
    ...     '18M2D19M')]
    [0, 0, 37]

    '''
    aligned = re.findall('([0-9]+)M', cigar)
    return sum(int(i) for i in aligned)


def _convert(s):
    field, t, v = s.split(':')
    if t == 'i':
        v = int(v)
    elif t == 'f':
        v = float(v)
    return v
