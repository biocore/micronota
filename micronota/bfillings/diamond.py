# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import stat, remove
from os.path import join, basename, splitext, exists
from logging import getLogger
import re

import pandas as pd
from skbio import read, Sequence
from dumpling import (
    check_choice, Dumpling, OptionParam, Parameters)

import string
import random

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


def run(out_dir, query, db, fmt='sam', aligner='blastp', **kwargs):
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
    run_blast(daa=daa, db=db, query=query, aligner=aligner, **kwargs)
    run_view(daa=daa, out=out, fmt=fmt)


class FeatureAnnt:
    '''
    Attributes
    ----------
    dat : list of str
        list of file path to databases.
    cache : list of skbio.Sequence
        list of skbio.Sequences to search against.
    '''
    def __init__(self, dat, out_dir, tmp_dir=None, cache=None):
        super().__init__(dat, out_dir, tmp_dir)
        self.cache = cache
        self.dat = dat

    def _annotate_fp(self, fp, aligner='blastp', evalue=0.001, cpus=1,
                     outfmt='sam', params=None):
        '''Annotate the sequences in the file.'''

        if self.has_cache() and not self.cache.is_empty():
            self.cache.build()
            dbs = [self.cache.db] + self.dat
        else:
            dbs = self.dat

        seqs = []
        found = set()
        res = pd.DataFrame()
        logger = getLogger(__name__)
        for db in dbs:
            out_prefix = splitext(basename(db))[0]
            daa_fp = join(self.out_dir, '%s.daa' % out_prefix)
            out_fp = join(self.out_dir, '%s.diamond' % out_prefix)
            self.run_blast(fp, daa_fp, db, aligner=aligner,
                           evalue=evalue, cpus=cpus, params=params)
            self.run_view(daa_fp, out_fp, params={'--outfmt': outfmt})
            # res = res.append(self.parse_tabular(out_fp))
            if outfmt == 'tab':
                res = res.append(
                    self._filter_best(self.parse_tabular(out_fp)))
            elif outfmt == 'sam':
                res = res.append(
                    self._filter_id_cov(self.parse_sam(out_fp)))

            # save to a tmp file the seqs that do not hit current database
            new_fp = join(self.tmp_dir, '%s.fa' % out_prefix)
            found = found | set(res.index)
            with open(new_fp, 'w') as f:
                for seq in read(fp, format='fasta'):
                    if seq.metadata['id'] not in found:
                        seq.write(f, format='fasta')
            logger.info('Number of diamond hits: %d' % len(res.index))

            # no seq left
            if stat(new_fp).st_size == 0:
                break
            else:
                fp = new_fp
        if outfmt == 'sam' and self.has_cache():
            for x in res.index:
                seqs.append(
                    Sequence(res.loc[x, 'sseq'],
                             metadata={'id': res.loc[x, 'sseqid']}))

        # Update cache (inplace)
        if self.has_cache():
            self.cache.update(seqs)
            self.cache.close()
        return res


def parse_tabular(diamond_res, column='bitscore'):
    '''Parse the tabular output of diamond blastp/blastx.

    Parameters
    ----------
    diamond_res : str
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
    df = pd.read_table(diamond_res, names=columns)
    return df


def parse_sam(diamond_res):
    '''Parse the output of diamond blastp/blastx.

    Parameters
    ----------
    diamond_res : str
        file path

    Returns
    -------
    pandas.DataFrame
        The best matched records for each query sequence.
    '''
    columns = ['qseqid', 'sseqid', 'pident', 'qlen', 'mismatch',
               'qstart', 'sstart', 'evalue', 'bitscore', 'sseq']
    df = pd.DataFrame(columns=columns)
    try:
        seqs = read(diamond_res, format='sam')
    except StopIteration:
        return df
    for i, seq in enumerate(seqs):
        sseq = str(seq)

        qseqid = seq.metadata['QNAME']
        sseqid = seq.metadata['RNAME']
        pident = seq.metadata['ZI']
        qlen = seq.metadata['ZL']
        mismatch = seq.metadata['CIGAR']
        qstart = seq.metadata['POS']
        sstart = seq.metadata['ZS']
        evalue = seq.metadata['ZE']
        bitscore = seq.metadata['ZR']
        row = pd.Series([qseqid, sseqid, pident,
                         qlen, mismatch,
                         qstart, sstart,
                         evalue, bitscore, sseq],
                        index=columns)
        df.loc[i] = row
    return df


def _filter_best(df, column='evalue'):
    # pick the rows that have highest bitscore for each qseqid
    # df_max = df.groupby('qseqid').apply(
    #     lambda r: r[r[column] == r[column].max()])
    if column == 'evalue':
        idx = df.groupby('qseqid')[column].idxmin()
    elif column == 'bitscore':
        idx = df.groupby('qseqid')[column].idxmax()
    df_best = df.loc[idx]
    df_best.set_index('qseqid', drop=True, inplace=True)
    return df_best


def _filter_ident_overlap(df, pident=90, cov=80):
    '''Filter away the hits using the same UniRef clustering standards.'''
    select_id = df.pident >= pident
    aligned_length = df.mismatch.apply(_compute_aligned_length)
    select_cov = ((aligned_length * 100 / df.qlen >= cov) &
                  (aligned_length * 100 / df.sseq.apply(len) >= cov))
    # if qlen * 100 / len(row.sequence) >= 80:
    df_filtered = df[select_id & select_cov]
    df_filtered.set_index('qseqid', drop=True, inplace=True)
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


class DiamondCache:
    '''
    Attributes
    ----------
    out_dir : str
        output directory file path
    fname : str
        fasta file to store cached sequences
    db : str
        diamond database to store cached sequences
    maxSize : int
        maxinum size of DiamondCache
    seqs : list of skbio.Sequence
        list of sequence objects
    '''
    def __init__(self, seqs=None, maxSize=200000, out_dir=""):
        self.out_dir = out_dir
        self.fname = self._generate_random_file()  # substitute for tempfile
        self.fasta = join(out_dir, '%s.fasta' % self.fname)
        self.db = join(out_dir, '%s.dmnd' % self.fname)
        self.maxSize = maxSize
        if seqs is None:
            self.seqs = []
        else:
            self.seqs = seqs

    def _generate_random_file(self, N=10):
        s = ''.join(random.SystemRandom().choice(
                string.ascii_uppercase + string.digits) for _ in range(N))
        # look for unique filepath
        while exists(join(self.out_dir, s)):
            s = ''.join(random.SystemRandom().choice(
                    string.ascii_uppercase + string.digits) for _ in range(N))
        return s

    def dbname(self):
        return self.db.name

    def is_empty(self):
        return (self.seqs is None) or len(self.seqs) == 0

    def build(self, params=None):
        if self.is_empty():
            return
        for seq in self.seqs:
            seq.write(self.fasta, format='fasta')
        make_db(self.fasta, self.db, params)

    def update(self, seqs):
        """
        Parameters
        ----------
        seqs : list of skbio.Sequence
           List of sequences to update the cache.
        """
        self.seqs = seqs + self.seqs
        self.seqs = self.seqs[:self.maxSize]

    def close(self):
        # Remove files if they exist
        # They won't be present if the cache is empty
        if exists(self.fasta):
            remove(self.fasta)
        if exists(self.db):
            remove(self.db)
