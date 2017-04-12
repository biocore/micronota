r'''
Utility functionality
=====================

.. currentmodule:: micronota.util

This module (:mod:`micronota.util`) provides various utility functionality,
'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase
from sqlite3 import connect

from skbio import read, write, Sequence


def convert(in_f, in_fmt, out_f, out_fmt):
    '''convert between file formats

    Parameters
    ----------
    in_fmt : str
        input file format
    out_fmt : str
        output file format
    in_f : str
        input file path
    out_f: str
        output file path
    '''
    for obj in read(in_f, format=in_fmt):
        write(obj, format=out_fmt, into=out_f)


def _filter_sequence_ids(in_fp, out_fp, ids, negate=False):
    '''Filter away the seq with specified IDs.'''
    with open(out_fp, 'w') as out:
        for seq in read(in_fp, format='fasta', constructor=Sequence):
            seq_id = seq.metadata['id']
            if seq_id not in ids:
                write(seq, format='fasta', into=out)

    def filter_ident_overlap(df, pident=90, overlap=80):
        '''Filter away the hits using the same UniRef clustering standards.

        Parameters
        ----------
        df : ``pandas.DataFrame``
            it must have columns of 'pident', 'gaps', and 'slen'
        pident : ``Numeric``
            minimal percentage of identity
        overlap : ``Numeric``
            minimal percentage of overlap for subject sequences.

        Returns
        -------
        ``pandas.DataFrame``
            The data frame only containing hits that pass the thresholds.
        '''
        select_id = df.pident >= pident
        overlap_length = df.length - df.gaps
        select_overlap = overlap_length * 100 / df.slen >= overlap
        # if qlen * 100 / len(row.sequence) >= 80:
        df_filtered = df[select_id & select_overlap]
        # df_filtered.set_index('qseqid', drop=True, inplace=True)
        return df_filtered


class _DBTest(TestCase):
    def _test_eq_db(self, db1, db2):
        '''Test if two database files have the same contents.'''
        with connect(db1) as o, connect(db2) as e:
            co = o.cursor()
            ce = e.cursor()
            # compare each table
            for table, in e.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"):
                co.execute('SELECT * FROM %s' % table)
                ce.execute('SELECT * FROM %s' % table)
                self.assertEqual(co.fetchall(), ce.fetchall())


def split(is_another, construct=None, ignore=None, **kwargs):
    '''Return a function to yield record.

    Parameters
    ----------
    is_another : callable
        accept a string and return a bool indicating if the current
        line starts a new record or not.
    construct : callable (optional)
        accept a string and return a modified each input line. Default
        is to strip the whitespaces at the ends. Do nothing if it is
        ``None``.
    ignore : callable (optional)
        A callable to ignore the line if it returns ``True``. Do nothing
        if it is ``None``.
    kwargs : dict
        optional key word arguments passing to ``is_another``

    Returns
    -------
    function
        a function that accepts a file-object-like and yields record
        one by one from it.

    '''
    def parse(stream):
        lines = []
        for line in stream:
            if ignore is not None and ignore(line):
                continue
            if is_another(line, **kwargs):
                if lines:
                    yield lines
                    lines = []
            if construct is not None:
                line = construct(line)
            lines.append(line)
        if lines:
            yield lines
    return parse


def split_head(line, is_head=lambda line: line.startswith('>')):
    '''

    Examples
    --------
    >>> import io
    >>> s = """>seq1
    ... ATGC
    ... >seq2
    ... A
    ... T
    ... """
    >>> f = io.StringIO(s)
    >>> gen = split(split_head, construct=lambda x: x.strip())
    >>> list(gen(f))
    [['>seq1', 'ATGC'], ['>seq2', 'A', 'T']]
    '''
    if is_head(line):
        return True
    else:
        return False


class SplitterTail:
    r'''Split a stream of lines into record delimited by its tail of each record.

    Parameters
    ----------
    is_tail : callable
        to return a bool indicating if the current line concludes current record

    Examples
    --------
    Let's create a file object that has sequences separated by "//" at
    the end of each record (similar to multiple GenBank records in a file):

    >>> import io
    >>> s = """seq1
    ... AT
    ... //
    ... seq2
    ... ATGC
    ... //
    ... """
    >>> f = io.StringIO(s)

    And then we can yield each sequence with this code:

    >>> splitter = SplitterTail(lambda x: x == '//\n')
    >>> gen = split(splitter, construct=lambda x: x.strip())
    >>> list(gen(f))
    [['seq1', 'AT', '//'], ['seq2', 'ATGC', '//']]
    '''
    def __init__(self, is_tail):
        self.is_tail = is_tail
        self._flag = False

    def __call__(self, line):
        if self.is_tail(line):
            self._flag = True
            return False
        else:
            if self._flag:
                self._flag = False
                return True


class SplitterID:
    r'''Split a stream of lines into record based on every line of each record.

    Parameters
    ----------
    identify : callable
        accept a string and return a value. By comparing the value to
        that of its previous line, it judges if the current line is
        like its previous line belonging to the same record

    Examples
    --------
    >>> from pprint import pprint
    >>> import io
    >>> s = """##gff-version 3
    ... ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001
    ... ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002
    ... ctg124\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003
    ... """
    >>> f = io.StringIO(s)
    >>> splitter = SplitterID(lambda x: x.split('\t')[0])
    >>> gen = split(splitter, lambda x: x.strip(), lambda x: x.startswith('#'))
    >>> pprint(list(gen(f)))
    [['ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001',
      'ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002'],
     ['ctg124\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003']]

    '''
    def __init__(self, identify):
        self._ident = None
        self.identify = identify

    def __call__(self, line):
        ident = self.identify(line)
        if self._ident == ident:
            return False
        else:
            self._ident = ident
            return True
