# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

import pandas as pd


logger = getLogger(__name__)


def parse_diamond_uniref(fn, pident=90):
    ''''''
    def parse(columns=['qseqid', 'qlen', 'sseqid', 'slen',
                       'pident', 'length', 'gaps', 'evalue', 'bitscore',
                       'qstart', 'qend', 'sstart', 'send']):
        logger.debug('Parsing diamond search against UniRef')
        df = pd.read_table(fn, names=columns)
        df_filtered = filter_ident_overlap(df, pident)
        return {k: v for k, v in zip(df_filtered.qseqid, df_filtered.sseqid)}
    return parse


parsers = {'parse_diamond_uniref90':
           parse_diamond_uniref(pident=90, fn='diamond_uniref90.m12'),
           'parse_diamond_uniref50':
           parse_diamond_uniref(pident=50, fn='diamond_uniref50.m12')}


def filter_ident_overlap(df, pident=90, overlap=80):
    '''Filter away the hits using the same UniRef clustering standards.

    Parameters
    ----------
    df : ``pandas.DataFrame``
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
