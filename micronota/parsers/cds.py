# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger
from collections import defaultdict

import pandas as pd
from skbio import read, write, Protein


logger = getLogger(__name__)


def _filter_proteins(in_fp, out_fp, ids):
    '''Filter away the seq with specified IDs.'''
    with open(out_fp, 'w') as out:
        for seq in read(in_fp, format='fasta', constructor=Protein):
            seq_id = seq.metadata['id']
            if seq_id not in ids:
                write(seq, format='fasta', into=out)


def _add_protein_annotation(imd_dict, seq_hits):
    ''''''
    for seq_id, imd in imd_dict.items():
        hits = seq_hits[seq_id]
        for intvl in imd._intervals:
            md = intvl.metadata
            if md['type'] == 'CDS':
                # this md is parsed from prodigal output, which
                # has ID like "1_1", "1_2" for genes
                idx = md['ID'].split('_')[1]
                if idx in hits:
                    md['db_xref'] = hits[idx]


def _parse_seq_id(hits):
    '''prodigal outputs faa file with seq id of
    'gi|556503834|ref|NC_000913.3|_3224'. need to get the input seq id
    and the index for the protein seq

    '''
    d = defaultdict(dict)
    for row in hits.itertuples():
        seq_id, i = row[1].rsplit('_', 1)
        d[seq_id][i] = row[2]
    return d


def parse_diamond_uniref(fn, pident=90):
    ''''''
    def parse(columns=['qseqid', 'qlen', 'sseqid', 'slen',
                       'pident', 'length', 'gaps', 'evalue', 'bitscore',
                       'qstart', 'qend', 'sstart', 'send']):
        logger.debug('Parsing diamond search against UniRef')
        df = pd.read_table(fn, names=columns)
        df_filtered = filter_ident_overlap(df, pident)
        return df_filtered.loc[:, ['qseqid', 'sseqid']]
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
