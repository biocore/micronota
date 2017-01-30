# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

import pandas as pd
from skbio import read, write, Sequence

from ..database._util import query, format_xref


logger = getLogger(__name__)


def _add_cds_metadata(seqs, cds_metadata):
    '''Add metadata to all the CDS interval features.'''
    for seq_id, seq in seqs.items():
        hits = cds_metadata[seq_id]
        for intvl in seq.interval_metadata._intervals:
            md = intvl.metadata
            md_type = md.get('type', '')
            if md_type == 'CDS':
                # this md is parsed from prodigal output, which
                # has ID like "seq1_1", "seq1_2" for genes
                idx = md['ID'].split('_')[1]
                if idx in hits:
                    md.update(hits[idx])


def _fetch_cds_metadata(hit_f, db):
    '''Get metadata for the protein sequences matching any reference.

    Prodigal outputs faa file with seq id like
    'gi|556503834|ref|NC_000913.3|_3224'. It is needed to split to get
    the input seq id and the index for the protein seq

    '''
    hits = pd.read_table(hit_f)
    ref = hits.columns[1]
    for row in hits.itertuples():
        seq_id, i = row[1].rsplit('_', 1)
        accn = row[2]
        hit = '{0}:{1}'.format(ref, accn)
        md = format_xref(query(db, ref, accn))
        if 'db_xref' in md:
            md['db_xref'].append(hit)
        else:
            md['db_xref'] = [hit]
        yield seq_id, i, md


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
