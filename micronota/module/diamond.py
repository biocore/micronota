# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict
from sqlite3 import connect

import pandas as pd

from . import BaseMod
from ..database._util import query, format_xref


class Module(BaseMod):
    def __init__(self, directory, file_patterns=None):
        if file_patterns is None:
            file_patterns = {'hit': 'diamond.hit'}
        super().__init__(directory, file_patterns)

    def parse(self,
              metadata,
              columns=['qseqid', 'qlen', 'sseqid', 'slen',
                       'pident', 'length', 'gaps', 'evalue', 'bitscore',
                       'qstart', 'qend', 'sstart', 'send'],
              db='UniRef'):
        df = pd.read_table(self.files['hit'], names=columns)
        if db == 'UniRef':
            db = 'uniprot'
            # strip out the prefix 'UniRefXX_' from the IDs.
            df.sseqid = df.sseqid.apply(lambda s: s.split('_', 1)[-1])
        self.result = self._fetch_cds_metadata(df, db, metadata)

    @staticmethod
    def _fetch_cds_metadata(hits, db, metadata):
        '''Get metadata for the protein sequences matching any reference.

        Prodigal outputs faa file with seq id like
        'gi|556503834|ref|NC_000913.3|_3224'. It is needed to split to get
        the input seq id and the index for the protein seq

        '''
        protein = defaultdict(defaultdict)
        if metadata is None:
            for row in hits.itertuples():
                seq_id, i = row.qseqid.rsplit('_', 1)
                accn = row.sseqid
                hit = '{0}:{1}'.format(db, accn)
                protein[seq_id][i] = {'db_xref': hit}
            return protein
        else:
            with connect(metadata) as c:
                for row in hits.itertuples():
                    seq_id, i = row.qseqid.rsplit('_', 1)
                    accn = row.sseqid
                    hit = '{0}:{1}'.format(db, accn)
                    md = format_xref(query(c, db, accn))
                    if 'db_xref' in md:
                        md['db_xref'].append(hit)
                    else:
                        md['db_xref'] = [hit]
                    protein[seq_id][i] = md

            return protein
