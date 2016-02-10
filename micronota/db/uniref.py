r'''
UniRef
======

.. currentmodule:: micronota.db.uniref

This module create databases from UniRef for micronota usage. Using UniRef
for functional assignment of coding genes, we can get the clustering of
the annotated proteins for free. Using UniRef100 as the reference database,
we can easily collapse the clusters down to the similarity levels
of 90% or 50%. For those UniRef records from UniProKB, we also transfer
the metadata associated with those UniProKB records.

UniProtKB (UniProt Knowledgebase) [#]_ contains two parts, UniProtKB/Swiss-Prot
and UniProtKB/TrEMBL. The first part contains protein records that are
manually annotated and reviewed while the second part is done computationally
and not manually reviewed.

UniParc (UniProt Archive) [#]_ is also part of UniProt project. It is a
is a comprehensive and non-redundant database that contains most of the
publicly available protein sequences in the world.

The UniProt Reference Clusters (UniRef) [#]_ provide clustered sets (UniRef100,
UniRef90 and UniRef50 clusters) of sequences from the UniProtKB
and selected UniParc records, in order to obtain complete coverage of sequence
space at several resolutions (100%, >90% and >50%) while hiding redundant
sequences (but not their descriptions) from view.

Release
-------
:version: 2016_01
:date:    01/20/2016
:link:    http://www.uniprot.org/downloads

Files
-----
* UniRef files

  UniRef fasta files clustered at the similarity levels of 50, 90 or 100.

  * uniref50.fasta.gz [#]_

  * uniref90.fasta.gz [#]_

  * uniref100.fasta.gz [#]_

* UniProtKB files

  The .dat.gz files contain the UniProtKB records in variant format of EMBL.
  The records are divided into taxa groups [#]_.

  * uniprot_sprot_archaea.dat.gz

  * uniprot_sprot_bacteria.dat.gz

  * uniprot_trembl_archaea.dat.gz

  * uniprot_trembl_bacteria.dat.gz


Reference
---------
.. [#] http://www.uniprot.org/help/uniparc
.. [#] http://www.uniprot.org/help/uniprotkb
.. [#] http://www.ncbi.nlm.nih.gov/pubmed/17379688
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/README
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/README
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/README
.. [#] ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/README
'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, expanduser
from sqlite3 import connect

from skbio import read, Sequence


def prepare_db(out_d, prefix='uniref_sprot', force=False,
               server='ftp.uniprot.org',
               path='pub/databases/uniprot/uniref/uniref90',
               fasta='uniref90.fasta.gz'):
    '''
    '''
    prepare_metadata(join(expanduser('~'), 'Documents', 'uniprot_sprot.dat.gz'),
                     '%s.db' % join(out_d, prefix))


def prepare_metadata(in_fp, fp):
    '''
    Returns
    -------
    int
        The number of records processed.
    '''
    with connect(fp) as conn:
        conn.execute("DROP TABLE IF EXISTS uniref")
        conn.execute('''CREATE TABLE uniref (
                            id       TEXT    NOT NULL,
                            tag      TEXT    NOT NULL,
                            value    BLOB    NOT NULL,
                            transfer BOOLEAN NOT NULL,
                        CHECK (transfer IN (0,1)));''')

        for n, seq in enumerate(
                read(in_fp, format='embl', constructor=Sequence), 1):
            md = seq.metadata
            id = md['AC']
            insert = '''INSERT INTO uniref (id, tag, value, transfer)
                        VALUES (?,?,?,?);'''
            conn.execute(insert, (id, 'OX', md['OX'], 0))
            conn.execute(insert, (id, 'PE', md['PE'], 0))
            kingdom = md['OC'].split('; ', 1)[0]
            conn.execute(insert, (id, 'KD', kingdom, 0))
        # don't forget to index the column to speed up query
        conn.execute('CREATE INDEX id ON uniref (id);')
        conn.commit()

    return n
