r'''
UniRef
======

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

from skbio import read


def prepare_db(out_d, prefix='uniref', force=False):
    '''
    '''
    prepare_metadata(join(expanduser('~'), 'uniref'), prefix)


def prepare_metadata(in_fp, fp):
    '''
    '''
    with connect(fp) as conn:
        conn.execute("DROP TABLE IF EXISTS uniref")
        conn.execute('''CREATE TABLE uniref (
                            id       TEXT    NOT NULL,
                            tag      TEXT    NOT NULL,
                            value    BLOB    NOT NULL,
                            transfer BOOLEAN NOT NULL,
                        CHECK (transfer IN (0,1)));''')

        for seq in read(in_fp, format='embl'):
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
