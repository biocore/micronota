# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


def query(c, db, accn):
    '''Query with accession number of a reference db.

    Parameters
    ----------
    c : ``sqlite3.Connection``
        connection to the db that has the entry metadata
    db : str
        reference db name (eg uniprot, tigrfam, kegg, etc)
    accn : str
        the accession number in the ref db

    Returns
    -------
    dict
        The cross-ref IDs to other db and product of the accession.
    '''
    info = {}
    tables = {i[0] for i in c.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    tables.discard(db)
    db1 = db + '_'
    db2 = '_' + db
    for table in tables:
        # find all the junction tables that is linked to the ref table
        if db1 in table:
            other = table.replace(db1, '')
        elif db2 in table:
            other = table.replace(db2, '')
        else:
            continue
        query_xref = '''SELECT {1}.accn, t.name FROM {1}
                        INNER JOIN {2} j ON j.{1}_id = {1}.id
                        INNER JOIN {0} t ON j.{0}_id = t.id
                        WHERE t.accn = ?;'''.format(db, other, table)
        # 'K9NBS6'
        res = [i for i in c.execute(query_xref, (accn,))]
        if res:
            if 'product' not in info:
                name = [i[1] for i in res]
                info['product'] = name[0]
            info[other] = [i[0] for i in res]
    return info


def format_xref(d):
    '''format the metadata

    Parameters
    ----------
    d : dict
        metadata for a gene obtained from ``query`` function

    Returns
    -------
    dict
        the formatted metadata convenient for output
    '''
    refs = ('TIGRFAM', 'eggNOG', 'KEGG', 'Pfam')
    # GO id already in the form of "GO:1290834". only need to add to xref
    xref = d.pop('GO', [])
    # For other ref dbs, ref db name need to prefix their IDs
    for ref in refs:
        xref.extend(['{0}:{1}'.format(ref, i) for i in d.pop(ref, [])])
    if xref:
        d['db_xref'] = xref
    return d
