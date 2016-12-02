from sqlite3 import connect


def query(sql_fp, db, accn):
    '''Query with accession number of a reference db.

    Parameters
    ----------
    sql_fp : str
        sqlite file path
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
    with connect(sql_fp) as c:
        tables = {i[0] for i in c.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        tables.discard(db)
        db1 = db + '_'
        db2 = '_' + db
        for table in tables:
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
            if res and 'product' not in info:
                name = [i[1] for i in res]
                info['product'] = name[0]
            info[other] = [i[0] for i in res]
    return info


def _format_xref(d):
    old_ec, new_ec = 'ec_number', 'EC_number'
    if old_ec in d:
        d[new_ec] = d.pop(old_ec)
    replaces = [('tigrfam', 'TIGRFAM'),
                ('eggnog', 'eggNOG'),
                ('kegg', 'KO'),
                ('pfam', 'Pfam')]
    # GO id already in the form of "GO:1290834". only need to add to xref
    xref = d.pop('go', [])
    # For other ref dbs, ref db name need to prefix their IDs
    for old, new in replaces:
        xref.extend(['{0}:{1}'.format(new, i) for i in d.pop(old, [])])
    d['db_xref'] = xref
    return d
