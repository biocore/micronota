# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger
from sqlite3 import connect, IntegrityError
from xml.etree import ElementTree as ET


logger = getLogger(__name__)


def add_metadata(xml_fh, db_fp):
    '''Add to the database the metadata of records from the UniProt xml file.

    Parameters
    ----------
    in_fps : list of str
        The gzipped files of either UniProtKB Swiss-Prot or TrEMBLE.
    db_fp : str
        The output database file. See ``Notes``.

    Returns
    -------
    int
        The number of records processed.

    '''
    logger.info('Adding UniProt metadata to db')

    # this is the namespace for uniprot xml files.
    ns_map = {'xmlns': 'http://uniprot.org/uniprot',
              'xsi': 'http://WWW.w3.org/2001/XMLSchema-instance'}
    entry_tag = '{{{ns}}}{tag}'.format(ns=ns_map['xmlns'], tag='entry')
    paths = {'ec': './/xmlns:ecNumber',  # E.C. number
             'go': './xmlns:dbReference[@type="GO"]',  # GO
             'tigrfam': './xmlns:dbReference[@type="TIGRFAMs"]'}
    inserts = {}
    with connect(db_fp) as conn:
        c = conn.cursor()
        # The INTEGER PRIMARY KEY column created is simply an
        # alias for ROWID or _ROWID_ or OID.
        # You can't ignore this column because ROWID can't server
        # as foreign key.
        c.execute('CREATE TABLE IF NOT EXISTS uniprot ('
                  ' id   INTEGER PRIMARY KEY,'
                  ' accn TEXT  UNIQUE,'
                  ' name TEXT  NOT NULL);')
        insert = ('INSERT OR IGNORE INTO uniprot (id, accn, name)'
                  ' VALUES (?,?,?);')

        for other_db in paths:
            ct, it, clt, ilt = _cross_ref_table(other_db)
            c.execute(ct)
            c.execute(clt)
            inserts[other_db] = [it, ilt]
        for n, entry in enumerate(_parse_xml(xml_fh, entry_tag), 1):
            try:
                # get the primary accession number
                accn = entry.find('./xmlns:accession', ns_map).text
                # get the protein product name
                name = entry.find('.//xmlns:fullName', ns_map).text
            except AttributeError as e:
                raise Exception(
                    'failed to get accession and name '
                    'for record %d' % n) from e
            # customize with more informative error msg
            try:
                c.execute(insert, (None, accn, name))
            except IntegrityError as e:
                raise Exception(
                    'failed to insert {}'.format((accn, name))) from e
            uniprot_id = c.execute(
                'SELECT id FROM uniprot WHERE accn = ?;',
                (accn, )).fetchone()[0]
            for other_db, path in paths.items():
                it, ilt = inserts[other_db]
                select = 'SELECT id FROM %s WHERE accn = ?;' % other_db
                for elem in entry.findall(path, ns_map):
                    if other_db == 'ec':
                        other_accn = elem.text
                    else:
                        other_accn = elem.attrib['id']
                    c.execute(it, (None, other_accn, None))
                    row_id = c.execute(select, (other_accn, )).fetchone()[0]
                    c.execute(ilt, (uniprot_id, row_id))

        conn.commit()
    return n


def _cross_ref_table(name):
    '''sqlite3 statement to create table for cross ref database.

    name: other database name
    '''
    create = ('CREATE TABLE IF NOT EXISTS {} ('
              ' id   INTEGER PRIMARY KEY,'
              ' accn TEXT  UNIQUE,'
              ' name TEXT);').format(name)
    insert = ('INSERT OR IGNORE INTO {} (id, accn, name)'
              ' VALUES (?,?,?);').format(name)
    # junction table from uniprot to other annotation databases
    create_link_table = ('CREATE TABLE IF NOT EXISTS uniprot_{0} ('
                         ' uniprot_id INTEGER,'
                         ' {0}_id     INTEGER,'
                         ' PRIMARY KEY (uniprot_id, {0}_id),'
                         ' FOREIGN KEY (uniprot_id) REFERENCES uniprot(id),'
                         ' FOREIGN KEY ({0}_id) REFERENCES {0}(id));').format(
                             name)
    insert_link_table = ('INSERT OR IGNORE INTO uniprot_{0} (uniprot_id, {0}_id)'
                         ' VALUES (?,?);').format(name)
    return create, insert, create_link_table, insert_link_table


def _parse_xml(xml_fh, tag):
    '''Return the elem with specified tag.

    Parameters
    ----------
    xml_fh : xml file object or file path
    ns_map : dict
        namespace map
    '''
    # it is very important to set the events to 'end'; otherwise,
    # elem would be an incomplete record.
    for event, elem in ET.iterparse(xml_fh, events=['end']):
        if elem.tag == tag:
            yield elem
            # this is necessary for garbage collection
            elem.clear()


def _process_entry(elem, ns_map, path):
    '''Return the text content of an element matching the path.

    Parameters
    ----------
    elem : xml.etree.ElementTree.Element
    ns_map : dict
        namespace map
    '''
    match = elem.find(path, ns_map)
    if match is None:
        raise ValueError('Cannot find the path {} in elem: \n {}'.format(
            path,
            '\n'.join('%s:%s' % (child.tag, child.text) for child in elem)))

    return match.text


def query_xref(db_fp, accn):
    '''Query with UniProt accession number.'''
    tables = ['ec', 'go', 'tigrfam']
    with connect(db_fp) as conn:
        c = conn.cursor()
        c.execute('SELECT name FROM uniprot WHERE accn = ?;', (accn,))
        name = c.fetchone()[0]
        for table in tables:
            query_xref = '''SELECT {0}.accn FROM {0}
                        INNER JOIN uniprot_{0} j ON j.{0}_id = {0}.id
                        INNER JOIN uniprot u ON j.uniprot_id = u.id
                        WHERE u.accn = ?;'''.format(table)
            # 'K9NBS6'
            xref = []
            for i in c.execute(query_xref, (accn,)):
                xref.append(i[0])
            yield table, xref

