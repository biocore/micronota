# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from logging import getLogger
from sqlite3 import connect, IntegrityError
from xml.etree import ElementTree as ET
from os.path import join, basename


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
              'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
    entry_tag = '{{{ns}}}{tag}'.format(ns=ns_map['xmlns'], tag='entry')

    with connect(db_fp) as conn:
        c = conn.cursor()
        # The INTEGER PRIMARY KEY column created is simply an
        # alias for ROWID or _ROWID_ or OID.
        # You can't ignore this column because ROWID can't server
        # as foreign key.
        c.execute('''CREATE TABLE IF NOT EXISTS uniprot (
                        id           INTEGER PRIMARY KEY,
                        primary_accn TEXT  UNIQUE,
                        name         TEXT  NOT NULL);''')
        uniprot_insert = '''INSERT INTO uniprot (id, primary_accn, name)
                            VALUES (?,?,?);'''

        c.execute('''CREATE TABLE IF NOT EXISTS ec (
                        id         INTEGER PRIMARY KEY,
                        ec_number  TEXT  UNIQUE,
                        name       TEXT  );''')
        ec_insert = '''INSERT OR IGNORE INTO ec (id, ec_number, name)
                       VALUES (?,?,?);'''

        c.execute('''CREATE TABLE IF NOT EXISTS uniprot_ec (
            id      INTEGER PRIMARY KEY,
            uniprot INTEGER,
            ec      INTEGER,
            FOREIGN KEY (uniprot) REFERENCES uniprot(id),
            FOREIGN KEY (ec) REFERENCES ec(id) );''')
        uniprot_ec_insert = '''INSERT INTO uniprot_ec (id, uniprot, ec)
                               VALUES (?,?,?);'''


        for n, entry in enumerate(_parse_xml(xml_fh, entry_tag), 1):
            # get the primary accession number
            accn = _process_entry(entry, ns_map, './xmlns:accession')
            # get the protein name
            name = _process_entry(entry, ns_map, './/xmlns:fullName')
            # customize with more informative error msg
            try:
                c.execute(uniprot_insert, (None, accn, name))
            except IntegrityError as e:
                raise Exception(
                    'failed trying to insert {}'.format((accn, name))) from e
            uniprot_id = c.lastrowid
            # find all E.C. numbers
            for ec_number in entry.findall('.//xmlns:ecNumber', ns_map):
                c.execute(ec_insert, (None, ec_number.text, None))
                ec_id = c.lastrowid
                c.execute(uniprot_ec_insert, (None, uniprot_id, ec_id))

        conn.commit()
    return n


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
