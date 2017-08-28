# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from tempfile import NamedTemporaryFile
from unittest import main

from skbio.util import get_data_path

from micronota.util import _DBTest
from micronota.database.uniprot import add_metadata, _parse_xml


class UniProtTests(_DBTest):
    def setUp(self):
        self.uniprot = [get_data_path(i) for i in ['uniprot_sprot', 'uniprot_trembl']]
        self.uniprot_xml = ['%s.xml.gz' % i for i in self.uniprot]
        self.uniprot_db = ['%s.sqlite' % i for i in self.uniprot]
        self.uniprot_n = [6, 6]
        self.uniprot_all = get_data_path('uniprot.sqlite')

    def test_parse_xml(self):
        tag = '{http://uniprot.org/uniprot}entry'
        for xml, n in zip(self.uniprot_xml, self.uniprot_n):
            with gzip.open(xml) as fh:
                entries = _parse_xml(fh, tag=tag)
                for i, s in enumerate(entries, 1):
                    self.assertEqual(s.tag, tag)
                # there are 6 entries in the each xml file
                self.assertEqual(i, n)

    def test_add_metadata(self):
        for xml, n, db in zip(self.uniprot_xml, self.uniprot_n, self.uniprot_db):
            with gzip.open(xml) as fh:
                with NamedTemporaryFile() as tmp:
                    count = add_metadata(fh, tmp.name)
                    self._test_eq_db(tmp.name, db)
                    self.assertEqual(n, count)

    def test_add_metadata_both(self):
        n = 0
        with NamedTemporaryFile() as tmp:
            for xml in self.uniprot_xml:
                with gzip.open(xml) as fh:
                    n += add_metadata(fh, tmp.name)
            self._test_eq_db(tmp.name, self.uniprot_all)
            self.assertEqual(n, sum(self.uniprot_n))


if __name__ == '__main__':
    main()
