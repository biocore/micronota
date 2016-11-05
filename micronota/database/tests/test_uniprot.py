# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import gzip
from os.path import join, dirname, exists
from os import remove
from tempfile import mkstemp
from unittest import main
from shutil import rmtree

from micronota.util import _DBTest, _get_named_data_path
from micronota.database.uniprot import add_metadata, _parse_xml


class UniprotTests(_DBTest):
    def setUp(self):
        _, self.tmpf = mkstemp()
        self.uniprot = [_get_named_data_path(i) for i in
                        ['uniprot_sprot', 'uniprot_trembl']]
        self.uniprot_xml = ['%s.xml.gz' % i for i in self.uniprot]
        self.uniprot_db = ['%s.sqlite' % i for i in self.uniprot]
        self.uniprot_n = [6, 6]
        self.uniprot_all = _get_named_data_path('uniprot.sqlite')

    def test_parse_xml(self):
        tag = '{http://uniprot.org/uniprot}entry'
        for xml, n in zip(self.uniprot_xml, self.uniprot_n):
            with gzip.open(xml) as fh:
                entries = _parse_xml(fh, tag=tag)
                for i, s in enumerate(entries, 1):
                    self.assertEqual(s.tag, tag)
                # there are 6 entries in the input file
                self.assertEqual(i, n)

    def test_add_metadata(self):
        with gzip.open(self.uniprot_xml[0]) as fh:
            n = add_metadata(fh, self.tmpf)
        self._test_eq_db(self.tmpf, self.uniprot_db[0])
        self.assertEqual(n, self.uniprot_n[0])

    def test_add_metadata_both(self):
        n = 0
        for xml, i in zip(self.uniprot_xml, self.uniprot_n):
            with gzip.open(xml) as fh:
                n += add_metadata(fh, self.tmpf)
        self._test_eq_db(self.tmpf, self.uniprot_all)
        self.assertEqual(n, sum(self.uniprot_n))


if __name__ == '__main__':
    main()
