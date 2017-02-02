# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import main, TestCase
from sqlite3 import connect

from skbio.util import get_data_path

from micronota.database._util import query, format_xref


class Tests(TestCase):
    def setUp(self):
        self.db = get_data_path('uniprot.sqlite')

    def test_query_empty(self):
        with connect(self.db) as c:
            obs = query(c, 'uniprot', 'foo')
            exp = {}
            self.assertEqual(exp, obs)

    def test_query(self):
        with connect(self.db) as c:
            obs = query(c, 'uniprot', 'F7PH07')
            exp = {'eggNOG': ['arCOG02817', 'arCOG10292', 'COG0294'],
                   'KEGG': ['hti:HTIA_0709'],
                   'product': 'Dihydropteroate synthase protein',
                   'EC_number': ['2.5.1.15'],
                   'GO': ['GO:0004156', 'GO:0016787'],
                   'Pfam': ['PF06283']}
        self.assertEqual(exp, obs)

    def test_format_xref_empty(self):
        d = {'GO': [], 'TIGRFAM': [], 'eggNOG': [], 'KEGG': [], 'Pfam': []}
        obs = format_xref(d)
        exp = {}
        self.assertEqual(obs, exp)

    def test_format_xref(self):
        d = {'GO': ['GO:0016021', 'GO:0005886', 'GO:0030246', 'GO:0009405'],
             'KEGG': ['bmb:BruAb2_0497'], 'Pfam':['PF07886'],
             'EC_number': ['2.5.1.15'],
             'product': 'Lectin-like protein BA14k'}
        obs = format_xref(d)
        exp = {'db_xref': ['GO:0016021', 'GO:0005886', 'GO:0030246', 'GO:0009405',
                           'KEGG:bmb:BruAb2_0497', 'Pfam:PF07886'],
               'EC_number': ['2.5.1.15'],
               'product': 'Lectin-like protein BA14k'}

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
