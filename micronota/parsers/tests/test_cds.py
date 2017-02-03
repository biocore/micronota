# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main, mock
from os.path import join
from tempfile import mkdtemp

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata
from skbio import Protein, write, read

from micronota.parsers.cds import _fetch_cds_metadata, parse_diamond_uniref


class Tests(TestCase):
    def setUp(self):
        self.tmpd = mkdtemp()

    @mock.patch('micronota.parsers.cds.format_xref', return_value={'db_xref': ['GO:0006351']})
    def test_fetch_cds_metadata(self, mock_format_xref):
        f = io.StringIO('query\tuniprot\nNC_016822.1_3793\tA3DA47\n')
        obs = list(_fetch_cds_metadata(f, 'foo'))
        exp = [('NC_016822.1', '3793', {'db_xref': ['GO:0006351', 'uniprot:A3DA47']})]
        self.assertEqual(obs, exp)

    @mock.patch('micronota.parsers.cds.format_xref', return_value={'EC_number': ['2.5.1.15']})
    def test_fetch_cds_metadata_empty_db_xref(self, mock_format_xref):
        f = io.StringIO('query\tuniprot\nNC_016822.1_3793\tA3DA47\n')
        obs = list(_fetch_cds_metadata(f, 'foo'))
        exp = [('NC_016822.1', '3793', {'db_xref': ['uniprot:A3DA47'], 'EC_number': ['2.5.1.15']})]
        self.assertEqual(obs, exp)

    def test_parse_diamond_uniref(self):
        lines = ('seq_3\t330\tA3DA47\t329\t87.8\t329\t1\t2.0e-167\t570.5\t1\t329\t1\t328\n'
                 'seq_6\t100\tA3DA48\t129\t47.8\t100\t1\t2.0e-16\t57.5\t1\t99\t1\t128\n'
                 # this entry has overlap < 80%
                 'seq_7\t200\tA3DA49\t229\t97.8\t200\t100\t2.0e-67\t70.5\t1\t189\t1\t228\n'
                 'seq_9\t200\tA3DA49\t229\t97.8\t200\t10\t2.0e-67\t70.5\t1\t189\t1\t228\n')
        f = io.StringIO(lines)
        df = parse_diamond_uniref(f, 50)()
        with io.StringIO() as fh:
            df.to_csv(fh, sep='\t', index=False)
            obs = fh.getvalue()
            exp = ('query\tuniprot\n'
                   'seq_3\tA3DA47\n'
                   'seq_9\tA3DA49\n')
            self.assertEqual(obs, exp)

        f = io.StringIO(lines)
        df = parse_diamond_uniref(f, 97.8)()
        with io.StringIO() as fh:
            df.to_csv(fh, sep='\t', index=False)
            obs = fh.getvalue()
            exp = ('query\tuniprot\n'
                   'seq_9\tA3DA49\n')
            self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()

