# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join
from tempfile import mkdtemp
from shutil import rmtree

from skbio import write, read, Sequence
from skbio.metadata import IntervalMetadata
from skbio.util import get_data_path

from micronota.summarize import summarize, summarize_iter


class Tests(TestCase):
    def setUp(self):
        self.tmpd = mkdtemp()
        self.input_genbank_fp = get_data_path('GCF_000010365.1.gbff')

    def test_summarize(self):
        gb = Sequence.read(self.input_genbank_fp, format='genbank')
        # test length
        obs = summarize(gb, types=('length',))
        exp = [159662]
        self.assertListEqual(obs, exp)
        # test nucleotide frequencies
        obs = summarize(gb, types=('nuc_freq',))
        exp = [{'A': 66734, 'C': 13501, 'G': 12946, 'T': 66481}]
        self.assertListEqual(obs, exp)
        # test number of CDS'
        obs = summarize(gb, types=('CDS',))
        exp = [175]
        self.assertListEqual(obs, exp)
        # test numbers of ncRNA's, rRNA's, and tRNAs
        obs = summarize(gb, types=('ncRNA', 'rRNA', 'tRNA'))
        exp = [0, 2, 24]
        self.assertListEqual(obs, exp)
        # test invalid object
        self.assertRaises(AttributeError, summarize, 'not_a_sequence',
                          types=('CDS',))

    def test_summarize_iter(self):
        gb = Sequence.read(self.input_genbank_fp, format='genbank')
        # test length
        obs = list(summarize_iter([gb] * 3, types=('length',)))
        exp = [[159662]] * 3
        self.assertListEqual(obs, exp)

    def tearDown(self):
        rmtree(self.tmpd)


if __name__ == '__main__':
    main()
