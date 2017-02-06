# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata

from micronota.parsers.transtermhp import parse


class ParseTests(TestCase):
    def test_parse(self):
        fp = get_data_path('transtermhp.tt')
        exps = [('gi|556503834|ref|NC_000913.3|1', 12),
                ('gi|556503834|ref|NC_000913.3|2', 4)]
        for (sid, imd), exp in zip(parse(fp), exps):
            self.assertEqual(sid, exp[0])
            self.assertEqual(imd.num_interval_features, exp[1])

        # test the interval metadata from the 2nd sequence
        exp_imd = IntervalMetadata(None)
        exp_imd.add([(7857, 7876)],
                    metadata={'strand': '+', 'confidence': '95',
                              'sequence': 'TGCCCACGATTAAAG/GTGGCCGC/CCTG/GCGGTCAC/TTCTTTGAGAAAAGG',
                              'source': 'TransTermHP',
                              'ID': 'TERM_1', 'gene_id': '2_7',
                              'type': 'terminator'})
        exp_imd.add([(8919, 8958)],
                    metadata={'strand': '-', 'confidence': '100',
                              'sequence': 'AATGAGCCAGAATAA/GCTAAGGTTGAAGGGGC/TGGAAC/GCCCCTTCAACCTTAGC/AGTAGCGTGGGATGA',
                              'source': 'TransTermHP',
                              'ID': 'TERM_2', 'gene_id': '2_9',
                              'type': 'terminator'})
        exp_imd.add([(9258, 9273)],
                    metadata={'strand': '+', 'confidence': '89',
                              'sequence': 'GGCAGAAACAAAAAA/TCCCCG/GACT/CGGGGA/TTTATGTACAAGAGG',
                              'ID': 'TERM_3', 'gene_id': '2_9',
                              'source': 'TransTermHP',
                              'type': 'terminator'})
        exp_imd.add([(9258, 9273)],
                    metadata={'strand': '-', 'confidence': '100',
                              'sequence': 'GGCAGAAACAAAAAA/TCCCCG/GACT/CGGGGA/TTTATGTACAAGAGG',
                              'ID': 'TERM_4', 'gene_id': '2_9',
                              'source': 'TransTermHP',
                              'type': 'terminator'})
        self.assertEqual(exp_imd, imd)


if __name__ == '__main__':
    main()
