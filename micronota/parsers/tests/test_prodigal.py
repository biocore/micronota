# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import dirname

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata

from micronota.parsers.prodigal import parse


class ParseTests(TestCase):
    def test_parse(self):
        imd1 = IntervalMetadata(None)
        imd1.add(bounds=[(336, 2799)],
                 metadata={'phase': 0, 'strand': '+', 'partial': '00',
                           'rbs_motif': 'GGAG/GAGG', 'source': 'Prodigal_v2.60',
                           'rscore': '10.55', 'conf': '99.99', 'tscore': '3.83',
                           'gc_cont': '0.533', 'score': '326.87', 'ID': '1_1',
                           'sscore': '14.74', 'uscore': '0.36',
                           'start_type': 'ATG', 'type': 'CDS',
                           'cscore': '312.13', 'rbs_spacer': '5-10bp'})

        exp = (('NC_016822.1', imd1),)

        fn = 'prodigal.gff'
        fp = get_data_path(fn)
        gen = parse(dirname(fp), fn)

        for (exp_id, exp_imd), (obs_id, obs_imd) in zip(exp, gen):
            self.assertEqual(exp_id, obs_id)
            self.assertEqual(exp_imd, obs_imd)


if __name__ == '__main__':
    main()

