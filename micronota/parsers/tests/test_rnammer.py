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

from micronota.parsers.rnammer import parse


class ParseTests(TestCase):
    def test_parse(self):
        imd1 = IntervalMetadata(None)
        imd1.add(bounds=[(0, 2853)],
                 metadata={'source': 'RNAmmer-1.2', 'type': '23s_rRNA', 'strand': '-', 'score': '3222.8'})
        imd1.add(bounds=[(2924, 3040)],
                 metadata={'source': 'RNAmmer-1.2', 'type': '5s_rRNA', 'strand': '+', 'score': '80.8'})
        imd2 = IntervalMetadata(None)
        imd2.add(bounds=[(77272, 78834)],
                 metadata={'source': 'RNAmmer-1.2', 'type': '16s_rRNA', 'strand': '+', 'score': '1984.2'})

        exp = (('NZ_JXDA01000005.1', imd1),
               ('NZ_JXDA01000001.1', imd2))

        fp = get_data_path('rnammer.gff')
        gen = parse(fp)

        for (exp_id, exp_imd), (obs_id, obs_imd) in zip(exp, gen):
            self.assertEqual(exp_id, obs_id)
            self.assertEqual(exp_imd, obs_imd)


if __name__ == '__main__':
    main()
