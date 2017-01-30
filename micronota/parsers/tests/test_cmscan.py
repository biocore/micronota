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

from micronota.parsers.cmscan import parse


class ParseTests(TestCase):
    def test_parse(self):
        imd1 = IntervalMetadata(None)
        imd1.add(bounds=[(3588441, 3588818)],
                 metadata={'ncRNA_class': 'RNaseP_bact_a', 'type': 'ncRNA', 'strand': '-', 'db_xref': 'RF00010'})
        imd1.add(bounds=[(3355449, 3355633)],
                 metadata={'ncRNA_class': '6S', 'type': 'ncRNA', 'strand': '+', 'db_xref': 'RF00013'})
        imd2 = IntervalMetadata(None)
        imd2.add(bounds=[(85215, 85384)],
                 metadata={'ncRNA_class': 'U2', 'type': 'ncRNA', 'strand': '+', 'db_xref': 'RF00004'})
        imd3 = IntervalMetadata(None)
        imd3.add(bounds=[(8739, 8777)],
                 metadata={'ncRNA_class': 'U2', 'type': 'ncRNA', 'strand': '+', 'db_xref': 'RF00004'})
        exp = (('NC_016822.1', imd1),
               ('NC_016833.1', imd2),
               ('NC_016834.1', imd3))

        fp = get_data_path('cmscan.txt')
        gen = parse(fp)

        for (exp_id, exp_imd), (obs_id, obs_imd) in zip(exp, gen):
            self.assertEqual(exp_id, obs_id)
            self.assertEqual(exp_imd, obs_imd)


if __name__ == '__main__':
    main()

