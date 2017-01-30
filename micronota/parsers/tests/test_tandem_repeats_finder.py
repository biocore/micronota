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

from micronota.parsers.tandem_repeats_finder import parse


class ParseTests(TestCase):
    def test_parse(self):
        imd1 = IntervalMetadata(None)
        imd1.add(bounds=[(252151, 252184)],
                 metadata={'source': 'Tandem_Repeats_Finder', 'repeat': 'T', 'type': 'tandem_repeat'})
        imd1.add(bounds=[(261169, 261210)],
                 metadata={'source': 'Tandem_Repeats_Finder', 'repeat': 'CTCTGA', 'type': 'tandem_repeat'})
        imd2 = IntervalMetadata(None)
        imd2.add(bounds=[(172614, 172703)],
                 metadata={'source': 'Tandem_Repeats_Finder', 'repeat': 'AACAGCCGC', 'type': 'tandem_repeat'})
        exp = (('NC_016822.1', imd1), ('NC_016833.1', imd2))

        fp = get_data_path('tandem_repeats_finder.txt')
        gen = parse(fp)

        for (exp_id, exp_imd), (obs_id, obs_imd) in zip(exp, gen):
            self.assertEqual(exp_id, obs_id)
            self.assertEqual(exp_imd, obs_imd)


if __name__ == '__main__':
    main()

