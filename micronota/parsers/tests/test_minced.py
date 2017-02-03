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

from micronota.parsers.minced import parse


class ParseTests(TestCase):
    def test_parse(self):
        imd1 = IntervalMetadata(None)
        imd1.add(bounds=[(3175627, 3175900)],
                 metadata={'strand': '.', 'type': 'CRISPR', 'ID': 'CRISPR1', 'source': 'minced:0.2.0', 'score': '5'})
        exp = (('NC_016822.1', imd1),)

        fp = get_data_path('minced.gff')
        gen = parse(fp)

        for (exp_id, exp_imd), (obs_id, obs_imd) in zip(exp, gen):
            self.assertEqual(exp_id, obs_id)
            self.assertEqual(exp_imd, obs_imd)


if __name__ == '__main__':
    main()
