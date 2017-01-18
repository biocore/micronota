# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import main, TestCase
from os.path import dirname

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata

from micronota.parsers.aragorn import parse


class AragornTests(TestCase):
    def test_parse(self):
        imd1 = IntervalMetadata(None)
        imd1.add(bounds=[(237929, 238006)],
                 metadata={'strand': '+', 'type': 'tRNA', 'source': 'Aragorn', 'product': 'tRNA-Ile'})
        imd1.add(bounds=[(238048, 238124)],
                 metadata={'strand': '+', 'type': 'tRNA', 'source': 'Aragorn', 'product': 'tRNA-Ala'})
        imd2 = IntervalMetadata(None)
        imd2.add(bounds=[(4954141, 4954228)],
                 metadata={'strand': '-', 'type': 'tRNA', 'source': 'Aragorn', 'product': 'tRNA-Ser'})
        imd3 = IntervalMetadata(None)
        exp = (('NC_016822.1', imd1), ('NC_016833.1', imd2), ('NC_016834.1', imd3))

        fn = 'aragorn.txt'
        fp = get_data_path(fn)
        gen = parse(dirname(fp), fn)

        for (exp_id, exp_imd), (obs_id, obs_imd) in zip(exp, gen):
            self.assertEqual(exp_id, obs_id)
            self.assertEqual(exp_imd, obs_imd)


if __name__ == '__main__':
    main()

