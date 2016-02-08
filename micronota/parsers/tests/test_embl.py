# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from pprint import pprint
from skbio.util import get_data_path

from micronota.parsers.embl import _embl_to_protein, _embl_to_generator


class EmblIOTests(TestCase):
    def setUp(self):
        self.multi_fp = get_data_path('uniprot_multi.embl')


class ReaderTests(EmblIOTests):
    def test_embl_to_protein(self):
        a = get_data_path('uniprot_single.embl')
        embl = _embl_to_protein(a)
        pprint(embl)

    def test_embl_to_generator(self):
        for obs in _embl_to_generator(self.multi_fp):
            pprint(obs)


if __name__ == '__main__':
    main()
