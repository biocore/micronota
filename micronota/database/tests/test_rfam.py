# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkdtemp
from os.path import join, splitext, exists
import io

from skbio.util import get_data_path

from micronota.database.rfam import filter_models


class Tests(TestCase):
    def setUp(self):
        self.i = get_data_path('rfam.cm')

    def test_filter_models(self):
        with io.StringIO() as ofile, open(self.i) as ifile:
            filter_models(ifile, ofile)
            obs = ofile.getvalue()
        # starting from this line, no more records are filtered
        n = 2023
        with open(self.i) as ifile:
            exp = ''.join(ifile.readlines()[n:])
        self.assertEqual(obs, exp)

    def test_filter_models_negate(self):
        with io.StringIO() as ofile, open(self.i) as ifile:
            filter_models(ifile, ofile, True)
            obs = ofile.getvalue()
        # starting from this line, no more records are filtered
        n = 2023
        with open(self.i) as ifile:
            exp = ''.join(ifile.readlines()[:n])
        self.assertEqual(obs, exp)

    def test_filter_models_zero(self):
        with io.StringIO() as ofile, open(self.i) as ifile:
            filter_models(ifile, ofile, models=set())
            obs = ofile.getvalue()
        with open(self.i) as ifile:
            exp = ifile.read()
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
