# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from filecmp import cmp
from shutil import rmtree
from os.path import join, splitext
from unittest import TestCase, main
from subprocess import CalledProcessError

from skbio.util import get_data_path

from micronota.bfillings.minced import run
from micronota.util import _get_named_data_path


class MinCEDTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        self.negative_fps = [get_data_path(i)
                             for i in ['empty', 'whitespace_only']]

    def test_run_wrong_input(self):
        for fp in self.negative_fps:
            # 'empty' file raises JAVA (minced) error
            with self.assertRaises(CalledProcessError):
                run(fp, join(self.tmp_dir, 'foo'))

    def test_run(self):
        print(self.tmp_dir)
        # taken from MinCED test files
        fn = 'Aquifex_aeolicus_VF5.fna'
        query = _get_named_data_path(fn)
        exp = splitext(query)[0]
        obs = join(self.tmp_dir, splitext(fn)[0])
        params = [
            {'-searchWL': '8', 'gff': True, 'gffFull': False, 'spacers': False},
            {'-searchWL': '8', '-minNR': '3', 'gffFull': True, 'gff': False, 'spacers': False},
            {'gff': True, 'spacers': True, 'gffFull': False}]

        exp_fps = ['.gff', '.gffFull', '.gff']
        for param, s in zip(params, exp_fps):
            run(query, obs, **param)
            self.assertTrue(
                cmp(obs, exp + s))

    def tearDown(self):
        rmtree(self.tmp_dir)


if __name__ == '__main__':
    main()
