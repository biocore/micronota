# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from filecmp import cmp
from os.path import join
from shutil import rmtree, copy
from collections import namedtuple
from subprocess import CalledProcessError
from unittest import TestCase, main

from skbio.util import get_data_path

from micronota.bfillings.infernal import run_cmscan, run_cmpress
from micronota.util import _get_named_data_path


class InfernalTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()
        self.negative_fps = [get_data_path(i)
                             for i in ['empty', 'whitespace_only']]
        self.cm_fp = _get_named_data_path('test.cm')

    def tearDown(self):
        rmtree(self.tmp_dir)


class CMScanTests(InfernalTests):
    def setUp(self):
        super().setUp()
        suffices = ['cmscan']
        Files = namedtuple('Files', suffices)
        Case = namedtuple('Case', ['query', 'obs', 'exp'])
        self.cases = []
        for f in ['NC_018498.fna']:
            fns = ['{}.{}'.format(f, s) for s in suffices]
            obs_files = Files(*[join(self.tmp_dir, fn) for fn in fns])
            exp_files = Files(*[_get_named_data_path(fn) for fn in fns])
            self.cases.append(
                Case(_get_named_data_path(f), obs_files, exp_files))

    def test_run_cmscan_wrong_input(self):
        for fp in self.negative_fps:
            with self.assertRaises(CalledProcessError):
                run_cmscan(self.cm_fp, fp, join(self.tmp_dir, 'foo'))

    def test_run_cmscan(self):
        for case in self.cases:
            run_cmscan(self.cm_fp, case.query, case.obs.cmscan)
            with open(case.obs.cmscan) as obs, open(case.exp.cmscan) as exp:
                # skip comment lines as some contain running time info
                self.assertListEqual(
                    [i for i in exp.readlines() if not i.startswith('#')],
                    [j for j in obs.readlines() if not j.startswith('#')])


class CMPressTests(InfernalTests):
    def test_run_cmpress(self):
        # .i1i, ilm and i1p files are different from run to run. skip them.
        suffix = 'i1f'
        fp = copy(self.cm_fp, self.tmp_dir)
        run_cmpress(fp)
        fps = ['{}.{}'.format(i, suffix) for i in [fp, self.cm_fp]]
        self.assertTrue(cmp(*fps, shallow=False))
        # test not overwriting
        with self.assertRaises(CalledProcessError):
            run_cmpress(fp)


if __name__ == '__main__':
    main()
