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

from micronota.bfillings.hmmer import run, run_hmmpress
from micronota.util import _get_named_data_path


class HmmerTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()
        self.hmm_fp = _get_named_data_path('Pfam_B_1.hmm')
        self.negative_fps = [get_data_path(i)
                             for i in ['empty', 'whitespace_only']]

    def tearDown(self):
        rmtree(self.tmp_dir)


class HMMScanTests(HmmerTests):
    def setUp(self):
        super().setUp()
        suffices = ['tblout']
        Files = namedtuple('Files', suffices)
        Case = namedtuple('Case', ['query', 'obs', 'exp'])
        self.cases = []
        for f in ['Pfam_B_1.faa']:
            fns = ['{}.{}'.format(f, s) for s in suffices]
            obs_files = Files(*[join(self.tmp_dir, fn) for fn in fns])
            exp_files = Files(*[_get_named_data_path(fn) for fn in fns])
            self.cases.append(
                Case(_get_named_data_path(f), obs_files, exp_files))

    def test_run_wrong_input(self):
        for fp in self.negative_fps:
            with self.assertRaises(CalledProcessError):
                run(self.hmm_fp, fp, self.tmp_dir)

    def test_run(self):
        for case in self.cases:
            res = run(self.hmm_fp, case.query, self.tmp_dir)
            with open(res.params['out'].value) as obs, open(case.exp.tblout) as exp:
                # skip comment lines as some contain running time info
                self.assertListEqual(
                    [i for i in exp.readlines() if not i.startswith('#')],
                    [j for j in obs.readlines() if not j.startswith('#')])


class HmmpressTests(HmmerTests):
    def test_run_hmmpress(self):
        # .i1i, ilm and i1p files are different from run to run. skip them.
        suffix = 'h3f'
        fp = copy(self.hmm_fp, self.tmp_dir)
        run_hmmpress(fp)
        fps = ['{}.{}'.format(i, suffix) for i in [fp, self.hmm_fp]]
        self.assertTrue(cmp(*fps, shallow=False))
        # test not overwriting
        with self.assertRaises(CalledProcessError):
            run_hmmpress(fp)


if __name__ == "__main__":
    main()
