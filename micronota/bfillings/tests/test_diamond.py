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
from os.path import join
from collections import namedtuple
from unittest import TestCase, main
from subprocess import CalledProcessError

from skbio.util import get_data_path

from micronota.util import _get_named_data_path
from micronota.bfillings.diamond import (
    run_view, run_makedb, run_blast, run,
    parse_sam, filter_best, filter_ident_overlap)


class DiamondTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        self.faa = _get_named_data_path('db.faa')
        self.obs_db = join(self.tmp_dir, 'db.dmnd')
        self.exp_db = _get_named_data_path('db.dmnd')

        cases = [['blastp', 'test.faa']]
        # ['blastx', 'test.fna']]
        suffices = ['daa', 'tab', 'sam', 'sam_df', 'sam_best', 'sam_io5', 'sam_io12']
        Files = namedtuple('Files', suffices)
        Case = namedtuple('Case', ['aligner', 'query', 'obs', 'exp'])
        self.cases = []
        for b, f in cases:
            fns = ['{}.{}'.format(f, s) for s in suffices]
            obs_files = Files(*[join(self.tmp_dir, fn) for fn in fns])
            exp_files = Files(*[_get_named_data_path(fn) for fn in fns])
            self.cases.append(Case(b, _get_named_data_path(f), obs_files, exp_files))

        self.neg_fp = [get_data_path(i) for i in
                       ['empty', 'whitespace_only']]

    def tearDown(self):
        rmtree(self.tmp_dir)


class RunTests(DiamondTests):
    def test_run_makedb(self):
        run_makedb(self.faa, self.obs_db)
        with open(self.obs_db, 'rb') as obs, open(self.exp_db, 'rb') as exp:
            self.assertEqual(obs.read(), exp.read())

    def test_run_makedb_wrong_input(self):
        for i in self.neg_fp:
            with self.assertRaisesRegex(
                    CalledProcessError,
                    r'returned non-zero exit status 1'):
                run_makedb(i, self.obs_db)

    def test_blast_view(self):
        for case in self.cases:
            run_blast(aligner=case.aligner,
                      query=case.query, daa=case.obs.daa, db=self.exp_db,
                      tmpdir=self.tmp_dir)

            run_view(daa=case.obs.daa, out=case.obs.tab, fmt='tab')
            self.assertTrue(
                cmp(case.obs.tab, case.exp.tab, shallow=False))

            run_view(daa=case.obs.daa, out=case.obs.sam, fmt='sam')
            self.assertTrue(
                cmp(case.obs.sam, case.exp.sam, shallow=False))

    def test_run(self):
        for case in self.cases:
            run(self.tmp_dir, case.query, self.exp_db, aligner=case.aligner, tmpdir=self.tmp_dir)
            obs = join(self.tmp_dir, 'db.sam')
            self.assertTrue(
                cmp(obs, _get_named_data_path(case.exp.sam), shallow=False))


class ParsingTests(DiamondTests):
    def test_parse_sam(self):
        for case in self.cases:
            df = parse_sam(case.exp.sam)
            df.to_csv(case.exp.sam_df, sep='\t', index=False)
            df.to_csv(case.obs.sam_df, sep='\t', index=False)
            self.assertTrue(
                cmp(case.exp.sam_df, case.obs.sam_df, shallow=False))

    def test_filter_best(self):
        for case in self.cases:
            df = parse_sam(case.exp.sam)
            df_filter = filter_best(df)
            df_filter.to_csv(case.exp.sam_best, sep='\t', index=False)
            df_filter.to_csv(case.obs.sam_best, sep='\t', index=False)
            self.assertTrue(
                cmp(case.exp.sam_best, case.obs.sam_best, shallow=False))

    def test_filter_ident_overlap(self):
        for case in self.cases:
            df = parse_sam(case.exp.sam)
            df_filter = filter_ident_overlap(df, pident=90, overlap=5)
            df_filter.to_csv(case.exp.sam_io5, sep='\t', index=False)
            df_filter.to_csv(case.obs.sam_io5, sep='\t', index=False)
            self.assertTrue(
                cmp(case.exp.sam_io5, case.obs.sam_io5, shallow=False))

            df_filter = filter_ident_overlap(df, pident=90, overlap=12)
            df_filter.to_csv(case.exp.sam_io12, sep='\t', index=False)
            df_filter.to_csv(case.obs.sam_io12, sep='\t', index=False)
            self.assertTrue(
                cmp(case.exp.sam_io12, case.obs.sam_io12, shallow=False))


if __name__ == '__main__':
    main()
