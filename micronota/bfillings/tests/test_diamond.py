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

import numpy as np
import skbio
from skbio.util import get_data_path

from micronota.util import _get_named_data_path
from micronota.bfillings.diamond import (
    run_view, run_makedb, run_blast, run,
    DiamondCache)


class DiamondTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        self.faa = _get_named_data_path('db.faa')
        self.obs_db = join(self.tmp_dir, 'db.dmnd')
        self.exp_db = _get_named_data_path('db.dmnd')

        cases = [('blastp', 'test.faa', 'test.faa.daa', 'test.faa.tab', 'test.faa.sam'),
                 ('blastx', 'test.fna', 'test.fna.daa', 'test.fna.tab', 'test.fna.sam')]
        Case = namedtuple('Case', ['aligner', 'query', 'daa', 'tab', 'sam'])
        self.cases = [Case(*i) for i in cases]

        self.neg_fp = [get_data_path(i) for i in
                       ['empty', 'whitespace_only']]

    # def tearDown(self):
    #     rmtree(self.tmp_dir)


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
            obs_daa = join(self.tmp_dir, case.daa)
            # exp_daa = _get_named_data_path(case.daa)
            run_blast(aligner=case.aligner,
                      query=case.query, daa=obs_daa, db=self.exp_db,
                      tmpdir=self.tmp_dir)
            # self.assertTrue(cmp(obs_daa, exp_daa))

            obs_tab = join(self.tmp_dir, case.tab)
            exp_tab = _get_named_data_path(case.tab)
            run_view(daa=obs_daa, out=obs_tab, fmt='tab')
            self.assertTrue(cmp(obs_tab, exp_tab))

            obs_sam = join(self.tmp_dir, case.sam)
            exp_sam = _get_named_data_path(case.sam)
            run_view(daa=obs_daa, out=obs_sam, fmt='sam')
            self.assertTrue(cmp(obs_sam, exp_sam))

    def test_run(self):
        for case in self.cases:
            run(self.tmp_dir, case.query, self.exp_db, aligner=case.aligner, tmpdir=self.tmp_dir)
            obs = join(self.tmp_dir, 'db.sam')
            self.assertTrue(cmp(obs, _get_named_data_path(case.sam)))


class DiamondCacheTests(DiamondTests):
    def setUp(self):
        super().setUp()
        cases = [('blastp', 'WP_009885814.faa'),
                 ('blastx', 'WP_009885814.fna')]
        Test = namedtuple('Test', ['aligner', 'input', 'exp'])
        self.tests = [Test(i[0],
                           get_data_path(i[1]),
                           _get_named_data_path('%s.diamond' % i[1]))
                      for i in cases]

        seqs = skbio.read(_get_named_data_path('cache.faa'), format='fasta')
        self.cache = DiamondCache(list(seqs))

    def test_cache(self):
        np.random.seed(0)
        for test in self.tests:
            aligner, query, exp_fp = test.aligner, test.input, test.exp

            pred = FeatureAnnt([self.db], mkdtemp(dir=self.tmp_dir),
                               cache=self.cache)
            obs = pred(query, aligner=aligner)
            exp = pred._filter_best(pred.parse_tabular(exp_fp))
            self.assertSetEqual(set(exp['sseqid'].values),
                                set(obs['sseqid'].values))

    def test_cache_initialize(self):
        np.random.seed(0)
        for test in self.tests:
            aligner, query = test.aligner, test.input
            pred = FeatureAnnt([self.db], mkdtemp(dir=self.tmp_dir),
                               cache=DiamondCache())
            pred(query, aligner=aligner)
            self.assertTrue(len(pred.cache.seqs) > 0)

    def test_cache_empty_db(self):
        np.random.seed(0)
        for test in self.tests:
            aligner, query, exp_fp = test.aligner, test.input, test.exp
            pred = FeatureAnnt([], mkdtemp(dir=self.tmp_dir),
                               cache=self.cache)
            obs = pred(query, aligner=aligner)
            exp = pred._filter_best(pred.parse_tabular(exp_fp))
            self.assertEqual(exp['sseqid'].values, obs['sseqid'].values)


class ParsingTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()
        cases = ['WP_009885814.fna', 'WP_009885814.faa']
        Test = namedtuple('Test', ['input', 'exp', 'obs'])
        self.sam_tests = [Test(_get_named_data_path('%s.sam' % i),
                               _get_named_data_path('%s.txt' % i),
                               join(self.tmp_dir, '%s.txt' % i))
                          for i in cases]
        self.filter_tests = [Test(_get_named_data_path('%s.diamond' % i),
                                  _get_named_data_path('%s.best' % i),
                                  join(self.tmp_dir, '%s.best'))
                             for i in cases]
        self.filter_tests2 = [Test(_get_named_data_path('%s.sam' % i),
                                   _get_named_data_path('%s.idcov' % i),
                                   join(self.tmp_dir, '%s.idcov'))
                              for i in cases]

    def test_parse_sam(self):
        for test in self.sam_tests:
            df = FeatureAnnt.parse_sam(test.input)
            df.to_csv(test.obs, sep='\t', index=False)
            self.assertTrue(cmp(test.exp, test.obs, shallow=False))

    def test_filter_best(self):
        for test in self.filter_tests:
            df = FeatureAnnt.parse_tabular(test.input)
            df_filter = FeatureAnnt._filter_best(df)
            df_filter.to_csv(test.obs, sep='\t')
            self.assertTrue(cmp(test.exp, test.obs, shallow=False))

    def test_filter_id_cov(self):
        for test in self.filter_tests2:
            df = FeatureAnnt.parse_sam(test.input)
            df_filter = FeatureAnnt._filter_id_cov(df, pident=30, cov=92)
            df_filter.to_csv(test.obs, sep='\t')
            self.assertTrue(cmp(test.exp, test.obs, shallow=False))


if __name__ == '__main__':
    main()
