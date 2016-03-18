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
from os import getcwd
from os.path import join
from collections import namedtuple
from unittest import TestCase, main

import numpy as np
import skbio
from skbio.util import get_data_path
from burrito.util import ApplicationError

from micronota.util import _get_named_data_path
from micronota.bfillings.diamond import (
    DiamondMakeDB, make_db, FeatureAnnt,
    DiamondCache)


class DiamondTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()
        self.db_fa = _get_named_data_path('db.faa')
        self.db = _get_named_data_path('db.dmnd')
        self.neg_fp = [get_data_path(i) for i in
                       ['empty', 'whitespace_only']]

    def tearDown(self):
        rmtree(self.tmp_dir)


class DiamondMakeDBTests(DiamondTests):
    def test_base_command(self):
        c = DiamondMakeDB()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s' % (getcwd(), c._command))

    def test_make_db(self):
        fp = join(self.tmp_dir, 'db.dmnd')
        make_db(self.db_fa, fp)
        with open(fp, 'rb') as obs, open(self.db, 'rb') as exp:
            self.assertEqual(obs.read(), exp.read())

    def test_make_db_wrong_input(self):
        fp = join(self.tmp_dir, 'db.dmnd')
        for i in self.neg_fp:
            with self.assertRaisesRegex(
                    ApplicationError,
                    r'(Error reading file)|(Invalid input file format)'):
                make_db(i, fp)


class DiamondBlastTests(DiamondTests):
    def setUp(self):
        super().setUp()
        cases = [('blastp', 'WP_009885814.faa'),
                 ('blastx', 'WP_009885814.fna')]
        Test = namedtuple('Test', ['aligner', 'input', 'exp'])
        self.tests = [Test(i[0],
                           get_data_path(i[1]),
                           _get_named_data_path('%s.diamond' % i[1]))
                      for i in cases]

    def test_blast(self):
        for test in self.tests:
            pred = FeatureAnnt([self.db], mkdtemp(dir=self.tmp_dir))
            obs = pred(test.input, aligner=test.aligner)
            exp = pred._filter_best(pred.parse_tabular(test.exp))
            self.assertTrue(exp.equals(obs))

    def test_blast_wrong_input(self):
        pred = FeatureAnnt([self.db], self.tmp_dir)
        for i in self.neg_fp:
            for aligner in ['blastp', 'blastx']:
                with self.assertRaisesRegex(
                        ApplicationError,
                        r'(Error reading file)|(Invalid input file format)'):
                    pred(i, aligner=aligner)


class DiamondCacheTests(DiamondTests):
    def setUp(self):
        super().setUp()
        tests = ('blastp', 'WP_009885814.faa')
        self.blast = (tests[0], get_data_path(tests[1]),
                      _get_named_data_path('%s.diamond' % tests[1]))
        seqs = skbio.read(_get_named_data_path('cache.faa'), format='fasta')
        self.cache = DiamondCache(list(seqs))

    def test_cache(self):
        np.random.seed(0)
        aligner, query, exp_fp = self.blast
        pred = FeatureAnnt([self.db], mkdtemp(dir=self.tmp_dir),
                           cache=self.cache)
        obs = pred(query, aligner=aligner)
        exp = pred._filter_best(pred.parse_tabular(exp_fp))
        self.assertEqual(exp['sseqid'].values, obs['sseqid'].values)

    def test_cache_empty_db(self):
        np.random.seed(0)
        aligner, query, exp_fp = self.blast
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
            df_filter.to_csv(test.obs, sep='\t', index=False)
            self.assertTrue(cmp(test.exp, test.obs, shallow=False))

    def test_filter_id_cov(self):
        for test in self.filter_tests2:
            df = FeatureAnnt.parse_sam(test.input)
            df_filter = FeatureAnnt._filter_id_cov(df, pident=30, cov=92)
            df_filter.to_csv(test.obs, sep='\t', index=False)
            self.assertTrue(cmp(test.exp, test.obs, shallow=False))


if __name__ == '__main__':
    main()
