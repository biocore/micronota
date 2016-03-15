# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from shutil import rmtree
from os import getcwd
from os.path import join
from unittest import TestCase, main

from skbio.util import get_data_path
from burrito.util import ApplicationError

from micronota.util import _get_named_data_path
from micronota.bfillings.diamond import (
    DiamondMakeDB, make_db, FeatureAnnt,
    DiamondCache)
import pandas as pd
import pandas.util.testing as pdt
import numpy as np
import skbio


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
        tests = [('blastp', 'WP_009885814.faa'),
                 ('blastx', 'WP_009885814.fna')]
        self.blast = [
            (i[0], get_data_path(i[1]),
             _get_named_data_path('%s.diamond' % i[1]))
            for i in tests]

    def test_blast(self):
        for aligner, query, exp_fp in self.blast:
            pred = FeatureAnnt([self.db], mkdtemp(dir=self.tmp_dir))
            obs = pred(query, aligner=aligner)
            exp = pred.parse_tabular(exp_fp)
            self.assertTrue(exp.equals(obs))

    def test_blast_wrong_input(self):
        pred = FeatureAnnt([self.db], self.tmp_dir)
        for i in self.neg_fp:
            for aligner in ['blastp', 'blastx']:
                with self.assertRaisesRegex(
                        ApplicationError,
                        r'(Error reading file)|(Invalid input file format)'):
                    pred(i, aligner=aligner)


class TestDiamondCache(DiamondTests):
    def setUp(self):
        super().setUp()
        tests = ('blastp', 'WP_009885814.faa')
        self.blast = \
            (tests[0], get_data_path(tests[1]),
             _get_named_data_path('%s.diamond' % tests[1]))
        seqs = skbio.read(_get_named_data_path('cache.faa'), format='fasta')
        self.cache = DiamondCache(list(seqs))

    def test_cache(self):
        np.random.seed(0)
        aligner, query, exp_fp = self.blast
        pred = FeatureAnnt([self.db], mkdtemp(dir=self.tmp_dir),
                           cache=self.cache)
        obs = pred(query, aligner=aligner)
        exp = pred.parse_tabular(exp_fp)
        self.assertEquals(exp['sseqid'].values, obs['sseqid'].values)

    def test_cache_empty_db(self):
        np.random.seed(0)
        aligner, query, exp_fp = self.blast
        pred = FeatureAnnt([], mkdtemp(dir=self.tmp_dir),
                           cache=self.cache)
        obs = pred(query, aligner=aligner)
        exp = pred.parse_tabular(exp_fp)
        self.assertEquals(exp['sseqid'].values, obs['sseqid'].values)


class TestParseSam(TestCase):
    def setUp(self):
        tests = [('blastp', 'WP_009885814.faa'),
                 ('blastx', 'WP_009885814.fna')]
        self.blast = [
            (i[0],
             get_data_path(i[1]),
             _get_named_data_path(('%s.sam' % i[1])))
            for i in tests]

        self.exp = \
            pd.DataFrame({
                'sseqid': ['UniRef100_P47599', 'UniRef100_B2HPZ3',
                           'UniRef100_A4T166'],
                'evalue': [2.1e-229, 2.9e-58, 3.3e-57],
                'bitscore': [2009, 533, 524],
                'sequence': [
                            'MQSHKILVVNAGSSSIKFQLFNDKKQVLAKGLCERIFIDGFFKLEFNQK'
                            'KIEEKVQFNDHNLAVKHFLNALKKNKIITELSEIGLIGHRVVQGANYFT'
                            'DAVLVDTHSLAKIKEFIKLAPLHNKPEADVIEIFLKEIKTAKNVAVFDT'
                            'TFHTTIPRENYLYAVPENXEKNNLVRRYGFHGTSYKYINEFLEKKFNKK'
                            'PLNLIVCHLGNGASVCAIKQGKSLNTSMGFTPLEGLIMGTRSGDIDPAI'
                            'VSYIAEQQKLSCNDVVNELNKKSGMFAITGSSDMRDIFDKPEINDIAIK'
                            'MYVNRVADYIAKYLNQLSGEIDSLVFTGGVGENASYCVQLIIEKVASLG'
                            'FKTNSNLFGNYQDSSLISTNESKYQIFRVRTNEELMIVEDALRVSTNIK'
                            'K',
                            'ILVVNAGSSSIKFQLFNDKKQVLAKGLCERIFIDGFFKLEFNQKKIEEK'
                            'VQFNDHNLAVKHFLNALKKNKIITELSEIGLIGHRVVQGANYFTDAVLV'
                            'DTHSLAKIKEFIKLAPLHNKPEADVIEIFLKEIKTAKNVAVFDTTFHTT'
                            'IPRENYLYAVPENXEKNNLVRRYGFHGTSYKYINEFLEKKFNKKPLNLI'
                            'VCHLGNGASVCAIKQGKSLNTSMGFTPLEGLIMGTRSGDIDPAIVSYIA'
                            'EQQKLSCNDVVNELNKKSGMFAITGSSDMRDIFDKPEINDIAIKMYVNR'
                            'VADYIAKYLNQLSGEIDSLVFTGGVGENASYCVQLIIEKVASLGFKTNS'
                            'NLFGNYQDSSLISTNESKYQIFRVRTNEELMIVEDALRV',
                            'ILVVNAGSSSIKFQLFNDKKQVLAKGLCERIFIDGFFKLEFNQKKIEEK'
                            'VQFNDHNLAVKHFLNALKKNKIITELSEIGLIGHRVVQGANYFTDAVLV'
                            'DTHSLAKIKEFIKLAPLHNKPEADVIEIFLKEIKTAKNVAVFDTTFHTT'
                            'IPRENYLYAVPENXEKNNLVRRYGFHGTSYKYINEFLEKKFNKKPLNLI'
                            'VCHLGNGASVCAIKQGKSLNTSMGFTPLEGLIMGTRSGDIDPAIVSYIA'
                            'EQQKLSCNDVVNELNKKSGMFAITGSSDMRDIFDKPEINDIAIKMYVNR'
                            'VADYIAKYLNQLSGEIDSLVFTGGVGENASYCVQLIIEKVASLGFKTNS'
                            'NLFGNYQDSSLISTNESKYQIFRVRTNEELMI']
                })

    def test_parse_sam(self):
        for test in self.blast:
            df = FeatureAnnt.parse_sam(test[2])
            df = df.reindex_axis(sorted(df.columns), axis=1)
            exp = df.reindex_axis(sorted(self.exp.columns), axis=1)

            pdt.assert_frame_equal(df, exp)

    def test_parse_sam_best(self):
        for test in self.blast:
            df = FeatureAnnt.parse_sam(test[2], column='bitscore')
            df = df.reindex_axis(sorted(df.columns), axis=1)
            exp = df.reindex_axis(sorted(self.exp.columns), axis=1)

            pdt.assert_frame_equal(df, exp)


if __name__ == '__main__':
    main()
