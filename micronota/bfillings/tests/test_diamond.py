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

from micronota.bfillings.diamond import (
    DiamondMakeDB,
    make_db,
    FeatureAnnt)
from micronota.bfillings.util import _get_data_dir
import pandas as pd
import pandas.util.testing as pdt


class DiamondTests(TestCase):
    def setUp(self):
        self.db_fa = _get_data_dir()('db.faa')
        self.db = _get_data_dir()('db.dmnd')
        self.temp_dir = mkdtemp()
        self.out_dir = mkdtemp()
        self.neg_fp = map(
            get_data_path,
            ['empty', 'whitespace_only'])

    def tearDown(self):
        rmtree(self.temp_dir)


class DiamondMakeDBTests(DiamondTests):
    def test_base_command(self):
        c = DiamondMakeDB()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s' % (getcwd(), c._command))

    def test_make_db(self):
        fp = join(self.temp_dir, 'db.dmnd')
        res = make_db(self.db_fa, fp)
        with open(fp, 'rb') as obs, open(self.db, 'rb') as exp:
            self.assertEqual(obs.read(), exp.read())
        self.assertEqual(res, 0)

    def test_make_db_wrong_input(self):
        fp = join(self.temp_dir, 'db.dmnd')
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
            (i[0],
             get_data_path(i[1]),
             _get_data_dir()('%s.diamond' % i[1]))
            for i in tests]

    def test_blast(self):
        for aligner, query, exp in self.blast:
            res = search_protein_homologs(
                query, self.db, self.temp_dir, aligner)
            with open(res) as o, open(exp) as e:
                self.assertEqual(o.read(), e.read())

    def test_blast_wrong_input(self):
        for i in self.neg_fp:
            for aligner in ['blastp', 'blastx']:
                with self.assertRaisesRegex(
                        ApplicationError,
                        r'(Error reading file)|(Invalid input file format)'):
                    search_protein_homologs(
                        i, self.db, self.temp_dir, aligner)

class FeatAnntTests(TestCase):
    def setUp(self):
        #fat = FeatureAnnt(dat=[self.db],out_dir=self.out_dir)
        tests = [('blastp', 'WP_009885814.faa'),
                 ('blastx', 'WP_009885814.fna')]
        self.blast = [
            (i[0],
             get_data_path(i[1]),
             _get_data_dir()('%s.sam' % i[1]))
            for i in tests]

        self.exp = \
            pd.DataFrame({
                'sseqid' : ['UniRef100_P47599', 'UniRef100_B2HPZ3',
                            'UniRef100_A4T166'],
                'evalue' : [2.1e-229, 2.9e-58, 3.3e-57],
                'bitscore' : [2009, 533, 524],
                'sequence' :[
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
