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
    make_db, search_protein_homologs)
from micronota.bfillings.util import _get_data_dir


class DiamondTests(TestCase):
    def setUp(self):
        self.db_fa = _get_data_dir()('db.faa')
        self.db = _get_data_dir()('db.dmnd')
        self.temp_dir = mkdtemp()
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


if __name__ == '__main__':
    main()
