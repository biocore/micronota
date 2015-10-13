#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import shutil
from os import getcwd
from unittest import TestCase, main

from skbio.util import get_data_path
from burrito.util import ApplicationError

from micronota.bfillings.prodigal import Prodigal, predict_genes


class ProdigalTests(TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.positive_fps = list(map(get_data_path, [
            # modified from NC_018498.gbk
            'NC_018498_partial_1.gbk',
            'NC_018498_partial_1.gbk',
            'NC_018498_partial_2.gbk',
            ]))
        self.positive_params = [
            {'-p': 'meta'},
            {'-p': 'meta', '-f': 'gff'},
            {}]
        self.positive_prefices = [
            'prodigal_meta',
            'prodigal_meta',
            'prodigal_single']
        self.positive_suffices = [
            {'-o': 'gbk', '-a': 'faa', '-d': 'fna'},
            {'-o': 'gff', '-a': 'faa', '-d': 'fna'},
            {'-o': 'gbk', '-a': 'faa', '-d': 'fna'}]

        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only']))

    def test_base_command(self):
        c = Prodigal()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s' % (getcwd(), c._command))

        c.Parameters['-q'].on()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s -q' % (getcwd(), c._command))

        c.Parameters['-q'].off()
        c.Parameters['-p'].on('single')
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s -p single' % (getcwd(), c._command))

    def test_predict_genes_wrong_input(self):
        for fp in self.negative_fps:
            with self.assertRaisesRegex(
                    ApplicationError,
                    r'Sequence read failed \(file must be Fasta, '
                    'Genbank, or EMBL format\).'):
                predict_genes(fp, self.temp_dir, 'foo')

    def test_predict_genes(self):
        for fp, params, prefix, suffix in zip(self.positive_fps,
                                              self.positive_params,
                                              self.positive_prefices,
                                              self.positive_suffices):
            res = predict_genes(fp, self.temp_dir, prefix, params)
            self.assertEqual(res['ExitStatus'], 0)
            for i in ['-o', '-d', '-a']:
                fp = get_data_path('.'.join([prefix, suffix[i]]))
                with open(fp) as f:
                    self.assertEqual(f.read(), res[i].read())
                res[i].close()
            res['StdOut'].close()
            res['StdErr'].close()

    def tearDown(self):
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)


if __name__ == '__main__':
    main()
