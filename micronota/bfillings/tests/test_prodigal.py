#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import getcwd
from os.path import join
import tempfile
import shutil
from unittest import TestCase, main
from micronota.bfillings.prodigal import Prodigal, predict_genes
from skbio.util import get_data_path
from skbio import read
from burrito.util import ApplicationError


class ProdigalTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'prodigal_genome',
            'prodigal_metagenome']))
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
            with self.assertRaisesRegexp(
                    ApplicationError,
                    r'Sequence read failed \(file must be Fasta, '
                    'Genbank, or EMBL format\).'):
                predict_genes({'-i': fp})

    def test_predict_genes_meta(self):
        pass
