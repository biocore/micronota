# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from os import getcwd
from shutil import rmtree
from os.path import join
from unittest import TestCase, main
from functools import partial
from skbio.util import get_data_path
from burrito.util import ApplicationError

from micronota.bfillings.minced import MinCED, predict_crispr


class MinCEDTests(TestCase):
    def setUp(self):
        self.temp_dir = mkdtemp()
        self.get_minced_path = partial(
            get_data_path, subfolder=join('data', 'minced'))

        # taken from MinCED test files
        self.positive_fps = list(map(self.get_minced_path,
                                 ['Aquifex_aeolicus_VF5.fna',
                                  'Aquifex_aeolicus_VF5.fna',
                                  'Aquifex_aeolicus_VF5.fna']))
        # 'empty' file raises JAVA (minced) error
        self.negative_fps = list(map(get_data_path, [
            'whitespace_only',
            'empty']))
        self.positive_params = [
            {'-searchWL': '8'},
            {'-searchWL': '8', '-minNR': '3'},
            {}]
        self.positive_flags = [
            {'gff': True, 'gffFull': False, 'spacers': False},
            {'gffFull': True, 'gff': False, 'spacers': False},
            {'gff': True, 'spacers': True, 'gffFull': False}]
        self.positive_prefix = 'Aquifex_aeolicus_VF5'

    def test_base_command(self):
        c = MinCED()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s' % (getcwd(), c._command))

    def test_predict_crispr_wrong_input(self):
        for fp in self.negative_fps:
            with self.assertRaisesRegex(
                    ApplicationError,
                    r'Error constructing CommandLineAppResult.'):
                predict_crispr(fp, self.temp_dir, 'foo')

    def test_predict_crispr(self):
        for fp, params, flags in zip(self.positive_fps, self.positive_params,
                                     self.positive_flags):
            prefix = self.positive_prefix
            res = predict_crispr(fp, self.temp_dir,
                                 prefix, params,
                                 gff=flags['gff'],
                                 gffFull=flags['gffFull'],
                                 spac=flags['spacers'])
            self.assertEqual(res['ExitStatus'], 0)
            if flags['gff']:
                suffix = 'gff'
            elif flags['gffFull']:
                suffix = 'gffFull'
            else:
                suffix = 'crisprs'
            fp = self.get_minced_path('.'.join([prefix, suffix]))
            with open(fp) as f:
                self.assertEqual(
                     # skip comment lines as some contain runtime info
                     [i for i in f.readlines()
                      if not i.startswith('Time')],
                     [j for j in res['output'].readlines()
                      if not j.startswith('Time')])
            # SPACERS flag produces an *additional* OUT_spacers.fa file
            # other flags produce OUT.FLAG outputs
            if flags['spacers']:
                suffix = 'spacers.fa'
                fp = self.get_minced_path('_'.join([prefix, suffix]))
                with open(fp) as f:
                    self.assertEqual(f.read(), res['spacers'].read())
            res['StdOut'].close()
            res['StdErr'].close()

    def tearDown(self):
        # remove the tempdir and contents
        rmtree(self.temp_dir)

if __name__ == '__main__':
    main()
