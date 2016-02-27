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
from os.path import join
from unittest import TestCase, main
from functools import partial
from skbio.util import get_data_path
from burrito.util import ApplicationError

from micronota.bfillings.minced import MinCED, predict_crispr


class MinCEDTests(TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.get_minced_path = partial(
            get_data_path, subfolder=join('data', 'minced'))

        self.positive_fps = list(map(self.get_minced_path, [
            # taken from MinCED test files
            'Aquifex_aeolicus_VF5.fna'
            ]))
        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only']))
        self.positive_params = [
            {'-searchWL': '8'},
            {'-searchWL': '9', '-minNR': '2'},
            {'-minRL': '20', '-maxRL': '45'},
            {}]
        self.positive_prefices = [
            '-gff',
            '-spacers']

    def test_base_command(self):
        c = MinCED()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s' % (getcwd(), c._command))

    def test_predict_crispr_wrong_input(self):
        for fp in self.negative_fps:
            with self.assertRaisesRegex(
                    ApplicationError,
                    r'Sequence read failed \(file must be Fasta format\).'):
                predict_crispr(fp, self.temp_dir, 'foo')

    def test_predict_crispr(self):
        for fp, params, prefix in zip(self.positive_fps,
                                      self.positive_params,
                                      self.positive_prefices):
            res = predict_crispr(fp, self.temp_dir, prefix, params)
            self.assertEqual(res['ExitStatus'], 0)
            for i in ['.expected', '_spacers.fa']:
                fp = self.get_minced_path(''.join([prefix, i]))
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
