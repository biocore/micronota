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
from skbio.metadata import Feature
from burrito.util import ApplicationError

from micronota.bfillings.prodigal import (
    Prodigal, identify_features, _parse_faa)


class ProdigalTests(TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.get_prodigal_path = partial(
            get_data_path, subfolder=join('data', 'prodigal'))

        self.positive_fps = list(map(self.get_prodigal_path, [
            # modified from NC_018498.gbk
            'NC_018498_partial_1.gbk',
            'NC_018498_partial_1.gbk',
            'NC_018498_partial_2.gbk',
            ]))
        self.positive_params = [
            {'-p': 'meta'},
            {'-p': 'meta', '-f': 'gff'},
            {'-p': 'single'}]
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

        self.parse_fp = self.get_prodigal_path('output.faa')
        self.parse_exp = [
            {Feature(type_='CDS',
                     id='1_1',
                     right_partial_=False,
                     left_partial_=False,
                     location='686..1828',
                     translation='MKILINKSELNKILKKMNNVIISNNKIKPHHSYFLIEAKEKEINFYANNEYFSVKCNLNKYFLITSKSEPELKQILVPSR*',
                     note='"start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.236"',
                     rc_=False): [(685, 1828)],
             Feature(type_='CDS',
                     id='1_2',
                     location='1828..>2757',
                     translation='MNLYDLLELPTTASIKEIKIAYKRLAKRYHPDVNKLGSQTFVEINNAYSILSDPNQKEKYFNYKTQHFID',
                     right_partial_=True,
                     left_partial_=False,
                     note='"start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.271"',
                     rc_=False): [(1827, 2757)]},

            {Feature(type_='CDS',
                     id='2_1',
                     location='21577..22128',
                     right_partial_=False,
                     left_partial_=False,
                     translation='MKKTSPFILRRTKNKVLKELPKKIITDIYVELSEEHQKLYDKQKTDGLKEIKESDAKNALFDV*',
                     note='"start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.272"',
                     rc_=False): [(21576, 22128)]}]

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

    def test_identify_features_wrong_input(self):
        for fp in self.negative_fps:
            with self.assertRaisesRegex(
                    ApplicationError,
                    r'(Sequence read failed)|(no input sequences to analyze)'):
                identify_features(fp, self.temp_dir, 'foo')

    def test_identify_features(self):
        for fp, params, prefix, suffix in zip(self.positive_fps,
                                              self.positive_params,
                                              self.positive_prefices,
                                              self.positive_suffices):
            res = identify_features(fp, self.temp_dir, prefix, params)
            self.assertEqual(res['ExitStatus'], 0)
            for i in ['-o', '-d', '-a']:
                fp = self.get_prodigal_path('.'.join([prefix, suffix[i]]))
                with open(fp) as f:
                    self.assertEqual(f.read(), res[i].read())
                res[i].close()
            res['StdOut'].close()
            res['StdErr'].close()

    def test_parse_faa(self):
        obs = _parse_faa(self.parse_fp)
        for e, o in zip(self.parse_exp, obs):
            self.assertEqual(e, o)

    def tearDown(self):
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)


if __name__ == '__main__':
    main()
