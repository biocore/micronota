# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from shutil import rmtree
from os import getcwd, listdir
from os.path import join
from unittest import TestCase, main
from filecmp import cmp

from skbio.util import get_data_path
from skbio.metadata import Feature
from burrito.util import ApplicationError

from micronota.util import _get_named_data_path
from micronota.bfillings.prodigal import (
    Prodigal, FeaturePred)


class ProdigalTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        self.positive_fps = [_get_named_data_path(i) for i in [
            # modified from NC_018498.gbk
            'NC_018498_partial_1.gbk',
            'NC_018498_partial_1.gbk',
            'NC_018498_partial_2.gbk',
            ]]
        self.positive_params = [
            {'-p': 'meta'},
            {'-p': 'meta', '-f': 'gff'},
            {'-p': 'single'}]
        self.positive_suffices = [
            {'-o': 'gbk', '-a': 'faa', '-d': 'fna'},
            {'-o': 'gff', '-a': 'faa', '-d': 'fna'},
            {'-o': 'gbk', '-a': 'faa', '-d': 'fna'}]
        self.positive_outdir = [
            'test_1',
            'test_2',
            'test_3']

        self.negative_fps = [get_data_path(i) for i in [
            'empty',
            'whitespace_only']]

        self.parse_fp = _get_named_data_path('parse_test.faa')
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

    def test_pred_wrong_input_fp(self):
        pred = FeaturePred(None, self.tmp_dir)
        for fp in self.negative_fps:
            with self.assertRaisesRegex(
                    ApplicationError,
                    r'(Sequence read failed)|(no input sequences to analyze)'):
                pred(fp)

    def test_pred_run(self):
        for fp, params, outdir in zip(self.positive_fps,
                                      self.positive_params,
                                      self.positive_outdir):
            exp_d = _get_named_data_path(outdir)
            obs_d = join(self.tmp_dir, outdir)
            pred = FeaturePred(None, obs_d)
            res = pred.run(fp, params)
            self.assertEqual(res['ExitStatus'], 0)
            for f in listdir(exp_d):
                self.assertTrue(cmp(join(obs_d, f), join(exp_d, f)))
            res['StdOut'].close()
            res['StdErr'].close()

    def test_pred_parse_faa(self):
        pred = FeaturePred(None, self.tmp_dir)
        obs = pred._parse_faa(self.parse_fp)
        for e, o in zip(self.parse_exp, obs):
            self.assertEqual(e, o)

    def tearDown(self):
        # remove the tempdir and contents
        rmtree(self.tmp_dir)


if __name__ == '__main__':
    main()
