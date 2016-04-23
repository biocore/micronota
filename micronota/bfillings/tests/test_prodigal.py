# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from shutil import rmtree
from os import listdir
from os.path import join
from unittest import TestCase, main
from filecmp import cmp
from collections import namedtuple
from subprocess import CalledProcessError

from skbio.util import get_data_path
from skbio.metadata import Feature

from micronota.util import _get_named_data_path
from micronota.bfillings.prodigal import run, _parse_faa


class RunTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()
        self.negative_fps = [get_data_path(i) for i in
                             ['empty', 'whitespace_only']]
        Case = namedtuple('Case', ['input', 'kwargs', 'outdir'])
        self.cases = [Case(i, j, k) for i, j, k in
                      zip([_get_named_data_path(i) for i in
                           # modified from NC_018498.gbk
                           ['NC_018498_partial_1.gbk',
                            'NC_018498_partial_1.gbk',
                            'NC_018498_partial_2.gbk']],
                          [{'-p': 'meta', '-f': 'gbk'},
                           {'-p': 'meta'},
                           {'-p': 'single'}],
                          ['test_1',
                           'test_2',
                           'test_3'])]

    def test_wrong_input_fp(self):
        for fp in self.negative_fps:
            with self.assertRaisesRegex(
                    CalledProcessError,
                    'returned non-zero exit status'):
                run(self.tmp_dir, input=fp)

    def test_run(self):
        for case in self.cases:
            exp_d = _get_named_data_path(case.outdir)
            obs_d = join(self.tmp_dir, case.outdir)
            run(obs_d, input=case.input, **case.kwargs)
            for f in listdir(exp_d):
                self.assertTrue(
                    cmp(join(obs_d, f), join(exp_d, f), shallow=False))

    def tearDown(self):
        # remove the tempdir and contents
        rmtree(self.tmp_dir)


class ParseTests(TestCase):
    def setUp(self):
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

    def test_pred_parse_faa(self):
        obs = _parse_faa(self.parse_fp)
        for e, o in zip(self.parse_exp, obs):
            self.assertEqual(e, o)


if __name__ == '__main__':
    main()
