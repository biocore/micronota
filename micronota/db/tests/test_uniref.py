# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, dirname, exists
from os import remove
from tempfile import mkdtemp
from unittest import main
from shutil import rmtree

from micronota.util import _DBTest, _get_named_data_path
from micronota.db._uniref import create_metadata, sort_uniref
from micronota.db.uniref100 import prepare_db


class UnirefTests(_DBTest):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        self.db_fp = 'uniprotkb.db'
        self.obs_db_fp = join(self.tmp_dir, self.db_fp)
        self.exp_db_fp = _get_named_data_path(self.db_fp)
        self.uniprotkb = [_get_named_data_path('uniprot_sprot.xml.gz'),
                          _get_named_data_path('uniprot_trembl.xml.gz'),
                          12]
        self.d = dirname(self.exp_db_fp)

        self.uniref_fp = _get_named_data_path('uniref100.fasta.gz')
        self.uniref_res = [
            'Swiss-Prot_Archaea',
            'Swiss-Prot_Bacteria',
            'Swiss-Prot_Viruses',
            'Swiss-Prot_Eukaryota',
            'Swiss-Prot_other',
            'TrEMBL_Archaea',
            'TrEMBL_Bacteria',
            'TrEMBL_Viruses',
            'TrEMBL_Eukaryota',
            'TrEMBL_other',
            '_other']

    def test_prepare_metadata(self):
        n = create_metadata(self.uniprotkb[:2], self.obs_db_fp)
        self.assertEqual(n, self.uniprotkb[2])
        self._test_eq_db(self.obs_db_fp, self.exp_db_fp)

    def _test_eq(self):
        for fp in self.uniref_res:
            for suffix in ['fasta', 'dmnd']:
                fp = '.'.join([fp, suffix])
                obs = join(self.tmp_dir, fp)
                exp = _get_named_data_path(fp)
                if exists(exp):
                    with open(obs) as o, open(exp) as e:
                        self.assertEqual(o.read(), e.read())
                    remove(obs)

    def test_sort_uniref(self):
        sort_uniref(self.exp_db_fp, self.uniref_fp,
                    join(self.tmp_dir, 'uniref100'), 100)
        self._test_eq()

    def test_prepare_db(self):
        prepare_db(self.d, self.tmp_dir)
        self._test_eq()
        self._test_eq_db(self.obs_db_fp, self.exp_db_fp)

    def test_prepare_db_not_overwrite(self):
        with self.assertRaisesRegex(
                FileExistsError, r'The file .* exists.'):
            prepare_db(self.d, self.d)

    def tearDown(self):
        rmtree(self.tmp_dir)


if __name__ == '__main__':
    main()
