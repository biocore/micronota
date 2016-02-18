# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, dirname
from os import remove
from tempfile import mkdtemp
from unittest import main
from shutil import rmtree

from micronota.bfillings.util import _get_data_dir
from micronota.util import _DBTest
from micronota.db.uniref import (
    prepare_db, prepare_metadata, sort_uniref)


class UnirefTests(_DBTest):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        self.db_fp = 'uniprotkb.db'
        self.obs_db_fp = join(self.tmp_dir, self.db_fp)
        self.exp_db_fp = _get_data_dir()(self.db_fp)
        self.uniprotkb = [_get_data_dir()('uniprot_sprot.xml.gz'),
                          _get_data_dir()('uniprot_trembl.xml.gz'),
                          12]
        self.d = dirname(self.exp_db_fp)

        self.uniref_fp = _get_data_dir()('uniref100.fasta.gz')
        self.uniref_res = [
            'uniref100_Swiss-Prot_Archaea.fasta',
            'uniref100_Swiss-Prot_Bacteria.fasta',
            'uniref100_Swiss-Prot_Viruses.fasta',
            'uniref100_Swiss-Prot_Eukaryota.fasta',
            'uniref100_Swiss-Prot_other.fasta',
            'uniref100_TrEMBL_Archaea.fasta',
            'uniref100_TrEMBL_Bacteria.fasta',
            'uniref100_TrEMBL_Viruses.fasta',
            'uniref100_TrEMBL_Eukaryota.fasta',
            'uniref100_TrEMBL_other.fasta',
            'uniref100__other.fasta']

    def test_prepare_metadata(self):
        n = prepare_metadata(self.uniprotkb[:2], self.obs_db_fp)
        self.assertEqual(n, self.uniprotkb[2])
        self._test_eq_db(self.obs_db_fp, self.exp_db_fp)

    def _test_eq(self):
        for fp in self.uniref_res:
            obs = join(self.tmp_dir, fp)
            exp = _get_data_dir()(fp)
            with open(obs) as o, open(exp) as e:
                self.assertEqual(o.read(), e.read())
            remove(obs)

    def test_sort_uniref(self):
        sort_uniref(self.exp_db_fp, self.uniref_fp, self.tmp_dir)
        self._test_eq()

    def test_prepare_db(self):
        prepare_db(self.tmp_dir, self.d)
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
