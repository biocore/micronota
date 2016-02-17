# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from os import remove
from tempfile import mkstemp, mkdtemp
from unittest import TestCase, main
from sqlite3 import connect
import shutil

from micronota.bfillings.util import _get_data_dir
from micronota.db.uniref import (
    prepare_metadata, create_id_map, sort_uniref)


class UnirefTests(TestCase):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        _, self.obs_db_fp = mkstemp()
        self.exp_db_fp = _get_data_dir()('uniprokb.db')
        self.sprot = [2, _get_data_dir()('uniprot_sprot.xml.gz')]
        self.trembl = [2, _get_data_dir()('uniprot_trembl.xml.gz')]
        self.table_name = 'metadata'

        _, self.id_map_obs = mkstemp()
        self.id_map_exp = _get_data_dir()('id_map.db')
        self.id_map = [4, _get_data_dir()('id_map.txt.gz')]
        self.id_map_table = 'id_map'

        self.uniref_fp = _get_data_dir()('uniref100.fasta.gz')
        self.uniref_res = [
            'uniref100_sprot_archaea.fasta',
            'uniref100_sprot_bacteria.fasta',
            'uniref100_sprot_viruses.fasta',
            'uniref100_sprot_other.fasta',
            'uniref100_trembl_archaea.fasta',
            'uniref100_trembl_bacteria.fasta',
            'uniref100_trembl_viruses.fasta',
            'uniref100_trembl_other.fasta',
            'uniref100__other.fasta']

    def test_prepare_metadata(self):
        n_sprot = prepare_metadata(self.sprot[1], self.obs_db_fp)
        self.assertEqual(n_sprot, self.sprot[0])
        n_trembl = prepare_metadata(self.trembl[1], self.obs_db_fp)
        self.assertEqual(n_trembl, self.trembl[0])

        with connect(self.obs_db_fp) as o, connect(self.exp_db_fp) as e:
            co = o.cursor()
            co.execute('SELECT * from %s' % self.table_name)
            ce = e.cursor()
            ce.execute('SELECT * from %s' % self.table_name)
            self.assertCountEqual(co.fetchall(), ce.fetchall())

    def test_create_id_map(self):
        n = create_id_map(self.id_map[1], self.id_map_obs)
        self.assertEqual(n, self.id_map[0])
        with connect(self.id_map_obs) as o, connect(self.id_map_exp) as e:
            co = o.cursor()
            co.execute('SELECT * from %s' % self.id_map_table)
            ce = e.cursor()
            ce.execute('SELECT * from %s' % self.id_map_table)
            self.assertCountEqual(co.fetchall(), ce.fetchall())

    def test_sort_uniref(self):
        sort_uniref(self.exp_db_fp, self.uniref_fp, self.tmp_dir)
        for fp in self.uniref_res:
            obs = join(self.tmp_dir, fp)
            exp = _get_data_dir()(fp)
            with open(obs) as o, open(exp) as e:
                self.assertEqual(o.read(), e.read())
            remove(obs)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)
        remove(self.obs_db_fp)
        remove(self.id_map_obs)


if __name__ == '__main__':
    main()
