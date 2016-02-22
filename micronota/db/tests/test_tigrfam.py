# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkdtemp
from unittest import main
from os.path import dirname, join
from shutil import rmtree

from micronota.util import _DBTest
from micronota.bfillings.util import _get_data_dir
from micronota.db.tigrfam import prepare_db


class TigrfamTests(_DBTest):
    def setUp(self):
        self.tmp_dir = mkdtemp()

        self.db_fp = 'tigrfam_v15.0.db'
        self.obs_db_fp = join(self.tmp_dir, self.db_fp)
        self.exp_db_fp = _get_data_dir()(self.db_fp)
        self.d = dirname(self.exp_db_fp)

        self.pressed_fp = 'tigrfam_v15.0.hmm.h3f'
        self.obs_pressed_fp = join(self.tmp_dir, self.pressed_fp)
        self.exp_pressed_fp = _get_data_dir()(self.pressed_fp)

    def test_prepare_db(self):
        prepare_db(self.tmp_dir, self.d)
        self._test_eq_db(self.obs_db_fp, self.exp_db_fp)
        with open(self.obs_pressed_fp, 'rb') as o:
            with open(self.exp_pressed_fp, 'rb') as e:
                self.assertEqual(o.read(), e.read())

    def test_prepare_db_not_overwrite(self):
        with self.assertRaisesRegex(
                FileExistsError, r'The file .* exists.'):
            prepare_db(self.d, self.d)

    def tearDown(self):
        rmtree(self.tmp_dir)

if __name__ == '__main__':
    main()
