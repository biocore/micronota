#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mktemp
from unittest import TestCase, main
from os.path import dirname

from micronota.bfillings.util import _get_data_dir
from micronota.db.tigrfam import prepare_metadata


class TigrfamTests(TestCase):
    def setUp(self):
        self.obs_db_fp = mktemp()
        self.exp_db_fp = _get_data_dir()('tigrfam.db')
        self.d = dirname(self.exp_db_fp)

    def test_prepare_metadata(self):
        prepare_metadata(self.d, self.obs_db_fp)
        with open(self.obs_db_fp, 'rb') as o, open(self.exp_db_fp, 'rb') as e:
            self.assertEqual(o.read(), e.read())

if __name__ == '__main__':
    main()
