# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from configparser import ConfigParser

from skbio.util import get_data_path

from micronota.util import _create_config


class ConfigParserTests(TestCase):

    def test_create_config(self):
        cfg_exp = ConfigParser()
        cfg_exp['DEFAULT']['db_path'] = 'mn_db'

        cfg_obs = _create_config(get_data_path('default.cfg'))

        self.assertEqual(cfg_exp, cfg_obs)

    def test_create_config_overwrite(self):
        cfg_obs = _create_config(None)
        cfg_obs.read(get_data_path('default.cfg'))
        # overwrite "db_path"
        cfg_obs.read(get_data_path('param.cfg'))

        cfg_exp = ConfigParser()
        cfg_exp['DEFAULT']['db_path'] = 'db'
        cfg_exp.add_section('prodigal')
        cfg_exp['prodigal']['-t'] = '1'
        self.assertEqual(cfg_obs, cfg_exp)

if __name__ == '__main__':
    main()
