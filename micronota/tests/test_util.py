# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, mock
from configparser import ConfigParser

from skbio.util import get_data_path

from micronota.util import _create_config


class ConfigParserTests(TestCase):
    def setUp(self):
        self.cfg_fps = list(
            map(get_data_path,
                ['default.cfg',  # test default setting
                 'global.cfg',   # test ~/.micronota.conf
                 'local.cfg']))  # test --config

        # mock the _CONFIG_PATH
        self.patcher1 = mock.patch('micronota.util._CONFIG_PATH',
                                   self.cfg_fps[1])
        self.patcher2 = mock.patch('micronota.util._DB_PATH',
                                   '{d}')
        # the mocks will be effective in all tests.
        self.patcher1.start()
        self.patcher2.start()

    # re-mock to not read global config file
    @mock.patch('micronota.util._CONFIG_PATH', '')
    def test_create_config_default(self):
        cfg_exp = ConfigParser()
        cfg_exp.read(self.cfg_fps[0])

        cfg_obs = _create_config(None)

        self.assertEqual(cfg_exp, cfg_obs)

    def test_create_config_global(self):
        '''test overwriting by global cfg.'''
        cfg_obs = _create_config(None)

        cfg_exp = ConfigParser()
        cfg_exp.read(self.cfg_fps[:2])
        self.assertEqual(cfg_obs, cfg_exp)

    def test_create_config_local(self):
        '''test overwriting by local cfg.'''
        cfg_obs = _create_config(self.cfg_fps[2])

        cfg_exp = ConfigParser()
        cfg_exp.read(self.cfg_fps[:3])
        self.assertEqual(cfg_obs, cfg_exp)

    def tearDown(self):
        self.patcher1.stop()
        self.patcher2.stop()


if __name__ == '__main__':
    main()
