# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, mock
from configparser import ConfigParser
from os.path import dirname

from skbio.util import get_data_path

from micronota.config import Configuration


class ConfigurationTests(TestCase):
    def setUp(self):
        self.cfg_fps = list()
        self.misc_fp = get_data_path('misc.cfg')
        self.misc_fp_local = get_data_path('misc_local.cfg')
        self.param_fp = get_data_path('param.cfg')
        self.param_fp_local = get_data_path('param_local.cfg')
        self.dir = dirname(self.misc_fp)

    def test_global_config(self):
        '''Test the global setting override the default.'''
        with mock.patch('click.get_app_dir', return_value=self.dir):
            obs = Configuration()
            exp = ConfigParser(allow_no_value=True)
            exp.read(self.misc_fp)
            self.assertEqual(exp['general']['db_dir'], obs.db_dir)
            self.assertEqual(exp['feature'], obs.features)
            self.assertEqual(exp['cds'], obs.cds)
            exp = ConfigParser(allow_no_value=True)
            exp.read(self.param_fp)
            self.assertEqual(exp, obs.param)

    def test_local_config(self):
        '''Test the local settings override the default and global.'''
        with mock.patch('click.get_app_dir', return_value=self.dir):
            obs = Configuration(misc_fp=self.misc_fp_local,
                                param_fp=self.param_fp_local)
            exp = ConfigParser(allow_no_value=True)
            exp.read(self.misc_fp_local)
            self.assertEqual(exp['general']['db_dir'], obs.db_dir)
            self.assertEqual(exp['feature'], obs.features)
            self.assertEqual(exp['cds'], obs.cds)
            exp = ConfigParser(allow_no_value=True)
            exp.read([self.param_fp, self.param_fp_local])
            self.assertEqual(exp, obs.param)


if __name__ == '__main__':
    main()
