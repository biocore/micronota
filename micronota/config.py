# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, exists, expanduser, abspath, dirname
from configparser import ConfigParser
from importlib.util import find_spec
from logging.config import fileConfig

import click


class Configuration(object):
    def __init__(self, param_fp=None, log_fp=None, misc_fp=None):
        ''''''
        self.features = iter({})
        self.cds = iter({})
        root_dir = abspath(dirname(__name__))
        # where global conf files are
        app_dir = click.get_app_dir(root_dir)
        # where default conf files are
        dat_dir = 'support_files'

        misc_fn = 'misc.cfg'
        log_fn = 'log.cfg'
        param_fn = 'param.cfg'

        default_fp = join(root_dir, dat_dir, misc_fn)
        global_fp = join(app_dir, misc_fn)
        if misc_fp is not None:
            self._set_config(misc_fp)
        elif exists(global_fp):
            self._set_misc_config(global_fp)
        else:
            self._set_misc_config(default_fp)

        default_fp = join(root_dir, dat_dir, log_fn)
        global_fp = join(app_dir, log_fn)
        if log_fp is not None:
            self._set_log_config(log_fp)
        elif exists(global_fp):
            self._set_log_config(global_fp)
        else:
            self._set_log_config(default_fp)

        default_fp = join(root_dir, dat_dir, param_fn)
        global_fp = join(app_dir, param_fn)
        fps = [default_fp, global_fp, param_fp]
        self._set_param_config(fps)

        # check all specified tools are wrapped
        self._check_avail()

    @staticmethod
    def _read_config(f, **kwargs):
        config = ConfigParser(**kwargs)
        config.read(f)
        return config

    def _set_misc_config(self, fp):
        '''read in the config
        '''
        config = self._read_config(fp, allow_no_value=True)
        self.db_dir = expanduser(config['general']['db_dir'])
        if 'feature' in config:
            self.features = config['feature']
        if 'cds' in config:
            self.cds = config['cds']

    def _set_param_config(self, fps):
        config = self._read_config(fps, allow_no_value=True)
        self.param = config

    def _set_log_config(self, fp):
        fileConfig(fp)

    def _check_avail(self):
        mod = 'micronota'
        submod = 'bfillings'
        for i in self._tools:
            found = find_spec('.'.join([mod, submod, i]))
            if found is None:
                raise NotImplementedError('%s not implemented.' % i)
