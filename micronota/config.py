# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, exists, expanduser, abspath, dirname
from configparser import ConfigParser
from logging.config import fileConfig
from importlib.util import find_spec
from itertools import chain

import click


class Configuration(object):
    def __init__(self, misc_fp=None, param_fp=None, log_fp=None):
        '''
        Parameters
        ----------
        param_fp : str
            config file path for param
        log_fp : str
            config file path for log
        misc_fp : str
            config file path for other
        '''
        self._pkg = 'micronota'
        self.features = iter({})
        self.cds = iter({})
        root_dir = abspath(dirname(__file__))
        # where global conf files are
        self.app_dir = click.get_app_dir(self._pkg)
        # where default conf files are
        dat_dir = 'support_files'

        misc_fn = 'misc.cfg'
        log_fn = 'log.cfg'
        param_fn = 'param.cfg'

        default_fp = join(root_dir, dat_dir, misc_fn)
        global_fp = join(self.app_dir, misc_fn)
        if misc_fp is not None:
            self._set_misc_config(misc_fp)
        elif exists(global_fp):
            self._set_misc_config(global_fp)
        else:
            self._set_misc_config(default_fp)

        default_fp = join(root_dir, dat_dir, log_fn)
        global_fp = join(self.app_dir, log_fn)
        if log_fp is not None:
            self._set_log_config(log_fp)
        elif exists(global_fp):
            self._set_log_config(global_fp)
        else:
            self._set_log_config(default_fp)

        default_fp = join(root_dir, dat_dir, param_fn)
        global_fp = join(self.app_dir, param_fn)
        fps = [default_fp, global_fp]
        if param_fp is not None:
            fps.append(param_fp)
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
        mod = 'bfillings'
        tools = chain(self.param, self.features)
        for i in tools:
            if i == 'DEFAULT':
                continue
            found = find_spec('.'.join([self._pkg, mod, i]))
            if found is None:
                raise NotImplementedError('%s not implemented.' % i)

    def __repr__(self):
        from sys import platform, version
        lines = []

        info = dict()
        info['system'] = {
            'OS': platform,
            'Python': version}
        info['micronota'] = {
            'db directory': self.db_dir,
            'config directory': self.app_dir}

        info['parameters'] = dict()
        for tool in self.param:
            if tool == 'DEFAULT':
                continue
            l = []
            for opt in self.param[tool]:
                v = self.param[tool][opt]
                l.append('%s:%s' % (opt, v))
            info['parameters'][tool] = ' '.join(l)

        info['features'] = dict()
        for tool in self.features:
            info['features'][tool] = self.features[tool]

        info['cds'] = dict()
        for tool in self.cds:
            info['cds'][tool] = self.cds[tool]

        for k1 in sorted(info):
            lines.append(k1)
            lines.append('=' * len(k1))
            info_ = info[k1]
            max_len = max([len(i) for i in info_])
            for k2 in sorted(info_):
                line = "{key:<{width}}: {value}".format(
                    width=max_len, key=k2, value=info_[k2])
                lines.append(line)
            lines.append('\n')
        return '\n'.join(lines)
