# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import walk
from os.path import join, exists, expanduser, abspath, dirname, basename
from configparser import ConfigParser
from logging.config import fileConfig
from importlib.util import find_spec
from itertools import chain
from collections import OrderedDict

import click


class Configuration(object):
    '''micronota configuration container.

    Parameters
    ----------
    param_fp : str
        config file path for param.
    log_fp : str
        config file path for log.
    misc_fp : str
        config file path for other settings.

    Attributes
    ----------
    features : dict
        {tool: database}. Both key and value are string.
    cds : dict
        {tool: database}.
    param : dict of dict
        {tool: {param: value}}.
    db_dir : str
        database directory.
    db : dict
        database name and their abs path
    app_dir : str
        directory for micronota data files. It is different in
        different OS.
    '''
    def __init__(self, misc_fp=None, param_fp=None, log_fp=None):
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

        fp = join(self.app_dir, misc_fn)
        if misc_fp is not None:
            self._misc_fp = abspath(misc_fp)
        elif exists(fp):
            self._misc_fp = fp
        else:
            self._misc_fp = join(root_dir, dat_dir, misc_fn)
        self._set_misc_config(self._misc_fp)

        fp = join(self.app_dir, log_fn)
        if log_fp is not None:
            self._log_fp = abspath(log_fp)
        elif exists(fp):
            self._log_fp = fp
        else:
            self._log_fp = join(root_dir, dat_dir, log_fn)
        self._set_log_config(self._log_fp)

        default_fp = join(root_dir, dat_dir, param_fn)
        global_fp = join(self.app_dir, param_fn)
        self._param_fps = [default_fp, global_fp]
        if param_fp is not None:
            self._param_fps.append(param_fp)
        self._set_param_config(self._param_fps)

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
        tools = chain(self.param, self.features, self.cds)
        for i in tools:
            if i == 'DEFAULT':
                continue
            found = find_spec('.'.join([self._pkg, mod, i]))
            if found is None:
                raise NotImplementedError('%s not implemented.' % i)

    @property
    def db_dir(self):
        return self._db_dir

    @db_dir.setter
    def db_dir(self, db_dir):
        db = {}
        self._db_dir = db_dir
        for dirpath, dirnames, filenames in walk(db_dir):
            if not dirnames:
                db[basename(dirpath)] = dirpath
        self.db = db

    def __repr__(self):
        from sys import platform, version
        lines = []

        info = OrderedDict()
        info['system'] = {
            'OS': platform,
            'Python': version}

        info['micronota'] = OrderedDict([
            ('database directory', self.db_dir),
            ('global config directory', self.app_dir),
            ('general config file', self._misc_fp),
            ('log config file', self._log_fp),
            ('param config file', ','.join(self._param_fps))])

        info['databases'] = self.db

        info['features'] = dict()
        for tool in self.features:
            info['features'][tool] = self.features[tool]

        info['cds'] = dict()
        for tool in self.cds:
            info['cds'][tool] = self.cds[tool]

        info['parameters'] = dict()
        for tool in self.param:
            if tool == 'DEFAULT':
                continue
            l = []
            for opt in self.param[tool]:
                v = self.param[tool][opt]
                l.append('%s:%s' % (opt, v))
            info['parameters'][tool] = ' '.join(l)

        for k1 in info:
            lines.append(k1.upper())
            lines.append('=' * len(k1))
            info_ = info[k1]
            if info_:
                max_len = max([len(i) for i in info_])
                for k2 in info_:
                    if info_[k2] is None:
                        line = k2
                    else:
                        line = "{key:<{width}}: {value}".format(
                            width=max_len, key=k2, value=info_[k2])
                    lines.append(line)
            lines.append('\n')
        return '\n'.join(lines)
