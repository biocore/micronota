# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import abspath, dirname, join, expanduser
from configparser import ConfigParser
from importlib.util import find_spec
import logging
import re

from . import app_dir


class ConfigParsingError(Exception):
    ''''''
    pass


class Configuration(object):
    def __init__(self, local_fp):
        ''''''
        data_dir = 'support_files'
        config_fp = 'micronota.cfg'
        default_fp = join(dirname(abspath(__name__)), data_dir, config_fp)
        global_fp = join(app_dir, config_fp)
        workflow, param, log = self._read([default_fp, global_fp, local_fp])
        self._tools = []
        self.workflow = self._set_workflow(workflow)
        self.param = self._set_param(param)
        self._set_log(log)
        self._check_avail()

    @staticmethod
    def _read(fps):
        '''Return a list of 3 lists of strings.

        The order in fps matters. The later cfg can override the previous.'''
        s = [[], [], []]
        pattern = re.compile(r'^#{9,} +\w+ +#{9,}\n', re.MULTILINE)
        for fp in fps:
            with open(fp) as f:
                s = f.read()
                items = re.split(pattern, s)
                try:
                    _, wf, param, log = items
                except ValueError:
                    raise ConfigParsingError(
                        'There should be 3 parts in the config file')
                s[0].append(wf)
                s[1].append(param)
                s[2].append(log)

    @staticmethod
    def _read_config(ss):
        '''Read in list of config strings.'''
        config = ConfigParser()
        for s in ss:
            config.read_string(s)
        return config

    def _set_workflow(self, s):
        '''read in the config for workflow.

        also check if the implementation is available.
        '''
        config = self._read_config(s)
        config['general']['db_path'] = expanduser(config['general']['db_path'])
        secs = ['feature', 'cds']
        for s in secs:
            sec = config[s]
            for k in sec:
                self._tools.append(k)
                try:
                    v = sec.getboolean(k)
                    sec[k] = v
                except ValueError:
                    v = sec[k]
        return config

    def _set_param(self, s):
        config = self._read_config(s)
        for sec in config:
            if sec not in self._tools:
                self._tools.append(sec)
        return config

    def _set_log(self, s):
        config = self._read_config(s)
        logging.config.fileConfig(config)

    def _check_avail(self):
        mod = dirname(__name__)
        submod = 'bfillings'
        for i in self._tools:
            found = find_spec('.'.join([mod, submod, i]))
            if found is None:
                raise NotImplementedError('%s not implemented.' % i)
