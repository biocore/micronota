r'''
Utility functionality
=====================

.. currentmodule:: micronota.util

This module (:mod:`micronota.util`) provides various utility functionality,
including config config, unit-testing convenience function.

'''


# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from sys import platform, version
from os.path import join, expanduser
from configparser import ConfigParser


_HOME = expanduser('~')

_CONFIG_PATH = join(_HOME, '.micronota.config')


def _create_config():
    '''Return ``ConfigParser`` object.
    '''
    config = ConfigParser(allow_no_value=True,
                          strict=True)
    # Make the parser case sensitive; the default is not.
    config.optionxform = str
    # set the default key-value pairs
    config['DEFAULT']['db_path'] = join(_HOME, 'micronota_db')
    return config


def get_config_info(config):
    '''Return the micronota config info.'''
    info = dict()
    info['system'] = {
        'OS': platform,
        'python': version}
    info['micronota config'] = {
        'config file': _CONFIG_PATH,
        'database folder': config['DEFAULT']['db_path']}
    info['micronota database'] = list_db()
    return info


def list_db():
    '''Return database info.'''
    return dict()
