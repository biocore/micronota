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
from os import remove
from os.path import join, expanduser, exists
from configparser import ConfigParser
from contextlib import contextmanager
from ftplib import FTP
from urllib.parse import urlparse


_HOME = expanduser('~')

_CONFIG_PATH = join(_HOME, '.micronota.config')


def _create_config(fp):
    '''Return ``ConfigParser`` object.
    '''
    config = ConfigParser(allow_no_value=True,
                          strict=True)
    # Make the parser case sensitive; the default is not.
    config.optionxform = str

    # 1. set the default key-value pairs
    config['DEFAULT']['db_path'] = join(_HOME, 'micronota_db')

    # 2. read the default config file
    if exists(_CONFIG_PATH):
        config.read(_CONFIG_PATH)

    # 3. read the provided config file
    if fp:
        config.read(fp)

    return config


def get_config_info(config):
    '''Return the micronota config info.

    Parameters
    ----------
    config : ConfigParser

    Returns
    -------
    dict
        all the information about micronota setup.
    '''
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


def _overwrite_file(fp, overwrite=False, append=True):
    if exists(fp):
        if overwrite:
            remove(fp)
        elif append:
            return
        else:
            raise FileExistsError(
                'The file %s exists. User "force" to overwrite it.' % fp)


class FTPHandler(object):
    def __init__(self, parsed_url):
        ''''''
        self._ftp = FTP(parsed_url.netloc)
        self._ftp.login(parsed_url.username, parsed_url.password)


    def download(self, src, dst):
        ''''''

    def close(self):
        ''''''
        self._ftp.close()


@contextmanager
def get_remote_handler(url):
    url_parsed = urlparse(url)
    if url_parsed.scheme == 'ftp':
        co = FTPHandler(url_parsed)
    else:
        raise Exception('Invalid protocol: %s' % url_parsed.scheme)
    yield co
    co.close()

