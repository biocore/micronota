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


from os.path import join, expanduser
from configparser import ConfigParser


_HOME = expanduser('~')


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
