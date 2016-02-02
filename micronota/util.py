r'''
Utility functionality
=====================

.. currentmodule:: micronota.util

This module (:mod:`micronota.util`) provides various utility functionality,
including config parser, unit-testing convenience function.

'''


# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from os.path import expanduser, join
from configparser import ConfigParser


_config_path = join(expanduser("~"), '.micronota.config')


def parse_config(fp=_config_path):
    '''
    '''
    cfg = ConfigParser.read(fp)
    return cfg
