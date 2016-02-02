# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mktemp
from configparser import ConfigParser
from os.path import dirname

from micronota.util import parse_config


class ConfigParserTests(TestCase):
    def setUp(self):
        self.cfg_path = mktemp()
        self.cfg = {}

    def test_parse_config(self):
        pass


if __name__ == '__main__':
    main()
