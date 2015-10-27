# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import join
import urllib.request

_
def download_db():
    '''Download database.
    '''

    with urllib.request.urlopen('http://python.org/') as response:
        html = response.read()
    pass


