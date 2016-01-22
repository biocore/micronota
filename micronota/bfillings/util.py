# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial
from os.path import join, basename, splitext
from inspect import stack

from skbio.util import get_data_path


def _get_parameter(constructor, s, prefix='-', **kwargs):
    name = s.lstrip(prefix)
    i = len(s) - len(name)
    return constructor(Prefix=s[:i], Name=s[i:], **kwargs)


def _get_data_dir():
    # get caller's file name
    caller_fn = basename(stack()[1][1])
    # remove file suffix and prefix of "test_"
    d = splitext(caller_fn)[0][5:]
    return partial(get_data_path, subfolder=join('data', d))
