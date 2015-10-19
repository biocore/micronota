# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


def _get_parameter(constructor, s, prefix='-', **kwargs):
    name = s.lstrip(prefix)
    i = len(s) - len(name)
    return constructor(Prefix=s[:i], Name=s[i:], **kwargs)
