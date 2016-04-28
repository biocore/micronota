# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from dumpling import OptionParam, ArgmntParam


_scan_params = [
    OptionParam('--tblout', name='out', help='save parseable table of hits to file'),
    # set default to 1 instead of all available cores.
    OptionParam('--cpu', name='cpus', value=1, help='number of parallel CPU workers to use for multithreads'),
    ArgmntParam(name='db', help='HMM/CM database file'),
    ArgmntParam(name='query', help='input sequence to scan')]
