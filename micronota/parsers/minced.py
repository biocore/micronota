# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join

from skbio.io import read
from skbio.metadata import IntervalMetadata


def parse(interval_metadata, out_dir, seq_fn):
    '''Parse the annotation and add it to interval metadata.

    Parameters
    ----------
    interval_metadata : dict
        key is seq id and value is the interval metadata for the seq.
    out_dir : str
        the output dir
    seq_fn : str
        input seq file name. With `out_dir` and `seq_fn`, you should
        be able to create the file name output from this prediction
        tool to parse the annotaton.

    '''
    fp = join(out_dir, 'minced', '%s.gff' % seq_fn)
    list(read(fp, format='gff3', interval_metadata=interval_metadata))
