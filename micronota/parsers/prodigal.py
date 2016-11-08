# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from logging import getLogger

from skbio.io import read


logger = getLogger(__name__)


def parse(interval_metadata, out_dir):
    '''Parse the annotation and add it to interval metadata.

    Parameters
    ----------
    interval_metadata : dict
        key is seq id and value is the interval metadata for the seq.
    out_dir : str
        the output dir
    '''
    logger.debug('Parsing prodigal prediction')
    fp = join(out_dir, 'prodigal.gff')
    list(read(fp, format='gff3', interval_metadata_dict=interval_metadata))
