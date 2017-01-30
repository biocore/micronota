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


def parse(fp='prodigal.txt'):
    '''Parse the annotation and add it to interval metadata.

    Parameters
    ----------
    fn : str
        the file name from prodigal prediction

    Yield
    -----
    tuple of str and IntervalMetadata
        seq_id and interval metadata
    '''
    logger.debug('Parsing prodigal prediction')
    return read(fp, format='gff3')
