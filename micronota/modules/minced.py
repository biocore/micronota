# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

from skbio.io import read


logger = getLogger(__name__)


def parse(fp='minced.gff'):
    '''Parse the annotation and add it to interval metadata.

    Parameters
    ----------
    fp : str
        file path from minced prediction

    '''
    logger.debug('Parsing minced prediction')
    return read(fp, format='gff3')
