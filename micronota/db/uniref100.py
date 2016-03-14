# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

from ._uniref import _prepare


def prepare_db(downloaded, out_d='uniref',
               uniref_url='ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz',
               force=False):
    logger = getLogger(__name__)
    logger.info('Preparing UniRef100 database')

    _prepare(downloaded, out_d, uniref_url, 100, force)
