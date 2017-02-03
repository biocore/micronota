# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from os.path import join
from logging import getLogger

import click

from ..database.uniprot import add_metadata


logger = getLogger(__name__)


@click.command(hidden=True)
@click.argument('infile', type=click.Path(), nargs=-1)
@click.argument('outpath', type=click.Path(),  nargs=1)
@click.pass_context
def cli(ctx, infile, outpath):
    '''Create UniProt protein cross-ref database.'''
    n = 0
    for fp in infile:
        if fp.endswith('.gz'):
            with gzip.open(fp) as f:
                n += add_metadata(f, outpath)
        else:
            with open(fp) as f:
                n += add_metadata(f, outpath)
    logger.info('Parsed %d records from %r' % (n, infile))
