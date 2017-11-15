# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from logging import getLogger

import click

from ..database.rfam import filter_models


logger = getLogger(__name__)


# keep this command hidden from help msg
# this option is only available in click v7
@click.command()
@click.option('--operation', type=click.Choice(['bacteria', 'archaea', 'eukarya', 'default']),
              default='default', required=True,
              help='Keep bacteria/archaea/eukarya rRNA models or filter away tRNA, tmRNA, rRNA models (default)')
@click.argument('infile', type=click.File('r'), nargs=1)
@click.argument('outfile', type=click.Path(),  nargs=1, default='rfam_filtered.cm')
@click.pass_context
def cli(ctx, operation, infile, outfile):
    '''Create rfam database for micronota usage.'''
    kingdom_models = {'bacteria': {('RF00001', '5S_rRNA'),
                                   ('RF00177', 'SSU_rRNA_bacteria'),
                                   ('RF02541', 'LSU_rRNA_bacteria')},
                      'archaea': {('RF00001', '5S_rRNA'),
                                  ('RF01959', 'SSU_rRNA_archaea'),
                                  ('RF02540', 'LSU_rRNA_archaea')},
                      'eukarya': {('RF00001', '5S_rRNA'),
                                  ('RF00002', '5_8S_rRNA'),
                                  ('RF01960', 'SSU_rRNA_eukarya'),
                                  ('RF02543', 'LSU_rRNA_eukarya')}}
    if operation == 'default':
        with open(join(outfile), 'w') as out:
            filter_models(infile, out)
    else:
        with open(outfile, 'w') as out:
            filter_models(infile, out, negate=True, models=kingdom_models[kingdom])
