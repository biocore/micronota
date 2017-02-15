# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click
import os

from pkg_resources import resource_filename
import yaml

from ..workflow import integrate


@click.command()
@click.option('-o', '--out-dir', type=click.Path(file_okay=False),
              required=True,
              help='Output directory path.')
@click.option('--out-fmt', type=click.Choice(['gff3', 'genbank']),
              default='genbank',
              help='Output format for the annotated sequences.')
@click.option('--config', type=click.Path(exists=True, dir_okay=False),
              help='Config file for annotation workflow.')
@click.pass_context
def cli(ctx, out_dir, out_fmt, config):
    '''Integrate annotations into final output.'''
    if config is None:
        config = resource_filename('micronota', 'config.yaml')
    with open(config) as fh:
        cfg = yaml.load(fh)
    cfg['mode'] = 'draft'
    for fn in os.listdir(out_dir):
        if fn.endswith('.valid.fna'):
            integrate(cfg, out_dir, fn, out_fmt)
            break
