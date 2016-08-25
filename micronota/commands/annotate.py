# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click

from ..workflow import annotate


@click.command()
@click.option('-i', '--input_fp', type=click.Path(exists=True, dir_okay=False),
              required=True,
              help='Input file path.')
@click.option('--in_fmt', type=click.Choice(['fasta', 'genbank', 'gff3']),
              default='fasta',
              help='The format of input file. If it is gff3 format, it must contain seq in it.')
@click.option('-o', '--output_dir', type=click.Path(file_okay=False),
              required=True,
              help='Output directory path.')
@click.option('--out_fmt', type=click.Choice(['gff3', 'genbank']),
              default='genbank',
              help='Output format for the annotated sequences.')
@click.option('--cpus', type=int, default=1,
              help='Number of CPUs to use.')
@click.option('--force', is_flag=True,
              help='Force overwrite if the output directory exists')
@click.option('-d', '--dry_run', is_flag=True,
              help='Do not execute anything.')
@click.pass_context
def cli(ctx, input_fp, in_fmt, output_dir, out_fmt,
        cpus, force, dry_run):
    '''Annotate prokaryotic genomes.'''
    annotate(input_fp, in_fmt, output_dir, out_fmt,
             cpus, force, dry_run, ctx.parent)
