# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, splitext, basename
import click

from ..workflow import annotate


@click.command()
@click.option('-i', '--input_fp', type=str,
              required=True,
              help='Input file path.')
@click.option('--in_fmt', type=click.Choice(['fasta', 'genbank']),
              default='fasta',
              help='The format of input file.')
@click.option('-o', '--output_dir', type=str,
              required=True,
              help='Output directory path.')
@click.option('--out_fmt', type=click.Choice(['gff3', 'genbank']),
              default='genbank',
              help='Output format for the annotated sequences.')
@click.option('--cpus', type=int, default=1,
              help='Number of CPUs to use.')
@click.option('--kingdom',
              type=click.Choice(['Bacteria', 'Archaea', 'Viruses']),
              default='Bacteria',
              help='Kingdom of the input sequence organism.')
@click.pass_context
def cli(ctx, input_fp, in_fmt, output_dir, out_fmt, cpus, kingdom):
    '''Annotate prokaryotic genomes.'''
    annotate(input_fp, in_fmt, output_dir, out_fmt, cpus, ctx.parent.config)
