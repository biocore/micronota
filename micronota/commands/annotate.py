# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import basename, join
import os

import click

from ..workflow import annotate, validate_seq, integrate


@click.command()
@click.option('-i', '--input_fp', type=click.Path(exists=True, dir_okay=False),
              required=True,
              help='Input file path.')
@click.option('--in_fmt', type=click.Choice(['fasta', 'genbank', 'gff3']),
              default='fasta',
              help='The format of input file. If it is gff3 format, it must contain seq in it.')
@click.option('--min_len', type=int, default=500,
              help='Min seq length. Input sequence shorter than this will be filtered out.')
@click.option('-o', '--output_dir', type=click.Path(file_okay=False),
              required=True,
              help='Output directory path.')
@click.option('--out_fmt', type=click.Choice(['gff3', 'genbank']),
              default='genbank',
              help='Output format for the annotated sequences.')
@click.option('--gcode', type=int, default=11,
              help='Genetic code to predict ORFs.')
@click.option('--cpus', type=int, default=1,
              help='Number of CPUs to use.')
@click.option('--force', is_flag=True, default=False,
              help='Force overwrite if the output directory exists')
@click.option('-d', '--dry_run', is_flag=True,
              help='Do not execute anything.')
@click.option('--config', type=click.Path(exists=True, dir_okay=False),
              help='Config file for annotation workflow.')
@click.pass_context
def cli(ctx, input_fp, in_fmt, min_len, output_dir, out_fmt, gcode,
        cpus, force, dry_run, config):
    '''Annotate prokaryotic genomes.'''
    os.makedirs(output_dir, exist_ok=True)

    seq_fn = basename(input_fp)

    validate_seq(input_fp, in_fmt, min_len, join(output_dir, seq_fn))

    done = annotate(seq_fn, output_dir, gcode, cpus, force, dry_run, config)

    if done:
        # if snakemake finishes successfully
        integrate(output_dir, seq_fn, out_fmt)
