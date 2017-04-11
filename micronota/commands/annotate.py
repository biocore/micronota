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
@click.option('-i', '--in-seqs', type=click.Path(exists=True, dir_okay=False),
              required=True,
              help='Input sequence file (can be gzip file.')
@click.option('--in-fmt', type=click.Choice(['fasta', 'genbank', 'gff3']),
              default='fasta',
              help='The format of input file. If it is gff3 format, it must contain seq in it.')
@click.option('--min-len', type=int, default=500,
              help='Min seq length. Input sequence shorter than this will be filtered out.')
@click.option('-o', '--out-dir', type=click.Path(file_okay=False),
              required=True,
              help='Output directory.')
@click.option('--out-fmt', type=click.Choice(['gff3', 'genbank']),
              default='genbank',
              help='Output format for the annotated sequences.')
@click.option('--gcode', type=int, default=11,
              help='Genetic code to predict ORFs.')
@click.option('--kingdom', type=click.Choice(['bacteria', 'archaea']), default='bacteria',
              required=True,
              help='which Kingdom the sequences are from')
@click.option('--mode', type=click.Choice(['finished', 'draft', 'metagenome']),
              default='draft',
              help='Run the proper mode to annotate the input sequences '
                   '(finished genome, draft genome, or metagenome.')
@click.option('--task', type=str, nargs=-1,
              default=None,
              help='what annotation task(s) to run?')
@click.option('--cpu', type=int, default=1,
              help='Number of CPUs to use.')
@click.option('--force', is_flag=True, default=False,
              help='Force overwrite if the output directory exists')
@click.option('-d', '--dry-run', is_flag=True,
              help='Do not execute anything.')
@click.option('--config', type=click.Path(exists=True, dir_okay=False),
              help='Config file for annotation workflow.')
@click.pass_context
def cli(ctx, in_seqs, in_fmt, min_len, out_dir, out_fmt, gcode, kingdom, mode, task,
        cpu, force, dry_run, config):
    '''Annotate prokaryotic sequences.'''
    annotate(in_seqs, in_fmt, min_len,
             out_dir, out_fmt,
             gcode, kingdom, mode, task,
             cpu, force, dry_run, config)
