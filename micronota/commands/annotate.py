# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click

from ..workflow import annotate


def validate_gcode(ctx, param, value):
    try:
        rolls, dice = map(int, value.split('d', 2))
        return (dice, rolls)
    except ValueError:
        raise click.BadParameter('rolls need to be in format NdM')


@click.command()
@click.option('-i', '--in-seq', type=click.Path(exists=True, dir_okay=False),
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
@click.option('--gcode', type=int, default=None,
              help='Genetic code to predict ORFs. Default value depends on Kingdom. '
                   '11 for bacteria and archea; 1 for eukarya')
@click.option('--kingdom', type=click.Choice(['bacteria', 'archaea', 'eukarya']), default='bacteria',
              required=True,
              help='which Kingdom the sequences are from')
@click.option('--mode', type=click.Choice(['finished', 'draft', 'metagenome']),
              default='draft',
              help='Run the proper mode to annotate the input sequences '
                   '(finished genome, draft genome, or metagenome.')
@click.option('--cpu', type=int, default=1,
              help='Number of CPUs to use.')
@click.option('-f', '--force', is_flag=True, default=False,
              help='Force overwrite if the output directory exists')
@click.option('-d', '--dry-run', is_flag=True,
              help='Do not execute anything.')
@click.option('--quality', type=bool, default=False,
              help='whether to compute the quality score for the sequence/annotation')
@click.option('--config', type=click.Path(exists=True, dir_okay=False),
              help='yaml file to config annotation workflow.')
@click.argument('task', type=str, nargs=-1)
@click.pass_context
def cli(ctx, in_seq, in_fmt, min_len, out_dir, out_fmt, gcode, kingdom, mode, task,
        cpu, force, dry_run, quality, config):
    '''Annotate genomic sequences.'''
    if gcode is None:
        if kingdom == 'eukarya':
            gcode = 1
        else:
            gcode = 11

    annotate(in_seq, in_fmt, min_len,
             out_dir, out_fmt,
             gcode, kingdom, mode, task,
             cpu, force, dry_run, quality, config)
