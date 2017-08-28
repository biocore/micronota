# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click

from ..util import check_seq


@click.command()
@click.option('-i', '--in-file',
              type=click.Path(exists=True, dir_okay=False),
              required=True,
              help='Input file path.')
@click.option('--in-fmt', type=click.Choice(['fasta', 'genbank', 'gff3']),
              default='fasta',
              help='The format of input file. If it is gff3 format, it must contain seq in it.')
@click.option('-o', '--out-file',
              required=True,
              type=click.Path(exists=False, dir_okay=False),
              help='Output file path.')
@click.option('--out-fmt', type=click.Choice(['fasta', 'genbank', 'gff3']),
              default='fasta',
              help='The format of output file.')
@click.option('-l', '--length',
              type=int,
              default=500,
              required=True,
              help='Minimum seq length to keep')
@click.pass_context
def cli(ctx, in_file, in_fmt, out_file, out_fmt, length):
    '''Filter and validate input sequences.'''
    with open(out_file, 'w') as out:
        for seq in filter_seq(in_file, in_fmt, lambda s: len(s) < length):
            write(seq, format='fasta', into=out)

