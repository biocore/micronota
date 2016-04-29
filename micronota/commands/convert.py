# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click

from ..util import convert


@click.command()
@click.option('-i', '--in_fmt',
              required=True,
              help='Input file format.')
@click.option('-o', '--out_fmt',
              required=True,
              help='Output file path.')
@click.argument('input', type=click.File('r'), nargs=1)
@click.argument('output', type=click.File('w'), nargs=1)
@click.pass_context
def cli(ctx, in_fmt, out_fmt, in_f, out_f):
    '''Format conversion.'''
    convert(in_fmt, out_fmt, in_f, out_f)
