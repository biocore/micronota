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
@click.option('-i', '--in-seq', type=click.Path(exists=True, dir_okay=False),
              required=True,
              help='Input sequence file (can be gzip file.')
@click.option('-o', '--out-file', type=click.Path(exists=False, dir_okay=False),
              required=True,
              help='Output annotation file.')
@click.option('-d', '--annot-dir', type=click.Path(file_okay=False),
              required=True,
              help='directory that has the outputs of annotation tools.')
@click.option('--out-fmt', type=click.Choice(['gff3', 'genbank']),
              default='gff3',
              help='Output format for the annotation file.')
@click.option('--protein-xref', type=click.Path(exists=True, dir_okay=False),
              help='sqlite file that stores protein cross-ref info.')
@click.pass_context
def cli(ctx, in_seq, out_file, annot_dir, out_fmt, protein_xref):
    '''Integrate annotations into final output.'''
    print(protein_xref is None)
    integrate(in_seq, annot_dir, protein_xref, out_file,  out_fmt)

