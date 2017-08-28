from skbio import read, write
import click

from ..util import filter_partial_genes


@click.command()
@click.option('-i', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('-o', type=click.Path(exists=False, dir_okay=False), required=True)
def cli(i, o):
    filter_partial_genes(i, o)

