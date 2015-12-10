# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click

from ..cli import cmd, AliasedGroup


@cmd.group(cls=AliasedGroup)
@click.pass_context
def cli(ctx):
    '''Database operations.'''
    pass


@cli.command('update')
@click.option('-i', '--input', type=click.STRING)
@click.option('-o', '--output', type=click.STRING)
@click.pass_context
def update_db():
    '''Update database.'''
    pass
