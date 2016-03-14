# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click


@click.command()
@click.pass_context
def cli(ctx):
    '''Print micronota configuration details.

    It includes system info, database info, and external dependencies, etc.
    '''
    click.echo(repr(ctx.parent.config))
