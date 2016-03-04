# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import join
from importlib import import_module

import click

from ..cli import cmd, AliasedGroup
from .. import db


@cmd.group(cls=AliasedGroup)
@click.pass_context
def cli(ctx):
    '''Database operations.'''
    pass


@cli.command('prepare')
@click.argument('databases', nargs=-1)
@click.option('-d', '--cache_dir', required=True,
              type=click.Path(file_okay=False),
              help=('The directory to cache the downloaded files so that file '
                    'do not need to be downloaded again if it exists there.'))
@click.option('-f', '--force', is_flag=True,
              help='Force overwrite.')
@click.pass_context
def create_db(ctx, databases, cache_dir, force):
    '''Prepare database.

    Download the files for the specified DATABASES and manipulate
    them as proper format for micronota.'''
    # this cmd is 2-level nested, so double "parent"
    grandparent_ctx = ctx.parent.parent
    config = grandparent_ctx.config
    func_name = 'prepare_db'

    for d in databases:
        submodule = import_module('.%s' % d, db.__name__)
        f = getattr(submodule, func_name)
        out_d = join(config.db_dir, d)
        makedirs(out_d, exist_ok=True)
        f(out_d, cache_dir, force=force)
