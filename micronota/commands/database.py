# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkgutil
import importlib

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
@click.option('-d', '--cache_dir', type=str, required=True,
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
    verbose = grandparent_ctx.params['verbose']
    func_name = 'prepare_db'
    out_d = config['DEFAULT']['db_path']
    if not os.path.exists(out_d):
        os.mkdir(out_d)
    if not databases:
        databases = []
        for importer, modname, ispkg in pkgutil.iter_modules(db.__path__):
            databases.append(modname)
    for d in databases:
        if verbose > 0:
            click.echo('Start creating %s database...' % d)
        submodule = importlib.import_module('.%s' % d, db.__name__)
        f = getattr(submodule, func_name)
        f(out_d, cache_dir, force=force)
    if verbose > 0:
        click.echo('Finished creating databases')
