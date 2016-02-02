# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click
from os import listdir
from os.path import abspath, join, dirname, splitext, exists

from .util import _create_config, _HOME


_CONFIG_PATH = join(_HOME, '.micronota.config')

_CONTEXT_SETTINGS = dict(
    # allow case insensitivity for the (sub)commands and their options
    token_normalize_func=lambda x: x.lower(),
    # set --help option for all (sub)commands
    help_option_names=['-h', '--help'])


class AliasedGroup(click.Group):
    '''Custom subclass of click.Group to enable alias for commands.

    This implements a subclass of click.Group that accepts a prefix
    for a command. If there were a (sub)command called "push", it would
    accept "pus" as an alias (as long as it is unique).

    This is borrowed from `click` example of alias.
    '''
    def _get_command(self, ctx, cmd_name):
        return click.Group.get_command(self, ctx, cmd_name)

    def get_command(self, ctx, cmd_name):
        # allow automatic abbreviation of the command.  "status" for
        # instance will match "st".  We only allow that however if
        # there is only one command.
        matches = [x for x in self.list_commands(ctx)
                   if x.lower().startswith(cmd_name.lower())]
        if not matches:
            return
        elif len(matches) == 1:
            return self._get_command(ctx, matches[0])
        ctx.fail('Too many matches: %s' % ', '.join(sorted(matches)))


class ComplexCLI(AliasedGroup):
    '''Custom subclass to load subcommands dynamically from a plugin folder.

    It looks in `commands` folder for subcommands.

    This is borrowed from `click` examples of complex.
    '''
    def list_commands(self, ctx):
        rv = []
        cmd_folder = abspath(join(dirname(__file__), 'commands'))
        for filename in listdir(cmd_folder):
            if filename.endswith('.py') and filename != '__init__.py':
                rv.append(splitext(filename)[0])
        rv.sort()
        return rv

    def _get_command(self, ctx, cmd_name):
        try:
            mod = __import__('micronota.commands.' + cmd_name,
                             None, None, ['cli'])
        except ImportError as err:
            print(err)
            return
        return mod.cli


@click.group(cls=ComplexCLI, context_settings=_CONTEXT_SETTINGS)
@click.option('-v', '--verbose', count=True,
              help=("Verbosity. You can use multiple v's, "
                    "e.g. -vv, to gradually increase verbosity."))
@click.option('-d', '--debug', is_flag=True,
              help='Debug mode. Save intermediate files. ')
@click.option('--info', is_flag=True,
              help=('Print micronota configuration details, database info, '
                    'and external dependencies.'))
@click.option('--config', default=None,
              help='The config file.')
@click.version_option()   # add --version option
@click.pass_context
def cmd(ctx, debug, verbose, info, config):
    '''Annotation pipeline for Bacterial and Archaeal (meta)genomes.

    It predicts features (ncRNA, coding genes, etc.) on the input sequences
    and assign functions to those features.

    It works best on long sequences, such as assembled contigs, draft genomes,
    or full genomes.

    For more info, please check out https://github.com/biocore/micronota.
    '''
    # read the default config file if it exist.
    ctx.config = _create_config()
    if exists(_CONFIG_PATH):
        ctx.config.read(_CONFIG_PATH)
    # if config file is provided, also read it in and it will overwrite
    # the previous settings.
    if config is not None:
        ctx.config.read(config)
