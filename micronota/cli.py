# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click
from os import listdir
from os.path import abspath, join, dirname, splitext


_CONTEXT_SETTINGS = dict(
    # allow case insensitivity for the command options
    token_normalize_func=lambda x: x.lower(),
    # set --help option for all (sub)commands
    help_option_names=['-h', '--help'])


class AliasedGroup(click.Group):
    '''Custom subclass of click.Group to enable alias for commands.

    This implements a subclass of click.Group that accepts a prefix
    for a command. If there were a (sub)command called "push", it would
    accept "pus" as an alias (as long as it is unique).

    It looks in `commands` folder for commands.

    This is borrowed from `click` examples of alias and complex.
    '''
    def list_commands(self, ctx):
        rv = []
        cmd_folder = abspath(join(dirname(__file__), 'commands'))
        for filename in listdir(cmd_folder):
            if filename.endswith('.py') and filename != '__init__.py':
                rv.append(splitext(filename)[0])
        return rv

    def get_command(self, ctx, cmd_name):
        # step one: if it matches the full name
        rv = click.Group.get_command(self, ctx, cmd_name)
        if rv is not None:
            return rv
        # step two: allow auto abbr. and lower/upper cases
        matches = [x for x in self.list_commands(ctx)
                   if x.lower().startswith(cmd_name.lower())]
        if not matches:
            return None
        elif len(matches) == 1:
            try:
                mod = __import__('micronota.commands.' + matches[0],
                                 None, None, ['cli'])
            except ImportError as err:
                print(err)
                return
            return mod.cli
        else:
            ctx.fail('Too many matches: %s' % ', '.join(sorted(matches)))


@click.group(cls=AliasedGroup, context_settings=_CONTEXT_SETTINGS)
@click.option('-v', '--verbose', count=True,
              help=("Verbosity. You can use multiple v's, "
                    "e.g. -vv, to gradually increase verbosity."))
@click.option('-d', '--debug', is_flag=True,
              help='Debug mode. Save intermediate files. ')
@click.option('--info', is_flag=True,
              help=('Print micronota configuration details, database info, '
                    'and external dependencies.'))
@click.version_option()   # add --version option
@click.pass_context
def cli(ctx, debug, verbose, info):
    pass
