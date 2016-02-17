# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click
from ..util import get_config_info


@click.command()
@click.pass_context
def cli(ctx):
    '''Print micronota configuration details.

    It includes system info, database info, and external dependencies, etc.
    '''
    info = get_config_info(ctx.parent.config)
    for key in sorted(info):
        click.echo()
        click.echo(key)
        click.echo('=' * len(key))
        info_of_key = info[key]
        try:
            max_len = max([len(i) for i in info_of_key])
        except ValueError:
            continue
        for k in sorted(info_of_key):
            v = ' '.join(info_of_key[k].split('\n'))
            line = "{key:<{width}}: {value}".format(
                width=max_len, key=k, value=v)
            click.echo(line)
