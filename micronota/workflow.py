# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import os
from os.path import splitext, basename, join
from logging import getLogger
from subprocess import run

from .util import _overwrite, convert
from . import _RULES_PATH, _CONFIG_PATH


def annotate(in_fp, in_fmt, out_dir, out_fmt,
             cpus, force, dry_run, config):
    '''Annotate the sequences in the input file.

    Parameters
    ----------
    in_fp : str
        Input file path
    in_fmt : str
        Input file format.
    out_dir : str
        Output file directory.
    out_fmt : str
        Output file format.
    cpus : int
        Number of cpus to use.
    force : boolean
        Force to overwrite.
    config : ``micronota.config.Configuration``
        Container for configuration options.
    '''
    snakefile = join(_RULES_PATH, 'Snakefile')
    configfile = join(_CONFIG_PATH, 'config.yaml')
    cmd = ['snakemake', '-p']
    if force:
        cmd.append('-F')
    if dry_run:
        cmd.append('-n')
    out = join(out_dir, basename(in_fp) + '.gff')
    # this needs to be put at the end
    cmd.append('--snakefile {} {} --configfile {} --config output_dir={} threads={} seq="{}"'.format(
        snakefile, out, configfile, out_dir, cpus, in_fp))
    print(' '.join(cmd))
    proc = run(' '.join(cmd), shell=True)


def parse_annotation(out_fp, in_fp, feature_res, annotate_res):
    '''Parse all the annotations and write to disk.'''
