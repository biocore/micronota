# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from pkg_resources import resource_filename
from os.path import basename, abspath
from logging import getLogger
from subprocess import run

from snakemake import snakemake


def annotate(in_fp, in_fmt, out_dir, out_fmt, gcode,
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
    logger = getLogger(__name__)
    logger.info('Running annotation pipeline')

    snakefile = resource_filename(__name__, 'rules/Snakefile')
    configfile = resource_filename(__name__, 'support_files/config.yaml')

    out = basename(in_fp) + '.gff'
    snakemake(
        snakefile,
        targets=[out],
        cores=cpus,
        # set work dir to output dir so simultaneous runs
        # doesn't interfere with each other.
        workdir=out_dir,
        printshellcmds=True,
        dryrun=dry_run,
        forcetargets=force,
        config={'seq': abspath(in_fp)},
        configfile=configfile,
        keep_target_files=True)


def parse_annotation(out_fp, in_fp, feature_res, annotate_res):
    '''Parse all the annotations and write to disk.'''
