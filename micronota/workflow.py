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

    cmd = ['snakemake', '-p', '--keep-target-files']
    if force:
        cmd.append('-F')
    if dry_run:
        cmd.append('-n')
    out = basename(in_fp) + '.gff'
    # this needs to be put at the end
    cmd.extend(['--snakefile', snakefile,
                '--cores', str(cpus),
                # set work dir to output dir so simultaneous runs
                # doesn't interfere with each other.
                workdir: config['output_dir']

                '--configfile', configfile,
                out,
                '--config',
                'output_dir=%s' % out_dir, 'seq=%s' % abspath(in_fp)])
    proc = run(cmd)
    return proc


def parse_annotation(out_fp, in_fp, feature_res, annotate_res):
    '''Parse all the annotations and write to disk.'''
