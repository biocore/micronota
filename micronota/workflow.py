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
from skbio.io import read, write


def annotate(in_fp, out_dir, gcode,
             cpus, force, dry_run, config):
    '''Annotate the sequences in the input file.

    Parameters
    ----------
    in_fp : str
        Input seq file name
    out_dir : str
        Output file directory.
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

    success = snakemake(
        snakefile,
        cores=cpus,
        # set work dir to output dir so simultaneous runs
        # doesn't interfere with each other.
        workdir=out_dir,
        printshellcmds=True,
        dryrun=dry_run,
        forcetargets=force,
        config={'seq': in_fp, 'genetic_code': gcode},
        configfile=configfile,
        keep_target_files=True)

    return success


def validate_seq(in_fp, in_fmt, min_len, out_fp):
    '''Validate input seq file.

    1. filter out short seq;
    2. validate seq IDs (no duplicates)
    3. convert to fasta format

    Parameters
    ----------
    in_fp : str
        input seq file path
    in_fmt : str
        the format of seq file
    min_len : int
        cutoff of seq len to filter
    out_fp : str
        output seq file path
    '''
    ids = set()
    with open(out_fp, 'w') as out:
        for seq in read(in_fp, format=in_fmt):
            if len(seq) < min_len:
                continue
            ident = seq.metadata['id']
            if ident in ids:
                raise ValueError(
                    'Duplicate seq IDs in your input file: {}'.format(ident))
            else:
                ids.add(ident)
            write(seq, format='fasta', into=out)


def integrate(out_fp, seq_fp, *kwargs):
    '''integrate all the annotations and write to disk.

    seq_fp : str
        input seq file path.
    out_fp : str
        annotated output file path.
    kwargs : dict
        keys are the formats and valuse are the files to parse.
    '''

