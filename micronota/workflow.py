# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from pkg_resources import resource_filename
from os.path import join, exists
from logging import getLogger
from importlib import import_module

from snakemake import snakemake
from skbio.io import read, write
import pandas as pd

from . import parsers
from .parsers.cds import _add_protein_annotation, _parse_seq_id


logger = getLogger(__name__)


def annotate(in_fp, out_dir, gcode, cpus, force, dry_run, config):
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
    dry_run : boolean
    config : config file for snakemake
    '''
    logger.info('Running annotation pipeline')
    snakefile = resource_filename(__name__, 'rules/Snakefile')
    if config is None:
        config = resource_filename(__name__, 'rules/config.yaml')

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
        configfile=config,
        keep_target_files=True,
        keep_logger=True)

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
    if exists(out_fp):
        return
    logger.info('Filtering and validating input sequences')
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


def integrate(out_dir, seq_fn, out_fmt='genbank'):
    '''integrate all the annotations and write to disk.

    seq_fn : str
        input seq file path.
    out_dir : str
        annotated output file path.
    out_fmt : str
        output format
    '''
    logger.info('Integrating annotation for output')
    imd = {}
    seqs = []
    for seq in read(join(out_dir, seq_fn), format='fasta'):
        imd[seq.metadata['id']] = seq.interval_metadata
        seqs.append(seq)

    hits = []
    for f in os.listdir(out_dir):
        if f.endswith('.hit'):
            hits.append(join(out_dir, f))
        elif f.endswith('.ok'):
            # parse the output of each feature prediction tool into
            # interval metadata
            tool = f.rsplit('.', 1)[0]
            submodule = import_module('.%s' % tool, parsers.__name__)
            f = getattr(submodule, 'parse')
            f(imd, out_dir)

    # add functional metadata to the protein-coding gene
    hit = pd.concat([pd.read_table(i) for i in hits], axis=0)
    hit_dict = _parse_seq_id(hit)
    _add_protein_annotation(imd, hit_dict)

    # write out the annotation
    with open(join(out_dir, '%s.gbk' % seq_fn), 'w') as out:
        for seq in seqs:
            write(seq, into=out, format='genbank')

