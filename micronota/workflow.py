# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from collections import defaultdict
from os.path import join, exists, basename
from logging import getLogger
from importlib import import_module

from pkg_resources import resource_filename
from snakemake import snakemake
from skbio.io import read, write
import yaml
import pandas as pd

from . import parsers
from .parsers.cds import _add_cds_metadata, _fetch_cds_metadata


logger = getLogger(__name__)


def annotate(in_fp, in_fmt, min_len, out_dir, out_fmt, gcode,
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
    dry_run : boolean
    config : config file for snakemake
    '''
    logger.info('Running annotation pipeline')

    os.makedirs(out_dir, exist_ok=True)

    seq_fn = basename(in_fp)

    validate_seq(in_fp, in_fmt, min_len, join(out_dir, seq_fn))

    if config is None:
        config = resource_filename(__package__, 'config.yaml')
    with open(config) as fh:
        cfg = yaml.load(fh)

    cfg['general']['seq'] = seq_fn
    cfg['general']['genetic_code'] = gcode

    snakefile = resource_filename(__package__, 'rules/Snakefile')
    success = snakemake(
        snakefile,
        cores=cpus,
        # set work dir to output dir so simultaneous runs
        # doesn't interfere with each other.
        workdir=out_dir,
        printshellcmds=True,
        dryrun=dry_run,
        forcetargets=force,
        config=cfg,
        # configfile=config,
        keep_target_files=True,
        keep_logger=True)

    if success:
        # if snakemake finishes successfully
        integrate(cfg['general'], out_dir, seq_fn, out_fmt)


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


def integrate(cfg, out_dir, seq_fn, out_fmt='genbank'):
    '''integrate all the annotations and write to disk.

    seq_fn : str
        input seq file path.
    out_dir : str
        annotated output file path.
    out_fmt : str
        output format
    '''
    logger.info('Integrating annotation for output')
    seqs = {}
    for seq in read(join(out_dir, seq_fn), format='fasta'):
        seqs[seq.metadata['id']] = seq

    hits = []
    for f in os.listdir(out_dir):
        if f.endswith('.hit'):
            # protein homologous map
            hits.append(join(out_dir, f))
        elif f.endswith('.ok'):
            # parse the output of each feature prediction tool into
            # interval metadata
            tool = f.rsplit('.', 1)[0]
            submodule = import_module('.%s' % tool, parsers.__name__)
            f = getattr(submodule, 'parse')
            for seq_id, imd in f(out_dir):
                seq = seqs[seq_id]
                imd._upper_bound = len(seq)
                seq.interval_metadata.merge(imd)

    # add functional metadata to the protein-coding gene
    # create defaultdict of defaultdict of dict
    cds_metadata = defaultdict(lambda : defaultdict(dict))
    for f in hits:
        for seq_id, idx, md in _fetch_cds_metadata(f, cfg['metadata']):
            cds_metadata[seq_id][idx].update(md)
    if cds_metadata:
        _add_cds_metadata(seqs, cds_metadata)

    # write out the annotation
    with open(join(out_dir, '%s.gbk' % seq_fn), 'w') as out:
        for _, seq in seqs.items():
            seq.interval_metadata.sort()
            write(seq, into=out, format='genbank')
