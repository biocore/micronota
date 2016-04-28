# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, basename, join
from os import makedirs
from importlib import import_module
from logging import getLogger

from . import bfillings
from .util import _overwrite, convert


def annotate(in_fp, in_fmt, out_dir, out_fmt,
             cpus, force, config):
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
    _overwrite(out_dir, overwrite=force)
    makedirs(out_dir, exist_ok=force)

    prefix = splitext(basename(in_fp))[0]
    fn = '{p}.{f}'.format(p=prefix, f=out_fmt)
    out_fp = join(out_dir, fn)
    if in_fmt != 'fasta':
        fna = join(out_dir, '{p}.fna'.format(prefix))
        convert(in_fmt, 'fasta', in_fp, fna)
        in_fp = fna

    feature_res = identify_features(in_fp, out_dir, cpus, force, config)

    pro_fp = None
    if 'prodigal' in feature_res:
        pro_fp = feature_res['prodigal'].params['-a']

    if pro_fp is not None:
        annotate_res = annotate_cds(pro_fp, out_dir, cpus, force, config)

    parse_annotation(out_fp, in_fp, feature_res, annotate_res)


def parse_annotation(out_fp, in_fp, feature_res, annotate_res):
    '''Parse all the annotations and write to disk.'''


def identify_features(seq, out_dir, cpus, force, config):
    '''Identify all the features for the input sequence.

    It runs through all the tasks specified in sequential order.

    Parameters
    ----------
    seq : str
        Input sequence file path.
    out_dir : str
        Output directory.
    cpus : int
        number of cpu cores
    config : ``micronota.config.Configuration``
        Container for configuration options.

    Returns
    -------
    dict
    '''
    logger = getLogger(__name__)
    logger.info('Running feature identification ...')

    res = {}

    for tool in config.features:
        submodule = import_module('.%s' % tool, bfillings.__name__)
        pred = getattr(submodule, 'run')

        params = {}
        if tool in config.param:
            params = config.param[tool]

        d = join(out_dir, tool)
        _overwrite(d, overwrite=force)
        makedirs(d, exist_ok=force)

        db = config.features[tool]
        if db is not None:
            db_path = config.db[db]
            dumpling = pred(seq, d, db=db_path, cpus=cpus, **params)
        else:
            dumpling = pred(seq, d, cpus=cpus, **params)

        res[tool] = dumpling
    return res


def annotate_cds(pro_fp, out_dir, cpus, force, config):
    '''Annotate coding domain sequences (CDS).

    Parameters
    ----------
    pro_fp : str
        input file path of protein seq
    out_dir : str
        Output directory.
    config : ``micronota.config.Configuration``
        Container for configuration options.
    cpus : int
        Number of CPUs to use.

    Returns
    -------
    dict
    '''
    logger = getLogger(__name__)
    logger.info('Running CDS annotation ...')

    for tool in config.cds:
        submodule = import_module('.%s' % tool, bfillings.__name__)
        pred = getattr(submodule, 'run')

        d = join(out_dir, tool)
        _overwrite(d, overwrite=force)
        makedirs(d, exist_ok=True)

        db = config.cds[tool]
        db_path = config.db[db]
        pred(db_path, d)


def _get_uniref_db(kingdom):
    dbs = ['Swiss-Prot_Bacteria',
           'Swiss-Prot_Archaea',
           'Swiss-Prot_Viruses',
           'Swiss-Prot_Eukaryota',
           'Swiss-Prot_other',
           'TrEMBL_Bacteria',
           'TrEMBL_Archaea',
           'TrEMBL_Viruses',
           'TrEMBL_Eukaryota',
           'TrEMBL_other',
           '_other']
    order = {'bacteria': range(11),
             'archaea': [1, 0, 2, 3, 4, 6, 5, 7, 8, 9, 10],
             'viruses': [2, 0, 1, 3, 4, 7, 5, 6, 8, 9, 10]}
    return [dbs[i] for i in order[kingdom.lower()]]
