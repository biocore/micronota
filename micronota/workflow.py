# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, basename, join, exists
from os import makedirs, stat
from importlib import import_module
from tempfile import NamedTemporaryFile
from logging import getLogger

from skbio.metadata import IntervalMetadata
from skbio import read, Sequence
import pandas as pd

from . import bfillings
from .util import _overwrite
from . bfillings.diamond import DiamondCache


def annotate(in_fp, in_fmt, out_dir, out_fmt,
             cpus, kingdom, force, config, cache=False):
    '''Annotate the sequences in the input file.

    Parameters
    ----------
    in_fp : file_handle
        Input file handler object.
    in_fmt : str
        Input file format.
    out_dir : str
        Output file directory.
    out_fmt : str
        Output file format.
    kingdom : int
        Kingdom index corresponding to database (i.e. virus, bacteria ...)
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

    # declare DiamondCache
    if cache:
        cache = DiamondCache()
    else:
        cache = None

    with open(out_fp, 'w') as out:
        # cache
        # submit slurm jobs
        for seq in read(in_fp, format=in_fmt):
            # dir for useful intermediate files for the current input seq
            # replace non alnum char with "_"
            seq_fn = ''.join(x if x.isalnum() else '_'
                             for x in seq.metadata['id'])
            seq_dir = join(out_dir, seq_fn)
            # identify all features specified
            im = identify_all_features(seq, seq_dir, config)
            # pass in and retrieve DiamondCache
            im, cache = annotate_all_cds(im, seq_dir, kingdom, config, cache=cache)

            seq.interval_metadata.concat(IntervalMetadata(im), inplace=True)
            seq.write(out, format=out_fmt)


def identify_all_features(seq, out_dir, config):
    '''Identify all the features for the input sequence.

    It runs through all the tasks specified in sequential order.

    Parameters
    ----------
    seq : skbio.Sequence
        Input sequence object.
    out_dir : str
        Output directory.
    config : ``micronota.config.Configuration``
        Container for configuration options.

    Returns
    -------
    dict :
        Dictionary of skbio.metadata.Feature objects.
    '''
    logger = getLogger(__name__)
    logger.info('Running feature identification.')
    im = dict()
    with NamedTemporaryFile('w+') as f:
        seq.write(f.name, format='fasta')
        for tool in config.features:
            db = config.features[tool]
            if db is not None:
                db = config.db[db]
            seq_dir = join(out_dir, tool)
            submodule = import_module('.%s' % tool, bfillings.__name__)
            cls = getattr(submodule, 'FeaturePred')
            obj = cls(db, seq_dir)
            if tool in config.param:
                params = config.param[tool]
            else:
                params = None
            im.update(next(obj(f.name, params=params)))
    return im


def annotate_all_cds(im, out_dir, kingdom, config, cpus=1, cache=None):
    '''Annotate coding domain sequences (CDS).

    Parameters
    ----------
    seq : skbio.Sequence
        Input sequence object.
    out_dir : str
        Output directory.
    config : ``micronota.config.Configuration``
        Container for configuration options.
    kingdom : str
        Kingdom (i.e. virus, bacteria ...) of the input sequence. It will
        be used to prioritize databases to search.
    cpus : int
        Number of CPUs to use.

    Returns
    -------
    im : skbio.metadata.IntervalMetadata
        Interval metadata object
    '''
    logger = getLogger(__name__)
    logger.info('Running CDS functional annotation.')
    id_key = 'id'
    res = pd.DataFrame()
    for tool in config.cds:
        d = join(out_dir, tool)
        makedirs(d, exist_ok=True)
        pro_fp = join(d, '%s.fa' % tool)

        # write the protein seq into a file
        _write_cds(
            pro_fp, im, id_key,
            lambda x: x['type_'] == 'CDS' and x[id_key] not in res.index)
        if stat(pro_fp).st_size == 0:
            break
        db = config.cds[tool]
        if db in ['uniref100', 'uniref90', 'uniref50']:
            db_dir = config.db[db]
            db_fp = [join(db_dir, i) for i in _get_uniref_db(kingdom)]
            # in case the db file is empty
            db_fp = [i for i in db_fp if exists('%s.dmnd' % i)]
        elif db == 'tigrfam':
            pass
        else:
            raise ValueError('Database %s is not available.' % db)

        submodule = import_module('.%s' % tool, bfillings.__name__)
        cls = getattr(submodule, 'FeatureAnnt')
        if tool == 'diamond':
            obj = cls(dat=db_fp, out_dir=d, cache=cache)
        else:
            obj = cls(dat=db_fp, out_dir=d)
        if tool in config.param:
            params = config.param[tool]
        else:
            params = None
        res_ = obj(pro_fp, cpus=cpus, params=params)
        res = res.append(res_)
    return _update(im, id_key, res), obj.cache


def _update(im, id_key, res):
    '''
    Parameters
    ----------
    im : dict passable to IntervalMetadata
    '''
    features = list(im)
    for feature in features:
        id = feature[id_key]
        if id in res.index:
            new_feature = feature.update(db_xref=res.loc[id, 'sseqid'])
            im[new_feature] = im.pop(feature)
    return im


def _write_cds(fp, im, id_key, select=lambda x: x['type_'] == 'CDS'):
    '''Return a fasta file of all the proteins of a sequence.

    Parameters
    ----------
    fp : str
        The output file path
    im : iterable of Feature
    id_key : str
        key in ``Feature`` to get its value as seq ID
    select : callable
        what CDS to write down.
    '''
    with open(fp, 'w') as f:
        for feature in im:
            if select(feature):
                pro = Sequence(
                    feature['translation'], {'id': feature[id_key]})
                pro.write(f, format='fasta')


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
