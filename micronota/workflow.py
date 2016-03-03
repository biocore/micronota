# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, basename, join, exists
from os import remove, makedirs
from importlib import import_module
from tempfile import NamedTemporaryFile
from logging import getLogger

from skbio.metadata import IntervalMetadata
from skbio import read, Sequence

from . import bfillings
from .util import _overwrite


def annotate(in_fp, in_fmt, out_dir, out_fmt,
             cpus, kingdom, dry, force, config):
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
    dry : boolean
        Run the real computation if False
    force : boolean
        Force to overwrite.
    config : ``micronota.config.Configuration``
        Container for configuration options.
    '''
    if not dry:
        _overwrite(out_dir, overwrite=force)
        makedirs(out_dir, exist_ok=force)
    prefix = splitext(basename(in_fp))[0]
    fn = '{p}.{f}'.format(p=prefix, f=out_fmt)
    out_fp = join(out_dir, fn)
    with open(out_fp, 'w') as out:
        for seq in read(in_fp, format=in_fmt):
            # dir for useful intermediate files for the current input seq
            # replace non alnum char with "_"
            seq_fn = ''.join(x if x.isalnum() else '_'
                             for x in seq.metadata['id'])
            seq_dir = join(out_dir, seq_fn)
            # identify all features specified
            im = identify_all_features(seq, seq_dir, config)
            seq.interval_metadata = IntervalMetadata(im)
            im = annotate_all_cds(seq, out_dir, kingdom, config)
            seq.interval_metadata.concat(im, inplace=True)
            seq.write(out, format=out_fmt)


def identify_all_features(seq, out_dir, config,
                          id_func='identify_features',
                          parse_func='parse_output'):
    '''Identify all the features for the sequence in the input file.

    It runs through all the tasks specified in sequential order.

    Parameters
    ----------
    seq : skbio.Sequence
        Input sequence object.
    out_dir : str
        Output directory.
    config : ``micronota.config.Configuration``
        Container for configuration options.
    id_func: str
        Function to id features.
    parse_func : str
        Parser function found in bfillings.

    Returns
    -------
    dict :
        Dictionary of skbio.metadata.Feature objects.
    '''
    logger = getLogger(__name__)
    logger.info('Running feature identification')
    im = dict()
    with NamedTemporaryFile('w+') as f:
        seq.write(f.name, format='fasta')
        for tool in config.features:
            db = config.features[tool]
            submodule = import_module('.%s' % tool, bfillings.__name__)
            id_f = getattr(submodule, id_func)
            parse_f = getattr(submodule, parse_func)
            if db is None:
                res = id_f(f.name, out_dir)
            else:
                db = join(config.db_dir, db)
                res = id_f(f.name, out_dir, db)
            im.update(next(parse_f(res)))
    return im


def annotate_all_cds(seq, out_dir, kingdom, config,
                     search_func='search_protein_homologs',
                     parse_func='parse_output'):
    '''Annotate coding domain sequences (CDS).

    Parameters
    ----------
    seq : skbio.Sequence
        Input sequence object.
    out_dir : str
        Output directory.
    config : ``micronota.config.Configuration``
        Container for configuration options.
    tasks : OrderDict-like
        Ordered dictionary of tools to run and their parameters.
    parse_func : str
        Parser function found in bfillings.
    kingdom : int
        Kingdom index corresponding to database (i.e. virus, bacteria ...)
    cpus : int
        Number of cpus to use.

    Returns
    -------
    im : skbio.metadata.IntervalMetadata
        Interval metadata object
    '''
    logger = getLogger(__name__)
    logger.info('Running CDS functional annotation')
    uniref_dbs = ['uniref100_Swiss-Prot_Bacteria',
                  'uniref100_Swiss-Prot_Archaea',
                  'uniref100_Swiss-Prot_Viruses',
                  'uniref100_Swiss-Prot_Eukaryota',
                  'uniref100_Swiss-Prot_other',
                  'uniref100_TrEMBL_Bacteria',
                  'uniref100_TrEMBL_Archaea',
                  'uniref100_TrEMBL_Viruses',
                  'uniref100_TrEMBL_Eukaryota',
                  'uniref100_TrEMBL_other',
                  'uniref100__other']

    uniref_dbs = [join(config.db_dir, i) for i in uniref_dbs]
    order = {'bacteria': range(11),
             'archaea': [1, 0, 2, 3, 4, 6, 5, 7, 8, 9, 10],
             'viruses': [2, 0, 1, 3, 4, 7, 5, 6, 8, 9, 10]}
    im = seq.interval_metadata
    id_old_cds = dict()
    id_new_cds = dict()
    for feature in im.query(type_='CDS'):
        id_old_cds[feature['id']] = feature

    for tool in config.cds:
        db = config.cds[tool]
        submodule = import_module('.%s' % tool, bfillings.__name__)
        search_f = getattr(submodule, search_func)
        parse_f = getattr(submodule, parse_func)

        if db == 'uniref':
            dbs = [uniref_dbs[i] for i in order[kingdom.lower()]]
            for db_fp in dbs:
                if not id_old_cds:
                    break
                if not exists('%s.dmnd' % db_fp):
                    continue
                tmp = NamedTemporaryFile('w+', delete=False)
                tmp.close()

                out = open(tmp.name, 'w')
                for i in id_old_cds:
                    pro = Sequence(id_old_cds[i]['translation'], {'id': i})
                    pro.write(out, format='fasta')

                res = search_f(tmp.name, db_fp, out_dir)
                hits = parse_f(res)
                for idx, row in hits.iterrows():
                    old_cds = id_old_cds.pop(idx)
                    if 'db_xref' in old_cds:
                        db_xref = ','.join(old_cds['db_xref'], row['sseqid'])
                    else:
                        db_xref = row['sseqid']
                    new_cds = old_cds.update(db_xref=db_xref)
                    im.update(old_cds, new_cds)

                out.close()
                remove(tmp.name)

            for _id in id_new_cds:
                im[id_new_cds[_id]] = im.pop(id_old_cds[_id])
    return im
