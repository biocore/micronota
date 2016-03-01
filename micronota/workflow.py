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

from skbio.metadata import IntervalMetadata
from skbio import read, Sequence

from . import bfillings
from .util import _DB_PATH


def annotate(in_fp, in_fmt, out_dir, out_fmt, kingdom, cpus, config):
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
    config : ConfigParser
        Container for configuration options.
    '''
    fn = splitext(basename(in_fp))[0]
    # store annotated seq file.
    makedirs(out_dir, exist_ok=True)
    out_fp = join(out_dir, '%s.gbk' % fn)
    out = open(out_fp, 'w')
    for seq in read(in_fp, format=in_fmt):
        # dir for useful intermediate files for the current input seq
        seq_tmp_dir = join(out_dir, seq.metadata['id'])
        for sec in config:
            if sec == 'GENERAL':
                continue
            if sec == 'FEATURE':
                im = identify_all_features(seq, seq_tmp_dir, config[sec])
                seq.interval_metadata = IntervalMetadata(im)
            if sec == 'CDS':
                im = annotate_all_cds(seq, out_dir, kingdom, config[sec])

        seq.interval_metadata.concat(im, inplace=True)
        seq.write(out, format=out_fmt)
    out.close()


def identify_all_features(seq, out_dir, tasks,
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
    tasks : OrderDict-like
        Ordered dictionary of tools to run and their parameters.
    parse_func : str
        Parser function found in bfillings.

    Returns
    -------
    dict :
        Dictionary of skbio.metadata.Feature objects.
    '''
    im = dict()
    with NamedTemporaryFile('w+') as f:
        seq.write(f.name, format='fasta')

        for tool in tasks:
            value = _get_task_value(tasks, tool)
            if not value:
                continue
            try:
                submodule = import_module('.%s' % tool, bfillings.__name__)
            except ImportError:
                raise ImportError('Annotation with %s seems not implemented.' % tool)
            id_f = getattr(submodule, id_func)
            parse_f = getattr(submodule, parse_func)
            if value:
                # value is just an boolean indicator
                res = id_f(f.name, out_dir)
            else:
                # value is database name
                res = id_f(f.name, out_dir, value)
            im.update(next(parse_f(res)))
    return im


def annotate_all_cds(seq, out_dir, kingdom, tasks,
                     search_func='search_protein_homologs',
                     parse_func='parse_output'):
    '''Annotate coding domain sequences (CDS).

    Parameters
    ----------
    seq : skbio.Sequence
        Input sequence object.
    out_dir : str
        Output directory.
    tasks : OrderDict-like
        Ordered dictionary of tools to run and their parameters.
    parse_func : str
        Parser function found in bfillings.
    kingdom : int
        Kingdom index corresponding to database (i.e. virus, bacteria ...)
    cpus : int
        Number of cpus to use.
    config : ConfigParser
        Container for configuration options.

    Returns
    -------
    im : skbio.metadata.IntervalMetadata
        Interval metadata object
    '''
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

    uniref_dbs = [join(_DB_PATH, i) for i in uniref_dbs]
    order = {'bacteria': range(11),
             'archaea': [1, 0, 2, 3, 4, 6, 5, 7, 8, 9, 10],
             'viruses': [2, 0, 1, 3, 4, 7, 5, 6, 8, 9, 10]}
    im = seq.interval_metadata
    id_old_cds = dict()
    id_new_cds = dict()
    for feature in im.query(type_='CDS'):
        id_old_cds[feature['id']] = feature

    for tool in tasks:
        value = _get_task_value(tasks, tool)
        if not value:
            continue

        submodule = import_module('.%s' % tool, bfillings.__name__)

        search_f = getattr(submodule, search_func)
        parse_f = getattr(submodule, parse_func)

        if value == 'uniref':
            dbs = [uniref_dbs[i] for i in order[kingdom.lower()]]
            for db in dbs:
                if not id_old_cds:
                    break
                db = join(_DB_PATH, db)
                if not exists('%s.dmnd' % db):
                    continue
                tmp = NamedTemporaryFile('w+', delete=False)
                tmp.close()

                out = open(tmp.name, 'w')
                for i in id_old_cds:
                    pro = Sequence(id_old_cds[i]['translation'], {'id': i})
                    pro.write(out, format='fasta')

                res = search_f(tmp.name, db, out_dir)
                hits = parse_f(res)
                for idx, row in hits.iterrows():
                    old_cds = id_old_cds.pop(idx)
                    new_cds = feature.update(db_xref=row['sseqid'])
                    im.update(old_cds, new_cds)

                out.close()
                remove(tmp.name)

            for _id in id_new_cds:
                im[id_new_cds[_id]] = im.pop(id_old_cds[_id])
    return im

def _get_task_value(tasks, tool):
    try:
        value = tasks.getboolean(tool)
    except ValueError:
        value = tasks[tool]
    return value
