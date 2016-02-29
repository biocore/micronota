# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, basename, join
from os import remove
from importlib import import_module
from tempfile import NamedTemporaryFile

from skbio.metadata import IntervalMetadata
from skbio import read, write, Sequence

from . import bfillings
from .util import _DB_PATH


def annotate(in_fp, in_fmt, out_dir, out_fmt, kingdom, cpus, config):
    '''Annotate the sequences in the input file.'''
    fn = splitext(basename(in_fp))[0]
    # store annotated seq file.
    out_fp = join(out_dir, '%s.gbk' % fn)
    out = open(out_fp, 'w')
    for seq in read(in_fp, format=in_fmt):
        # dir for useful intermediate files for the current input seq
        seq_tmp_dir = join(out_dir, seq.metadata['id'])
        im = {}
        for sec in config:
            if sec == 'GENERAL':
                continue
            if sec == 'FEATURE':
                im.update(identify_all_features(seq, seq_tmp_dir, config[sec]))
            if sec == 'CDS':
                annotate_all_cds(im, out_dir, kingdom, config[sec])

        seq.interval_metadata.concat(im, inplace=True)
        seq.write(out, format=out_fmt)
    out.close()


def identify_all_features(seq, out_dir, tasks,
                          id_func='identify_features',
                          parse_func='parse_output'):
    '''Identify all the features for the sequence in the input file.

    It runs thru all the tasks specified in sequential order.

    Parameters
    ----------
    tasks : OrderDict-like
    Returns
    -------
    '''
    im = dict()
    with NamedTemporaryFile('w+') as f:
        seq.write(f.name, format='fasta')

        for tool in tasks:
            try:
                value = tasks.getboolean(tool)
            except ValueError:
                value = tasks[tool]
            if not value:
                continue
            try:
                submodule = import_module('.%s' % tool, bfillings.__name__)
            except ImportError:
                raise ImportError('Annotation with %s seems not implemented.' % tool)
            id_f = getattr(submodule, id_func)
            parse_f = getattr(submodule, parse_func)
            if value is True:
                # value is just an boolean indicator
                res = id_f(f.name, out_dir)
            else:
                # value is database name
                res = id_f(f.name, out_dir, value)
            im.update(next(parse_f(res)))
    return im


def annotate_all_cds(im, out_dir, kingdom, tasks,
                     search_func='search_protein_homologs',
                     parse_func='parse_output'):
    '''Annotate CDS.'''
    uniref_dbs = ['uniref100_Swiss-Prot_Bacteria.fasta',
                  'uniref100_Swiss-Prot_Archaea.fasta',
                  'uniref100_Swiss-Prot_Viruses.fasta',
                  'uniref100_Swiss-Prot_Eukaryota.fasta',
                  'uniref100_Swiss-Prot_other.fasta',
                  'uniref100_TrEMBL_Bacteria.fasta',
                  'uniref100_TrEMBL_Archaea.fasta',
                  'uniref100_TrEMBL_Viruses.fasta',
                  'uniref100_TrEMBL_Eukaryota.fasta',
                  'uniref100_TrEMBL_other.fasta',
                  'uniref100__other.fasta']

    uniref_dbs = [join(_DB_PATH, i) for i in uniref_dbs]
    order = {'bacteria': range(11),
             'archaea': [1, 0, 2, 3, 4, 6, 5, 7, 8, 9, 10],
             'viruses': [2, 0, 1, 3, 4, 7, 5, 6, 8, 9, 10]}

    id_old_features = {}
    for feature in im:
        if feature['type_'] == 'CDS':
            id_old_features[feature['id']] = feature

    for tool in tasks:
        value = tasks[tool]
        submodule = import_module('.%s' % tool, bfillings.__name__)

        search_f = getattr(submodule, search_func)
        parse_f = getattr(submodule, parse_func)

        if value == 'uniref':
            dbs = [uniref_dbs[i] for i in order[kingdom.lower()]]

            def _uniref(id_old_features, dbs, id_new_features=None):
                if id_new_features is None:
                    id_new_features = dict()
                if not id_old_features or not dbs:
                    return id_new_features

                for i in id_old_features:
                    tmp = NamedTemporaryFile('w+', delete=False)
                    pro = Sequence(id_old_features[i]['translation'], {'id': i})
                    pro.write(tmp.name, format='fasta')
                db = dbs[0]
                res = search_f(tmp, db, out_dir)
                hits = parse_f(res)
                for idx, row in hits.iterrows():
                    feature = id_old_features.pop(idx)
                    id_new_features[idx] = feature.update(db_xref=row['sseqid'])

                    tmp.close()
                    remove(tmp.name)
                return _uniref(id_old_features, dbs[1:], id_new_features)
            id_new_features = _uniref(id_old_features, dbs)
            for id in id_new_features:
                im[id_new_features[id]] = im.pop(id_old_features[id])

        return im
