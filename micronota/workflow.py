# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, basename, join
from importlib import import_module

from skbio.sequence import IntervalMetadata
from skbio import read, DNA

from .util import _tmp_file
from . import bfillings


def annotate(in_fp, in_fmt, out_dir, out_fmt, cpus, config):
    '''Annotate the sequences in the input file.'''
    fn = splitext(basename(in_fp))[0]
    # store annotated seq file.
    out_fp = join(out_dir, '%s.gbk' % fn)
    # append the file for multiple sequences
    out = open(out_fp, 'a')
    sec = 'WORKFLOW'
    workflows = config[sec]['order'].split(' > ')
    for seq in read(in_fp, format=in_fmt):
        # interval_metadata
        im = {}
        with _tmp_file() as f:
            _, fp = f
            seq.write(fp, format='fasta')
            seq_id = seq.metadata['id']
            for k in workflows:
                tasks = config[sec][k].split(' > ')
                if k == 'features':
                    im_ = identify_all_features(fp, seq_id, tasks, config)
                if k == 'CDS':
                    annotate_all_cds()
        seq.interval_metadata = IntervalMetadata(im)
        seq.write(out, format=out_fmt)
    out.close()


def identify_all_features(
        fp, out_dir, tasks, config, func_name='identify_features'):
    im = {}
    for task in tasks:
        items = task.split(':')
        tool = items[0]
        submodule = import_module('.%s' % tool, bfillings.__name__)
        f = getattr(submodule, func_name)

        params = None
        if tool in config:
            params = config._sections[tool]
        print(out_dir)
        im.update(f(fp, out_dir, params=params))
    return im


def annotate_all_cds():
    '''Annotate CDS.'''
