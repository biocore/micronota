r'''
Prodigal
========
This wraps Prodigal for CDS prediction.


Notes
-----
It does not handle genes with introns or deal with frame shift.
It is not well tested on viral gene prediction.

'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from os import makedirs
from logging import getLogger
import re

from skbio import read
from skbio.metadata import Feature
from skbio.io.format.genbank import _parse_features
from dumpling import (
    check_choice,
    Dumpling, OptionParam, Parameters)

from ..parsers.embl import _parse_records


params = [
    OptionParam('-i', 'input', help='FASTA/Genbank input file (default reads from stdin).'),
    OptionParam('-a', help='File to store protein translations.'),
    OptionParam('-d', help='File to store nuc sequence of predicted gene.'),
    OptionParam('-s', help='Write all potential genes (with scores) to the selected file.'),
    OptionParam('-t', help=('Write a training file (if none exists); '
                            'otherwise, read and use the specified training file.')),
    OptionParam('-o', 'output', help='output file (default writes to stdout).'),
    OptionParam('-f', 'fmt', action=check_choice(('gbk', 'gff', 'sco')),
                help='output format (gbk, gff, or sco).  Default is gbk.'),
    OptionParam('-g', help='translation table to use (default 11).'),
    OptionParam('-m', help='Treat runs of N as masked sequence; do not build genes across them.'),
    OptionParam('-p', 'mode', action=check_choice(('single', 'meta')),
                help='Select procedure (single or meta).  Default is single.'),
    OptionParam('-c', help='Closed ends.  Do not allow genes to run off edges.'),
    OptionParam('-v', help='Print version number and exit.'),
    OptionParam('-q', help='Run quietly (suppress normal stderr output).'),
    OptionParam('-h', help='Print help menu and exit.'),
    OptionParam('-n', help='Bypass Shine-Dalgarno trainer and force a full motif scan.')]


def run(out_dir, **kwargs):
    '''Run prodigal for gene prediction.

    Notes
    -----
    It will create 3 output files with prefix of "prodigal" in "out_dir" folder:
      1. the annotation file (GFF3 file by default)
      2. the nucleotide sequences for each predicted gene
         with file suffix of .fna.
      3. the protein sequence translated from each gene
         with file suffix of .faa.

    Parameters
    ----------
    out_dir : str
        output dir
    kwargs : dict
        keyword arguments for Prodigal.

    Returns
    -------
    `Dumpling`
    '''
    logger = getLogger(__name__)

    makedirs(out_dir, exist_ok=True)

    prodigal = Dumpling('prodigal', params=Parameters(*params),
                        version='v2.6.3', url='https://github.com/hyattpd/Prodigal')
    # set default output to gff and run mode to draft genome
    prodigal.update(fmt='gff', mode='single')
    # update with kwargs
    prodigal.update(**kwargs)

    suffices = {
        '-a': 'faa',
        # output file of nucleotide sequences of genes
        '-d': 'fna',
        '-o': prodigal.params['fmt'].value}

    for flag in ['-a', '-d', '-o']:
        p = prodigal.params[flag]
        if p.is_off():
            p.on(join(out_dir, 'prodigal.{}'.format(suffices[flag])))

    logger.info('Running CDS prediction {}'.format(prodigal.command))
    prodigal(stdout=join(out_dir, 'prodigal.log'))
    return prodigal


def parse_result(self, res, which='-a'):
    '''Parse gene prediction result from ``Prodigal``.

    It is parsed into a dict consumable by
    ``skbio.metadata.IntervalMetadata``.

    Parameters
    ----------
    res : burrito.util.CommandLineAppResult
    which : which output to parse
    Returns
    -------
    ``skbio.metadata.IntervalMetadata``
    '''
    # make sure to move to the beginning of the file.
    if which == '-a':
        return _parse_faa(res[which])
    elif which == '-o':
        return _parse_records(res[which], _parse_single_record)


def _parse_single_record(chunks):
    '''Parse single record of Prodigal GenBank output.

    Parameters
    ----------
    chunks : list of str
        a list of lines of the record to parse.

    Yields
    ------
    dict passable to ``skbio.metadata.IntervalMetadata``
    '''
    # get the head line
    head = chunks[0]
    _, description = head.split(None, 1)
    pattern = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
    desc = {}
    for i in pattern.findall(description):
        k, v = i.split('=', 1)
        desc[k] = v
    return _parse_features(chunks[1:], int(desc['seqlen']))


def _parse_faa(faa):
    '''Parse the faa output of Prodigal.

    Notice that Prodigal do not predict genes with introns.

    Yields
    ------
    dict passable to ``skbio.metadata.IntervalMetadata``.
    '''
    pattern = (r'# +([0-9]+)'    # start
               ' +# +([0-9]+)'   # end
               ' +# +(-?1)'      # strand
               ' +# +ID=([0-9]+_[0-9]+);'
               'partial=([01]{2});'
               '(.*)')
    im = dict()
    i = 1
    for seq in read(faa, format='fasta'):
        desc = seq.metadata['description']
        matches = re.match(pattern, desc)
        start, end, strand, id, partial, misc = matches.groups()
        # ordinal number of the parent seq
        ordinal = int(id.split('_', 1)[0])
        if ordinal > i:
            yield im
            # reset
            i += 1
            im = dict()
        interval = []
        feature = dict()
        feature['translation'] = str(seq)
        feature['type_'] = 'CDS'

        # don't forget to convert 0-based
        interval = [(int(start)-1, int(end))]
        feature['note'] = '"%s"' % misc
        feature['id'] = id
        if partial[0] == '0':
            feature['left_partial_'] = False
        else:
            feature['left_partial_'] = True
            start = '<%s' % start
        if partial[1] == '0':
            feature['right_partial_'] = False
        else:
            feature['right_partial_'] = True
            end = '>%s' % end
        location = '{s}..{e}'.format(s=start, e=end)
        if strand == '-1':
            feature['rc_'] = True
            location = 'complement(%s)' % location
        elif strand == '1':
            feature['rc_'] = False
        else:
            raise ValueError('Inappropriate value for strand: %s' % strand)
        feature['location'] = location
        im[Feature(**feature)] = interval

    if im:
        # don't forget to return the last one if it is not empty.
        yield im
