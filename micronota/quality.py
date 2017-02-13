# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict
from subprocess import run
from os.path import join, dirname
from tempfile import NamedTemporaryFile
import logging

from pkg_resources import resource_filename
from skbio.io import write

logger = logging.getLogger(__name__)


def compute_gene_score(seq_fp):
    '''Compute the quality score based on essential genes found on sequences.

    Parameters
    ----------
    seq_fp : str
        file path of proteins identified from the DNA sequences

    Returns
    -------
    float
        the score computed from sequence stats.
    '''
    hmm_fp = resource_filename(__package__, 'data/pfam/102.hmm')
    with NamedTemporaryFile() as out:
        proc = run(["hmmscan", "--tblout", out.name, "-o", '/dev/null', hmm_fp, seq_fp])
        p102 = set()
        with open(out.name) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                items = line.split()
                pid = items[1].split('.')[0]
                p102.add(pid)
    score = max(1 - (102 - len(p102)) * 0.01, 0.1)
    return score


def compute_seq_score(seqs, contigs=True):
    '''Compute the quality score solely based on sequence stats.

    Parameters
    ----------
    seqs : Iterable of ``skbio.Sequence`` or its child
        input sequences
    contigs : bool, optional
        whether the input seqs are contigs or chromosomes

    Returns
    -------
    float
        the score computed from sequence stats.
    '''
    ncontigs, ngood, nbad, n10N = 0, 0, 0, 0
    for i, seq in enumerate(seqs, 1):
        pattern = '(%s+)' % ('N' * 10)
        n10N += len([_ for _ in seq.find_with_regex(pattern)])

        freq = seq.frequencies()
        ngood += sum(freq.pop(k, 0) for k in ('A', 'T', 'G', 'C', 'N'))

        nbad += sum(v for _, v in freq.items())

    if contigs is False:
        ncontigs = 1
    else:
        ncontigs = i

    score = ngood / (ngood + nbad + 10000 * (ncontigs - 1) + 10000 * n10N)

    return score


def compute_rrna_score(imds):
    '''Compute the quality score based on rRNA found on sequences.

    Parameters
    ----------
    imds : Iterable of ``skbio.metadata.IntervalMetadata``
    contigs : bool, optional
        whether the input seqs are contigs or chromosomes

    Returns
    -------
    float
        the score computed from sequence stats.
    '''
    rrnas = {'5s_rRNA': (100, 120),
             '16s_rRNA': (1450, 1700),
             '23s_rRNA': (2900, 3500)}
    scores = defaultdict(list)
    for imd in imds:
        for rrna in imd.query(metadata={'type': 'rRNA'}):
            start, end = rrna.bounds[0]
            length = end - start
            t = rrna.metadata['product']
            lower, upper = rrnas[t]
            if lower <= length <= upper:
                scores[t].append(0.3)
            elif lower * 0.5 <= length:
                scores[t].append(0.2)
            else:
                scores[t].append(0.1)
    score = 0.1 + sum([max(v) for v in scores.values()])
    return score


def compute_trna_score(imds):
    '''Compute the quality score based on tRNA found on sequences.

    Parameters
    ----------
    imds : Iterable of ``skbio.metadata.IntervalMetadata``
    contigs : bool, optional
        whether the input seqs are contigs or chromosomes

    Returns
    -------
    float
        the score computed from sequence stats.
    '''
    aa = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
          'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'}
    for imd in imds:
        for trna in imd.query(metadata={'type': 'tRNA'}):
            a = trna.metadata['product'].split('-')[-1]
            aa.discard(a)

    score = max(0.1, 1 - len(aa) * 0.1)
    return score
