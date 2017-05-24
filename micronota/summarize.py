# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


def summarize(obj, types=('length', 'nuc_freq', 'CDS', 'ncRNA', 'rRNA', 'tRNA', 'tandem_repeat', 'terminator', 'CRISPR')):
    '''Summarize the sequence and its annotation.

    Parameters
    ----------
    seq : ``Sequence`` or ``IntervalMetadata`` object

    Returns
    -------
    tuple
        summary stat
    '''


def summarize_iter(objs, types=('length', 'nuc_freq', 'CDS', 'ncRNA', 'rRNA', 'tRNA', 'tandem_repeat', 'terminator', 'CRISPR')):
    '''Summarize the sequences and their annotations in a genome or metagenome sample.

    Parameters
    ----------
    seq : iterable of ``Sequence`` or ``IntervalMetadata`` objects

    Yields
    ------
    tuple
        summary stat
    '''


def compute_condon_usage(cds, genetic_code=11):
    '''Compute the condon usage of the input sequence

    Parameters
    ----------
    cds : ``DNA`` or ``RNA``
        CDS sequence
    genetic_code : int
        genetic code/translation table

    Returns
    -------
    dict
    '''
