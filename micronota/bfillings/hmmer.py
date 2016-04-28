# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger
from os.path import join, splitext, basename

from dumpling import Dumpling, OptionParam, ArgmntParam, Parameters

from ._base import _scan_params


def run_hmmpress(hmm, force=False):
    '''Compress the HMM database.

    Parameters
    ----------
    hmm : str
        The file path to HMM database.
    force : boolean
        Whether to overwrite.'''
    params = [
        OptionParam('-f', name='force', help='force overwrite'),
        ArgmntParam(name='hmm', help='hmm file to press')]
    hmmpress = Dumpling('hmmpress', params=Parameters(*params))
    hmmpress.update(hmm=hmm, force=force)
    hmmpress()
    return hmmpress


def run(db, query, out_dir, **kwargs):
    '''Scan a fasta file against a covariance model database.

    Parameters
    ----------
    db : str
        The file path to HMM database.
    query : str
        Input fasta file.
    out_dir : str
        dir to store output file path of target hits table.
    kwargs : dict
        Other command line parameters for hmmscan. key is the option
        (e.g. "-T") and value is the value for the option (e.g. "50").
        If the option is a flag, set the value to None.

    Returns
    -------
    `Dumpling`
    '''
    logger = getLogger(__name__)
    prefix = splitext(basename(query))[0]
    out = join(out_dir, '{}.tblout'.format(prefix))
    hmmscan = Dumpling('hmmscan', params=Parameters(*_scan_params))
    hmmscan.update(query=query, db=db, out=out, **kwargs)
    logger.info('Running {}'.format(hmmscan.command))
    hmmscan()
    return hmmscan
