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


cmscan_params = [
    OptionParam('--rfam', help='Set all filter thresholds as if the search space were more than 20 Gb.'),
    OptionParam('--noali', value=True,
                help='Omit the alignment section from the main output. This can greatly reduce the output volume.')]
cmscan_params.extend(_scan_params)


def run_cmpress(cm, force=False):
    '''Compress the CM database.

    Parameters
    ----------
    cm : str
        The file path to CM database.
    force : boolean
        Whether to overwrite.'''
    params = [
        OptionParam('-F', name='force', help='force overwrite'),
        ArgmntParam(name='cm', help='cm file to press')]
    cmpress = Dumpling('cmpress', params=Parameters(*params))
    cmpress.update(cm=cm, force=force)
    cmpress()
    return cmpress


def run(db, query, out_dir, **kwargs):
    '''Scan a fasta file against a covariance model database.

    Parameters
    ----------
    db : str
        The file path to CM database.
    query : str
        Input fasta file.
    out_dir : str
        dir to store output file path of target hits table.
    kwargs : dict
        keyword arguments. command line parameters for cmscan.

    Returns
    -------
    `Dumpling`
    '''
    logger = getLogger(__name__)
    prefix = splitext(basename(query))[0]
    out = join(out_dir, '{}.tblout'.format(prefix))
    cmscan = Dumpling('cmscan', params=Parameters(*cmscan_params))
    cmscan.update(query=query, db=db, out=out, **kwargs)
    logger.info('Running {}'.format(cmscan.command))
    cmscan()
    return cmscan
