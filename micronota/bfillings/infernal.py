# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from burrito.parameters import FlagParameter, ValuedParameter

from .util import _get_parameter
from .model import ModelFetch, ModelPress, ModelScan


class CMScan(ModelScan):
    '''cmscan application controller.

    cmscan is used to search a sequence against a covariance model database.
    This wrapper is tested for INFERNAL 1.1.1 (July 2014)
    '''
    _command = "cmscan"
    _valued_path_options = ModelScan._valued_path_options
    _valued_nonpath_options = [
        # set filters to defaults used for a search space of size Mb
        '--FZ',
        # with --mid, set P-value threshold for HMM stages to  [0.02]
        '--Fmid',
        # set max allowed alnment mx size to Mb [df: autodetermined]
        '--mxsize',
        # set max allowed size of search DP matrices to Mb  [128.]
        '--smxsize',
        # set W (expected max hit len) as * cm->clen (model len)
        '--wcx',
        # configure CMs listed in file in glocal mode, others in local
        '--glist',
    ] + ModelScan._valued_nonpath_options

    _flag_options = [
        # Options controlling acceleration heuristics*:
        # turn all heuristic filters off (slow)
        '--max',
        # skip all HMM filter stages, use only CM (slow)
        '--nohmm',
        # skip first two HMM filter stages (SSV & Vit)
        '--mid',
        # default: run search space size-dependent pipeline  [default]
        '--default',
        # set heuristic filters at Rfam-level (fast)
        '--rfam',
        # use HMM only, don't use a CM at all
        '--hmmonly',
        # do not allow truncated hits at sequence termini
        '--notrunc',
        # allow truncated hits anywhere within sequences
        '--anytrunc',
        # turn off the NULL3 post hoc additional null model
        '--nonull3',
        # use scanning CM CYK algorithm, not Inside in final stage
        '--cyk',
        # align hits with CYK, not optimal accuracy
        '--acyk',
        # only search the top strand
        '--toponly',
        # only search the bottom strand
        '--bottomonly'
    ] + ModelScan._flag_options

    _parameters = {}
    _parameters.update({
        i: _get_parameter(
            ValuedParameter, i, Delimiter=' ', IsPath=True)
        for i in _valued_path_options})
    _parameters.update({
        i: _get_parameter(
            ValuedParameter, i, Delimiter=' ')
        for i in _valued_nonpath_options})
    _parameters.update({
        i: _get_parameter(FlagParameter, i)
        for i in _flag_options})


class CMFetch(ModelFetch):
    '''cmfetch application controller.

    cmfetch is used to get a covariance model by name or accession
    from a CM database.
    This wrapper is tested for INFERNAL 1.1.1 (July 2014).
    '''
    pass


class CMPress(ModelPress):
    '''cmpress application controller.

    cmpress is used to format a CM database into a binary format for cmscan.
    This wrapper is tested for INFERNAL 1.1.1 (July 2014).
    '''
    _command = 'cmpress'


def cmpress_cm(cm, force=False):
    '''Compress the CM database.

    Parameters
    ----------
    cm : str
        The file path to CM database.
    force : boolean
        Whether to overwrite.'''

    app = CMPress(InputHandler='_input_as_path')
    if force is True:
        app.Parameters['-F'].on()
    res = app(cm)
    return res


def cmscan_fasta(cm, in_fp, out_fp, evalue=0.01, cores=0, params=None):
    '''Scan a fasta file against a covariance model database.

    Parameters
    ----------
    cm : str
        The file path to CM database.
    in_fp : str
        Input fasta file. It can contain multiple sequences.
    out_fp : str
        Output file path of target hits table.
    cores : int
        Number of CPU cores. Default to zero, i.e. running in serial-only mode.
    evalue : float
        Default to 0.01. Threshold E-value.
    params : dict
        Other command line parameters for cmscan. key is the option
        (e.g. "-T") and value is the value for the option (e.g. "50").
        If the option is a flag, set the value to None.

    Returns
    -------
    burrito.util.CommandLineAppResult
        It contains opened file handlers of stdout, stderr, and the 3
        output files, which can be accessed in a dict style with the
        keys of "StdOut", "StdErr", "--tblout". The exit status
        of the run can be similarly fetched with the key of "ExitStatus".
    '''
    app = CMScan(InputHandler='_input_as_paths', params=params)
    app.Parameters['--incE'].on(evalue)
    app.Parameters['--cpu'].on(cores)
    app.Parameters['--tblout'].on(out_fp)
    return app([cm, in_fp])
