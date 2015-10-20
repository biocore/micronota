# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import remove, close
from tempfile import mkstemp

from skbio import DNA, RNA
from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import CommandLineApplication, ResultPath

from .util import _get_parameter
from .model import ModelFetch, ModelPress, ModelScan


class HMMScan(ModelScan):
    '''hmmscan application controller.

    hmmscan is used to search a sequence against a covariance model database.
    This wrapper is tested for HMMER 3.1b1 (May 2013)
    '''
    _command = "hmmscan"

    _valued_path_options = ModelScan._valued_path_options
    _valued_nonpath_options = [
        # MSV threshold: promote hits w/ P <= F1  [0.02]
        '--F1',
        # Vit threshold: promote hits w/ P <= F2  [1e-3]
        '--F2',
        # Fwd threshold: promote hits w/ P <= F3  [1e-5]
        '--F3',
        # turn off composition bias filter
        '--nobias'
    ] + ModelScan._valued_nonpath_options

    _flag_options = [
        # Options controlling acceleration heuristics*:
        # turn all heuristic filters off (slow)
        '--max',
        # skip all HMM filter stages, use only HMM (slow)
        '--nohmm',
        # skip first two HMM filter stages (SSV & Vit)
        '--mid',
        # default: run search space size-dependent pipeline  [default]
        '--default',
        # set heuristic filters at Rfam-level (fast)
        '--rfam',
        # use HMM only, don't use a HMM at all
        '--hmmonly',
        # do not allow truncated hits at sequence termini
        '--notrunc',
        # allow truncated hits anywhere within sequences
        '--anytrunc',
        # turn off the NULL3 post hoc additional null model
        '--nonull3',
        # use scanning HMM CYK algorithm, not Inside in final stage
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


class HMMFetch(ModelFetch):
    '''hmmfetch application controller.

    hmmfetch is used to get a covariance model by name or accession
    from a HMM database.
    This wrapper is tested for HMMER 3.1b1 (May 2013).
    '''
    pass


class HMMPress(ModelPress):
    '''hmmpress application controller.

    hmmpress is used to format a HMM database into a binary format for hmmscan.
    This wrapper is tested for HMMER 3.1b1 (May 2013).
    '''
    _command = 'hmmpress'


def hmmpress_hmm(hmm, force=False):
    '''Compress the HMM database.

    Parameters
    ----------
    hmm : str
        The file path to HMM database.
    force : boolean
        Whether to overwrite.'''

    app = HMMPress(InputHandler='_input_as_path')
    if force is True:
        app.Parameters['-f'].on()
    res = app(hmm)
    return res


def hmmscan_fasta(hmm, in_fp, out_fp, evalue=0.01, cores=0, params=None):
    '''Scan a fasta file against a covariance model database.

    Parameters
    ----------
    hmm : str
        The file path to HMM database.
    in_fp : str
        Input fasta file. It can contain multiple sequences.
    out_fp : str
        Output file path of target hits table.
    cores : int
        Number of CPU cores. Default to zero, i.e. running in serial-only mode.
    evalue : float
        Default to 0.01. Threshold E-value.
    params : dict
        Other command line parameters for hmmscan. key is the option
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
    app = HMMScan(InputHandler='_input_as_paths', params=params)
    app.Parameters['--incE'].on(evalue)
    app.Parameters['--cpu'].on(cores)
    app.Parameters['--tblout'].on(out_fp)
    return app([hmm, in_fp])
