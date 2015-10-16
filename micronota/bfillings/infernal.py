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


class CMScan(CommandLineApplication):
    '''cmscan application controller.

    cmscan is used to search a sequence against a covariance model database.
    This wrapper is tested for INFERNAL 1.1.1 (July 2014)
    '''
    _command = "cmscan"
    _suppress_stderr = False

    _valued_path_options = [
        # -o <f> Save the output to a file <f>. The default
        # is to write it to standard output.
        '-o',
        # save parseable table of hits to file
        '--tblout'
    ]
    _valued_nonpath_options = [
        # report sequences <= this E-value threshold in output [10.0]  (x>0)
        '-E',
        # report sequences >= this score threshold in output
        '-T',
        # consider sequences <= this E-value threshold as significant  [0.01]
        '--incE',
        # consider sequences >= this score threshold as significant
        '--incT',
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
        # assert query <seqfile> is in format: no autodetection
        '--qformat',
        # configure CMs listed in file in glocal mode, others in local
        '--glist',
        # number of parallel CPU workers to use for multithreads
        '--cpu'
    ]
    _flag_options = [
        # report extra information; mainly useful for debugging
        '--verbose',
        # prefer accessions over names in output
        '--acc',
        # dont output alignments, so output is smaller
        '--noali',
        # unlimit ASCII text output line width
        '--notextw',
        # use CM's GA gathering cutoffs as reporting thresholds
        '--cut_ga',
        # use CM's NC noise cutoffs as reporting thresholds
        '--cut_nc',
        # use CM's TC trusted cutoffs as reporting thresholds
        '--cut_tc',
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
    ]
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

    def _accept_exit_status(self, exit_status):
        return exit_status == 0

    def _get_result_paths(self, data):
        result = {}
        for i in self._valued_path_options:
            o = self.Parameters[i]
            if o.isOn():
                out_fp = self._absolute(o.Value)
                result[i] = ResultPath(Path=out_fp, IsWritten=True)
        return result


class CMFetch(CommandLineApplication):
    '''cmfetch application controller.

    cmfetch is used to get a covariance model by name or accession
    from a CM database.
    This wrapper is tested for INFERNAL 1.1.1 (July 2014).
    '''
    pass


class CMPress(CommandLineApplication):
    '''cmpress application controller.

    cmpress is used to format a CM database into a binary format for cmscan.
    This wrapper is tested for INFERNAL 1.1.1 (July 2014).
    '''
    _suppress_stderr = False
    _command = 'cmpress'
    _parameters = {
        # force overwrite
        '-F': FlagParameter(Prefix='-', Name='F')}

    def _accept_exit_status(self, exit_status):
        return exit_status == 0


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


def cmscan_sequence(cm, seq, evalue=0.01, cores=0, params=None):
    '''scan a `skbio.Sequence` against a covariance model database.

    Parameters
    ----------
    cm : str
        The file path to CM database.
    seq : nucleotide str, skbio.Sequence, or its child classes
        If it is a str, it must be able to read into `skbio.DNA` or
        `skbio.RNA`.
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
    Same data type as seq with added ncRNA annotation as positional_metadata.
    '''
    if isinstance(seq, str):
        try:
            seq = DNA(seq)
        except ValueError:
            seq = RNA(seq)
    app = CMScan(InputHandler='_input_as_paths', params=params)
    app.Parameters['--incE'].on(evalue)
    app.Parameters['--cpu'].on(cores)
    table_fd, table_fp = mkstemp()
    seq_fd, seq_fp = mkstemp()
    with open(seq_fp, 'w') as f:
        seq.write(f, format='fasta')
    app.Parameters['--tblout'].oin(table_fp)
    res = app([cm, seq_fp])
    # todo: read table_fp and merge it into seq object
    res['--tblout']
    # remove temp files
    remove(table_fp)
    remove(seq_fp)
    # don't forget to close the file descriptors
    close(table_fd)
    close(seq_fd)
    return seq
