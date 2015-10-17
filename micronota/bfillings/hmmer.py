# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import join, isabs

from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import CommandLineApplication, ResultPath
from .util import _get_parameter


class HMMPress(CommandLineApplication):
    '''hmmpress application controller.
    hmmpress is used to format a HMM database into a binary format for hmmscan.
    This wrapper is tested for INFERNAL 1.1.1 (July 2014).
    '''
    _suppress_stderr = False
    _command = 'hmmpress'
    _parameters = {
        # force overwrite
        '-F': FlagParameter(Prefix='-', Name='F')}

    def _accept_exit_status(self, exit_status):
        return exit_status == 0


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
        app.Parameters['-F'].on()
    res = app(hmm)
    return res


class HMMScan(CommandLineApplication):
    '''HMMER (version 3.1) application controller.'''
    _command = 'hmmscan'

    _input_handler = '_input_as_list'

    _valued_path_options = [
        # direct output to file, not stdout
        '-o',
        # save parseable table of per-sequence hits to file <s>
        '--tblout',
        # save parseable table of per-domain hits to file <s>
        '--domtblout'
        ]
    _valued_nonpath_options = [
        # set max width of ASCII text output lines  [120]  (n>=120)
        '--textw',

        ## options controlling reporting thresholds:
        # report models <= this E-value threshold in output  [10.0]  (x>0)
        '-E',
        # report models >= this score threshold in output
        '-T',
        # report domains <= this E-value threshold in output  [10.0]  (x>0)
        '--domE',
        # report domains >= this score cutoff in output
        '--domT',

        ## options controlling inclusion (significance) thresholds:
        # consider models <= this E-value threshold as significant
        '--incE',
        # consider models >= this score threshold as significant
        '--incT',
        # consider domains <= this E-value threshold as significant
        '--incdomE',
        # consider domains >= this score threshold as significant
        '--incdomT',

        ## options controlling acceleration heuristics:
        # MSV threshold: promote hits w/ P <= F1  [0.02]
        '--F1',
        # Vit threshold: promote hits w/ P <= F2  [1e-3]
        '--F2',
        # Fwd threshold: promote hits w/ P <= F3  [1e-5]
        '--F3',

        ## other expert options:
        # set # of comparisons done, for E-value calculation
        '-Z',
        # set # of significant seqs, for domain E-value calculation
        '--domZ',
        # set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
        '--seed',
        # assert input <seqfile> is in format <s>: no autodetection
        '--qformat',
        # number of parallel CPU workers to use for multithreads
        '--cpu'
        ]

    _flag_options = [
        ## options controlling output:
        # prefer accessions over names in output
        '--acc',
        #don't output alignments, so output is smaller
        '--noali',
        # unlimit ASCII text output line width
        '--notextw',

        ## options for model-specific thresholding:
        # use profiles GA gathering cutoffs to set all thresholding
        '--cut_ga',
        # use profiles NC noise cutoffs to set all thresholding
        '--cut_nc',
        # use profiles TC trusted cutoffs to set all thresholding
        '--cut_tc',

        ## options controlling acceleration heuristics:
        # Turn all heuristic filters off (less speed, more power)
        '--max',

        # turn off composition bias filter
        '--nobias',

        ## other expert options:
        # turn off biased composition score corrections
        '--nonull2']

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

    def _input_as_list(self, data):
        '''Takes the positional arguments as input in a list.
        The list input here should be [query_file_path, database_file_path]
        '''
        query, database = data
        if (not isabs(database)) \
                or (not isabs(query)):
            raise ApplicationError("Only absolute paths allowed.\n%s" %
                                   ', '.join(data))

        self._database = FilePath(database)
        self._query = FilePath(query)
        return ''

    def _accept_exit_status(self, exit_status):
        return exit_status == 0


def hmmscan_fasta(hmmdb, in_fp, out_fp, evalue=0.01, cores=0, params=None):
    '''Scan fasta file against a covariance model database.
    Parameters
    ----------
    hmmdb : str
        Path to the hmm database
    in_fp : str
        Input fasta file.
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
    app = HMMScan(InputHandler='_input_as_paths', params=params)
    app.Parameters['--incE'].on(evalue)
    app.Parameters['--cpu'].on(cores)
    app.Parameters['--tblout'].on(out_fp)
    return app([hmmdb, in_fp])
