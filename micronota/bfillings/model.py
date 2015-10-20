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


class ModelFetch(CommandLineApplication):
    '''Base class for application controller for hmmfetch and cmfetch
    '''
    pass


class ModelPress(CommandLineApplication):
    '''Base class for application controller for hmmpress and cmpress
    '''
    _suppress_stderr = False
    _command = ""
    _parameters = {
        # force overwrite (cmscan)
        '-F': FlagParameter(Prefix='-', Name='F'),
        # force overwrite (hmmer)
        '-f': FlagParameter(Prefix='-', Name='f')}

    def _accept_exit_status(self, exit_status):
        return exit_status == 0


class ModelScan(CommandLineApplication):
    '''Base class for application controller for hmmscan and cmscan
    '''
    _command = ""
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
        # assert query <seqfile> is in format: no autodetection
        '--qformat',
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
    ]

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
