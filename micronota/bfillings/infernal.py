# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import remove
from tempfile import mkstemp

from burrito.parameters import FlagParameter, ValuedParameter, FilePath
from burrito.util import CommandLineApplication, ResultPath

from .util import _get_parameter


class CMScan(CommandLineApplication):
    '''cmscan application controller.

    INFERNAL 1.1.1 (July 2014)'''
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
        '--cpu',
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

    def _get_result_paths(self,data):
        result = {}
        for i in self._valued_path_options:
            o = self.Parameters[i]
            if o.isOn():
                out_fp = self._absolute(o.Value)
                result[i] = ResultPath(Path=out_fp, IsWritten=True)
        return result


def cmscan_fasta():
    '''Scan fasta file with CM models.'''
    pass
