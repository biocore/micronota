# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import join

from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import CommandLineApplication, ResultPath


class MinCED(CommandLineApplication):
    '''MinCED (version 0.2.0) application controller.'''
    _command = 'minced'
    _valued_nonpath_options = [
        # Length of search window used to discover CRISPRs (range: 6-9).
        # Default: 8
        '-searchWL',
        # Minimum number of repeats a CRISPR must contain. Default: 3
        '-minNR',
        # Minimum length of the CRISPR repeats. Default: 23
        '-minRL',
        # Maximum length of the CRISPR repeats. Default: 47
        '-maxRL',
        # Minimum length of the CRISPR spacers. Default: 26
        '-minSL',
        # Maximum length of the CRISPR spacers. Default: 50
        '-maxSL',
    ]
    _flag_options = [
        # Output version information
        '--version',
        # Output this handy help message
        '-h',
        # Output summary results in gff format containing
        # only the positions of the CRISPR arrays. Default: false
        '-gff',
        # Output detailed results in gff format containing
        # positions of CRISPR arrays and all repeat units. Default: false
        '-gffFull',
        # Output a fasta formatted file containing the spacers. Default: false
        '-spacers',
    ]

    _parameters = {}
    _parameters.update({
        i: ValuedParameter(
            Prefix=i[0], Name=i[1:], Delimiter=' ')
        for i in _valued_nonpath_options})
    _parameters.update({
        i: FlagParameter(
            Prefix=i[0], Name=i[1:])
        for i in _flag_options})
    _suppress_stderr = False

    def _accept_exit_status(self, exit_status):
        return exit_status == 0

    # this needs re-writing to get results from out_dir
    def _get_result_paths(self, data):
        result = {}
        for i in self._flag_options:
            if i != '-h' and i != '--version':
                o = self.Parameters[i]
                if o.isOn():
                    out_fp = self._absolute(o.Value)
                    result[i] = ResultPath(Path=out_fp, IsWritten=True)
        return result


def predict_crispr(in_fp, out_dir, prefix, params=None):
    '''Predict CRISPRs for the input file.

    Notes
    -----
    It will create 1 or 2 output files, depending on the parameters:
      1. file containing CRIPSR information, including locations of CRISPRs
         and their sequence composition OR
      2. GFF file with short information on CRISPR locations
      3. (OPTIONAL; -spacers flag) Fasta file of predicted CRISPR spacers

    Parameters
    ----------
    in_fp : str
        input file path
    out_dir : str
        output directory
    prefix : str
        prefix of output file name (with suffix)
    params : dict
        Other command line parameters for MinCED. key is the option
        (e.g. "-searchWL") and value is the value for the option (e.g. 6).
        If the option is a flag, set the value to None.

    Returns
    -------
    burrito.util.CommandLineAppResult
        It contains opened file handlers of stdout, stderr, and the
        output files, which can be accessed in a dict style with the
        keys of "StdOut", "StdErr". The exit status
        of the run can be similarly fetched with the key of "ExitStatus".
    '''
    # create dir if not exist
    makedirs(out_dir, exist_ok=True)

    if params is None:
        params = {}

    if '-gff' in params:
        out_suffix = 'gff'
    elif '-gffFull' in params:
        out_suffix = 'gffFull'
    else:
        out_suffix = 'crisprs'

    out_fp = join(out_dir, '.'.join([prefix, out_suffix]))

    app = MinCED(InputHandler='_input_as_paths', params=params)
    return app([in_fp, out_fp])
