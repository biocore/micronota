# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import CommandLineApplication, ResultPath
from skbio import read


class Prodigal(CommandLineApplication):
    '''Prodigal (version 2.6.2) application controller.
    '''
    _command = 'prodigal'
    _valued_path_options = [
        # Write protein translations to the selected file.
        '-a',
        # Write nucleotide sequences of genes to the selected file.
        '-d',
        # Write all potential genes (with scores) to the selected file.
        '-s',
        # Write a training file (if none exists);
        # otherwise, read and use the specified training file.
        '-t',
        # Specify output file (default writes to stdout).
        '-o',
    ]
    _valued_nonpath_options = [
        # Select output format (gbk, gff, or sco).  Default is gbk.
        '-f',
        # Specify a translation table to use (default 11).
        '-g',
        # Specify FASTA/Genbank input file (default reads from stdin).
        '-i',
        # Treat runs of N as masked sequence; don't build genes across them.
        '-m',
        # Select procedure (single or meta).  Default is single.
        '-p',
    ]
    _flag_options = [
        # Closed ends.  Do not allow genes to run off edges.
        '-c',
        # Print version number and exit.
        '-v',
        # Run quietly (suppress normal stderr output).
        '-q',
        # Print help menu and exit.
        '-h',
        # Bypass Shine-Dalgarno trainer and force a full motif scan.
        '-n',
    ]

    _parameters = {}
    _parameters.update({
        i: ValuedParameter(
            Prefix=i[0], Name=i[1:], Delimiter=' ',
            IsPath=True if i in _valued_path_options else False)
        for i in _valued_path_options + _valued_nonpath_options})
    _parameters.update({
        i: FlagParameter(
            Prefix=i[0], Name=i[1:])
        for i in _flag_options})
    _suppress_stderr = False

    def _accept_exit_status(self, exit_status):
        if exit_status == 0:
            return True
        else:
            return False

    def _get_result_paths(self, data):
        result = {}
        for i in ['-a', '-d', '-t', '-o', '-s']:
            o = self.Parameters[i]
            if o.isOn():
                out_fp = self._absolute(o.Value)
                result[i] = ResultPath(Path=out_fp, IsWritten=True)
        return result


def predict_genes(params):
    '''Predict genes for the input file.

    Parameters
    ----------
    params : dict
        Command line parameters for Prodigal. key is the option (e.g. "-p")
        and value is the value for the option (e.g. "meta"). If the option
        is a flag, set the value to None.

    Returns
    -------
    generator of `skbio.DNA`
    '''
    # this converts the format specification of prodigal to that of skbio
    formats = {'gbk': 'genbank'}
    # default is genbank
    f_param = params.get('-o', 'gbk')
    f = formats.get(f_param, f_param)
    app = Prodigal(params=params)
    res = app()
    pred = res.get('-o', None)
    if pred is None:
        # the output directed to stdout
        pred = res['StdOut']
    return read(pred, format=f)
