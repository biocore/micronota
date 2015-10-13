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


class Prodigal(CommandLineApplication):
    '''Prodigal (version 2.6.2) application controller.
    '''
    _command = 'prodigal'
    _valued_path_options = [
        # Specify FASTA/Genbank input file (default reads from stdin).
        '-i',
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
            IsPath=True)
        for i in _valued_path_options})
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


def predict_genes(in_fp, out_dir, prefix, params=None):
    '''Predict genes for the input file.

    It will create 3 output:
      1. the annotation file in format of GFF3 or
         GenBank feature table.
      2. the nucleotide sequences for each predicted gene.
      3. the protein sequence translated from each gene.

    Parameters
    ----------
    in_fp : str
        input file path
    out_dir : str
        output directory
    prefix : str
        prefix of output file name (without suffix)
    params : dict
        Other command line parameters for Prodigal. key is the option
        (e.g. "-p") and value is the value for the option (e.g. "single").
        If the option is a flag, set the value to None.

    Returns
    -------
    burrito.util.CommandLineAppResult
    '''
    # create dir if not exist
    makedirs(out_dir, exist_ok=True)

    if params is None:
        params = {}

    # default is genbank
    f_param = params.get('-f', 'gbk')

    out_suffices = {
        '-a': 'faa',
        # Write nucleotide sequences of genes to the selected file.
        '-d': 'fna',
        '-o': f_param}
    for i in out_suffices:
        out_fp = join(out_dir, '.'.join([prefix, out_suffices[i]]))
        params[i] = out_fp

    params['-i'] = in_fp

    app = Prodigal(params=params)
    return app()
