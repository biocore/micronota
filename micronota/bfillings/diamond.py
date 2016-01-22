# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, basename
from tempfile import mkdtemp

from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import (
    ApplicationError, CommandLineApplication)

from .util import _get_parameter


_OPTIONS_FLAG = {
    i: _get_parameter(FlagParameter, i)
    for i in [
            # only show alignments of forward strand
            '--forwardonly']}

_OPTIONS_VALUE = {
    i: _get_parameter(ValuedParameter, i, Delimiter=' ')
    for i in [
            '--threads',
            # 11 Gap open penalty.
            '--in',
            '--block-size',
            '--gapopen',
            # 1 Gap extension penalty.
            '--gapextend',
            # BLOSUM62 Scoring matrix.
            '--matrix',
            # Enable SEG masking of low complexity segments in the query.
            # (yes/no). The default is no for blastp and yes for blastx.
            '--seg',
            # -k 25 The maximum number of target sequences per query to
            # keep alignments for.
            '--max-target-seqs',
            # Keep alignments within the given percentage range of the top
            # alignment score for a query (overrides –max-target-seqs option).
            '--top',
            # -e 0.001 Maximum expected value to keep an alignment.
            '--evalue',
            # Minimum bit score to keep an alignment. Setting this option
            # will override the --evalue parameter.
            '--min-score',
            # Trigger the sensitive alignment mode with a 16x9 seed
            # shape config.
            '--sensitive',
            # Dynamic programming band for seed extension.
            '--band',
            # Compression for output file (0=none, 1=gzip).]
            '--compress',
            # Format of output file. (tab = BLAST tabular format;
            # sam = SAM format)
            '--outfmt']}

_OPTIONS_PATH = {
    i: _get_parameter(ValuedParameter, i, Delimiter=' ', IsPath=True)
    for i in [
            '--db',
            # -q Path to query input file in FASTA or FASTQ format
            # (may be gzip compressed).
            '--query',
            '--tmpdir',
            # -a
            '--daa',
            # output file
            '--out']}

_PARAMETERS = {}
_PARAMETERS.update(_OPTIONS_FLAG)
_PARAMETERS.update(_OPTIONS_VALUE)
_PARAMETERS.update(_OPTIONS_PATH)
_PARAMETERS_BLAST = {
    i: _PARAMETERS[i]
    for i in
    ['--threads',
     '--db',
     '--query',
     '--tmpdir',
     '--daa',
     '--gapopen',
     '--gapextend',
     '--matrix',
     '--max-target-seqs',
     '--top',
     '--evalue',
     '--min-score',
     '--sensitive',
     '--band',
     '--compress']}


class Diamond(CommandLineApplication):
    '''diamond controller.'''
    _command = 'diamond'
    _suppress_stderr = False

    # This re-implementation is required for subcommands
    def _get_base_command(self):
        if self._subcommand is None:
            raise ApplicationError('_subcommand has not been set.')
        # prevent append multiple subcommand
        if not self._command.endswith(self._subcommand):
            self._command = self._command_delimiter.join(
                [self._command, self._subcommand])
        return super()._get_base_command()

    BaseCommand = property(_get_base_command)

    def _accept_exit_status(self, exit_status):
        return exit_status == 0


class DiamondMakeDB(Diamond):
    '''diamond makedb controller.'''
    _subcommand = 'makedb'
    _parameters = {
        i: _PARAMETERS[i]
        for i in
        ['--in', '--block-size', '--threads', '--db']}


class DiamondBlastp(Diamond):
    '''diamond blastp controller.'''
    _subcommand = 'blastp'
    _parameters = _PARAMETERS_BLAST


class DiamondBlastx(Diamond):
    '''diamond blastp controller.'''
    _subcommand = 'blastx'
    _parameters = _PARAMETERS_BLAST


class DiamondView(Diamond):
    '''diamond view controller.'''
    _subcommand = 'view'
    _parameters = {
        i: _PARAMETERS[i]
        for i in
        ['--daa', '--out', '--outfmt', '--forwardonly']}


def make_db(i, o, params=None):
    '''Format database from a fasta file.

    This is similar to running ``diamond makedb --in db.faa --db db``.

    Parameters
    ----------
    i : str
        Input path for the fasta file.
    o : str
        Output path for the formatted database file.
    params : dict
        Other command line parameters for diamond blastp. key is the option
        (e.g. "-T") and value is the value for the option (e.g. "50").
        If the option is a flag, set the value to None.
    Returns
    -------
    int
        The exit code of the command.
    '''
    app = DiamondMakeDB(InputHandler='_input_as_paths', params=params)
    app.Parameters['--in'].on(i)
    app.Parameters['--db'].on(o)
    res = app()
    res.cleanUp()
    return res['ExitStatus']


def search_protein_homologs(query, db, out_dir, aligner='blastp', outfmt='tab',
                            evalue=0.01, cores=0, params=None):
    '''Search a query seq against the database.

    diamond blastp --db db -q query.faa -o D17.faa.dmd -t tmp_dir -a

    Parameters
    ----------
    query : str
        The file path to the query sequence.
    db : str
        The file path to diamond formatted database.
    cores : int
        Number of CPU cores. Default to 0, i.e. use all available cores.
    evalue : float
        Default to 0.01. Threshold E-value.
    params : dict
        Other command line parameters for diamond blastp. key is the option
        (e.g. "-T") and value is the value for the option (e.g. "50").
        If the option is a flag, set the value to None.
    Returns
    -------
    str
        The file path of the blast result.
    '''
    suffix = basename(query)
    tmpd = mkdtemp(suffix='', prefix='diamond_', dir=out_dir)
    daa_fp = join(out_dir, '%s.daa' % suffix)
    if aligner == 'blastp':
        app = DiamondBlastp
    elif aligner == 'blastx':
        app = DiamondBlastx
    else:
        raise ValueError('unknown aliger')

    blast = app(InputHandler='_input_as_paths', params=params)
    blast.Parameters['--query'].on(query)
    blast.Parameters['--db'].on(db)
    blast.Parameters['--evalue'].on(evalue)
    blast.Parameters['--threads'].on(cores)
    blast.Parameters['--tmpdir'].on(tmpd)
    blast.Parameters['--daa'].on(daa_fp)
    blast_res = blast()
    blast_res.cleanUp()

    out_fp = join(out_dir, '%s.diamond' % suffix)
    view = DiamondView(InputHandler='_input_as_paths', params=params)
    view.Parameters['--daa'].on(daa_fp)
    view.Parameters['--out'].on(out_fp)
    view.Parameters['--outfmt'].on(outfmt)
    view_res = view()
    view_res.cleanUp()
    # print(app.BaseCommand)
    return out_fp
