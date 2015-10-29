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


class _DiamondCommon(CommandLineApplication):
    pass


class DiamondMakeDB(CommandLineApplication):
    _command = 'diamond makedb'
    _parameters = {
        '--in': ValuedParameter(
            Prefix='--', Name='in', Delimiter=' ', IsPath=True),
        '--db': ValuedParameter(
            Prefix='--', Name='db', Delimiter=' ', IsPath=True),
        '--block-size': ValuedParameter(
            Prefix='--', Name='block-size', Delimiter=' '),
        '--threads': ValuedParameter(
            Prefix='--', Name='threads', Delimiter=' ')}
    _suppress_stderr = False

    def _accept_exit_status(self, exit_status):
        return exit_status == 0

    pass


class _DiamondBlastpx(CommandLineApplication):
    _valued_nonpath_options = [
        # -p max Number of CPU threads.
        '--threads',
        # -d Path to DIAMOND database file
        '--db',
        # -q Path to query input file in FASTA or FASTQ format
        # (may be gzip compressed).
        '--query',
        # -a
        '--daa',
        # 11 Gap open penalty.
        '--gapopen',
        # 1 Gap extension penalty.
        '--gapextend',
        # BLOSUM62 Scoring matrix.
        '--matrix',
        # Enable SEG masking of low complexity segments in the query (yes/no).
        # The default is no for blastp and yes for blastx.
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
        '--min-score']
    pass


class DiamondBlastp(CommandLineApplication):
    _command = 'diamond blastp'

    pass
