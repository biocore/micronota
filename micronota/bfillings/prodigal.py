# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, basename, splitext
import logging
import re

from skbio import read
from skbio.metadata import Feature
from skbio.io.format.genbank import _parse_features
from burrito.parameters import FlagParameter, ValuedParameter
from burrito.util import CommandLineApplication, ResultPath

from ._base import IntervalMetadataPred
from ..parsers.embl import _parse_records


class Prodigal(CommandLineApplication):
    '''Prodigal (version 2.6.2) application controller.'''
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
        return exit_status == 0

    def _get_result_paths(self, data):
        result = {}
        for i in self._valued_path_options:
            if i != '-i':
                o = self.Parameters[i]
                if o.isOn():
                    out_fp = self._absolute(o.Value)
                    result[i] = ResultPath(Path=out_fp, IsWritten=True)
        return result


class FeaturePred(IntervalMetadataPred):
    def _identify_features_fp(self, fp, params=None):
        '''Predict genes for the input sequence with Prodigal.'''
        res = self.run(fp, params)
        return self.parse_result(res)

    def run(self, fp, params=None):
        '''Predict genes for the input file.

        Notes
        -----
        It will create 3 output files:
          1. the annotation file in format of GFF3 or
             GenBank feature table with .gbk suffix.
          2. the nucleotide sequences for each predicted gene
             with file suffix of .fna.
          3. the protein sequence translated from each gene
             with file suffix of .faa.

        Parameters
        ----------
        params : dict
            Other command line parameters for Prodigal. key is the option
            (e.g. "-p") and value is the value for the option (e.g. "single").
            If the option is a flag, set the value to None.

        Returns
        -------
        ``burrito.util.CommandLineAppResult``
            It contains opened file handlers of stdout, stderr, and the 3
            output files, which can be accessed in a dict style with the
            keys of "StdOut", "StdErr", "-o", "-d", "-a". The exit status
            of the run can be similarly fetched with the key of "ExitStatus".
        '''
        logger = logging.getLogger(__name__)

        if params is None:
            params = {}

        # default output is genbank
        f_param = params.get('-f', 'gbk')

        out_suffices = {
            '-a': 'faa',
            # output file of nucleotide sequences of genes
            '-d': 'fna',
            '-o': f_param}
        out_prefix = splitext(basename(fp))[0]
        for i in out_suffices:
            out_fp = join(self.out_dir,
                          '.'.join([out_prefix, out_suffices[i]]))
            params[i] = out_fp

        params['-i'] = fp
        app = Prodigal(params=params)
        logger.info('Running: %s' % app.BaseCommand)
        return app()

    def parse_result(self, res, which='-a'):
        '''Parse gene prediction result from ``Prodigal``.

        It is parsed into a dict consumable by
        ``skbio.metadata.IntervalMetadata``.

        Parameters
        ----------
        res : burrito.util.CommandLineAppResult
        which : which output to parse
        Returns
        -------
        ``skbio.metadata.IntervalMetadata``
        '''
        # make sure to move to the beginning of the file.
        if which == '-a':
            return self._parse_faa(res[which])
        elif which == '-o':
            return _parse_records(res[which], self._parse_single_record)

    @staticmethod
    def _parse_single_record(chunks):
        '''Parse single record of Prodigal GenBank output.

        Parameters
        ----------
        chunks : list of str
            a list of lines of the record to parse.

        Yields
        ------
        dict passable to ``skbio.metadata.IntervalMetadata``
        '''
        # get the head line
        head = chunks[0]
        _, description = head.split(None, 1)
        pattern = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
        desc = {}
        for i in pattern.findall(description):
            k, v = i.split('=', 1)
            desc[k] = v
        return _parse_features(chunks[1:], int(desc['seqlen']))

    @staticmethod
    def _parse_faa(faa):
        '''Parse the faa output of Prodigal.

        Notice that Prodigal do not predict genes with introns.

        Yields
        ------
        dict passable to ``skbio.metadata.IntervalMetadata``.
        '''
        pattern = (r'# +([0-9]+)'    # start
                   ' +# +([0-9]+)'   # end
                   ' +# +(-?1)'      # strand
                   ' +# +ID=([0-9]+_[0-9]+);'
                   'partial=([01]{2});'
                   '(.*)')
        im = dict()
        i = 1
        for seq in read(faa, format='fasta'):
            desc = seq.metadata['description']
            matches = re.match(pattern, desc)
            start, end, strand, id, partial, misc = matches.groups()
            # ordinal number of the parent seq
            ordinal = int(id.split('_', 1)[0])
            if ordinal > i:
                yield im
                # reset
                i += 1
                im = dict()
            interval = []
            feature = dict()
            feature['translation'] = str(seq)
            feature['type_'] = 'CDS'

            # don't forget to convert 0-based
            interval = [(int(start)-1, int(end))]
            feature['note'] = '"%s"' % misc
            feature['id'] = id
            if partial[0] == '0':
                feature['left_partial_'] = False
            else:
                feature['left_partial_'] = True
                start = '<%s' % start
            if partial[1] == '0':
                feature['right_partial_'] = False
            else:
                feature['right_partial_'] = True
                end = '>%s' % end
            location = '{s}..{e}'.format(s=start, e=end)
            if strand == '-1':
                feature['rc_'] = True
                location = 'complement(%s)' % location
            elif strand == '1':
                feature['rc_'] = False
            else:
                raise ValueError('Inappropriate value for strand: %s' % strand)
            feature['location'] = location
            im[Feature(**feature)] = interval

        if im:
            # don't forget to return the last one if it is not empty.
            yield im
