# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from abc import ABCMeta, abstractmethod
from tempfile import mkdtemp

from skbio import Sequence
from skbio.metadata import IntervalMetadata


class IntervalMetadataPred(metaclass=ABCMeta):
    '''
    Attributes
    ----------
    fp : str
        input file path of fasta seq.
    dat : str
        data file (eg database) needed to run the app.
    out_dir : str
        output directory
    tmp_dir : str
        temp directory
    '''
    def __init__(self, dat, out_dir, tmp_dir=None):
        self.dat = dat
        self.out_dir = out_dir
        # create dir if not exist
        makedirs(self.out_dir, exist_ok=True)
        if tmp_dir is None:
            self.tmp_dir = mkdtemp(prefix='tmp', dir=out_dir)
        else:
            self.tmp_dir = tmp_dir
        makedirs(self.tmp_dir, exist_ok=True)

    def identify_features(self, input, **kwargs):
        '''Identify features for the input.

        Parameters
        ----------
        input : ``skbio.Sequence`` or sequence file.

        Yield
        -----
        dict passable to ``skbio.metadata.IntervalMetadata``
        '''
        if isinstance(input, Sequence):
            return self._identify_features_seq(input, **kwargs)
        elif isinstance(input, str):
            return self._identify_features_fp(input, **kwargs)

    def _identify_features_seq(self, seq, **kwargs):
        '''Identify features on the input sequence.

        Parameters
        ----------
        seq : ``skbio.Sequence`` object

        Returns
        -------
        dict passable to ``skbio.metadata.IntervalMetadata``
        '''
        with NamedTemporaryFile('w+', self.tmp_dir) as f:
            seq.write(f)
            self._identify_features_fp(f.name, **kwargs)

    @abstractmethod
    def _identify_features_fp(self, fp, **kwargs):
        '''Identify features on the sequence in the input file.'''


class MetadataPred(metaclass=ABCMeta):
    '''
    Attributes
    ----------
    dat : list of str
        list of data files (eg database) needed to run the app.
    out_dir : str
        output directory
    tmp_dir : str
        temp directory
    '''
    def __init__(self, dat, out_dir, tmp_dir=None):
        self.dat = dat
        self.out_dir = out_dir
        # create dir if not exist
        makedirs(self.out_dir, exist_ok=True)
        if tmp_dir is None:
            self.tmp_dir = mkdtemp(prefix='tmp', dir=out_dir)
        else:
            self.tmp_dir = tmp_dir
        makedirs(self.tmp_dir, exist_ok=True)

    def annotate_features(self, input, **kwargs):
        '''
        Parameters
        ----------
        input : ``skbio.Sequence`` or sequence file.

        Returns
        -------
        pd.DataFrame
            row name should be the query seq id. each column is data
            of e-value, bitscore, etc. For protein sequences, a column
            named 'sseqid' is mandatory to record the seq id of the hit.
        '''
        if isinstance(input, Sequence):
            return self._annotate_seq(input, **kwargs)
        elif isinstance(input, str):
            return self._annotate_fp(input, **kwargs)

    def _annotate_seq(self, seq, **kwargs):
        '''Add metadata to the input seq.

        Assign the function, product, cross-reference, etc. info to the
        input sequence.

        Parameters
        ----------
        seq : ``skbio.Sequence`` object
        '''
        with NamedTemporaryFile('w+', self.tmp_dir) as f:
            seq.write(f)
            self._identify_features_fp(f.name, **kwargs)

    @abstractmethod
    def _annotate_fp(self, fp, **kwargs):
        '''Add metadata to the sequences in the input file.

        Parameters
        ----------
        fp : input file of sequences
        '''
