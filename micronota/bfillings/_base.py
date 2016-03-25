# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from abc import ABC, abstractmethod
from tempfile import mkdtemp, NamedTemporaryFile
from inspect import signature

from pandas import DataFrame
from skbio import Sequence


class SubclassImplementError(Exception):
    '''Raised when a subclass do not follow the enforcement.'''
    def __init__(self, cls, message=('This class definition violates '
                                     'the enforced rule of its parent class')):
        super().__init__('%s: %s' % (message, cls))


class IntervalMetadataPred(ABC):
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
    @classmethod
    def __subclasshook__(cls, C):
        '''Enforce the API of functions in child classes.'''
        if cls is IntervalMetadataPred:
            f = C.__dict__['_identify_fp']
            sig = signature(f)
            # enforce it to return dict
            if not issubclass(sig.return_annotation, dict):
                raise SubclassImplementError(C)
        return True

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

    def __call__(self, input, **kwargs) -> dict:
        '''Identify features for the input.

        Parameters
        ----------
        input : ``skbio.Sequence`` or sequence file.

        Yield
        -----
        dict passable to ``skbio.metadata.IntervalMetadata``
        '''
        if isinstance(input, Sequence):
            return self._identify_seq(input, **kwargs)
        elif isinstance(input, str):
            return self._identify_fp(input, **kwargs)

    def _identify_seq(self, seq, **kwargs):
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
    def _identify_fp(self, fp, **kwargs):
        '''Identify features on the sequence in the input file.'''

    def has_cache(self):
        return self.cache is not None


class MetadataPred(ABC):
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
    @classmethod
    def __subclasshook__(cls, C):
        '''Enforce the API of functions in child classes.'''
        if cls is MetadataPred:
            f = C.__dict__['_annotate_fp']
            sig = signature(f)
            # enforce it to return dict
            if not issubclass(sig.return_annotation, DataFrame):
                raise SubclassImplementError(C)
        return True

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

    def __call__(self, input, **kwargs):
        '''
        Parameters
        ----------
        input : list of ``skbio.Sequence`` or sequence file.

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
            self._annotate_fp(f.name, **kwargs)

    def has_cache(self):
        return self.cache is not None

    @abstractmethod
    def _annotate_fp(self, fp, **kwargs):
        '''Add metadata to the sequences in the input file.

        Parameters
        ----------
        fp : input file of sequences
        '''
