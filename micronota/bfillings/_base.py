# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from abc import ABCMeta, abstractmethod

from skbio import Sequence


class IntervalMetadataPred(metaclass=ABCMeta):
    @abstractmethod
    def identify_features_on_seq(self, seq, dat, **kwargs):
        '''Identify features on the input sequence.

        Parameters
        ----------
        seq : ``skbio.Sequence`` object

        Returns
        -------
        ``skbio.metadata.IntervalMetadata``
        '''
        if isinstance(seq, Sequence):
            pass

    @abstractmethod
    def identify_features(self, fp, dat, **kwargs):
        '''
        Parameters
        ----------
        fp : sequence file.

        Yield
        -----
        ``skbio.metadata.IntervalMetadata``
            The interval metadata for each sequence.
        '''

    @abstractmethod
    def _parse_result(self):
        '''Parse the result from running an application.'''


class MetadataPred(metaclass=ABCMeta):
    @abstractmethod
    def annotate_features_on_seq(self, seq, dat, **kwargs):
        '''Add metadata to the input seq.

        Assign the function, product, cross-reference, etc. info to the
        input sequence.

        Parameters
        ----------
        seq : ``skbio.Sequence`` object

        Returns
        -------
        ``skbio.metadata``
        '''

    @abstractmethod
    def annotate_features(self, fp, dat, **kwargs):
        '''Add metadata to the input seq.

        Assign the function, product, cross-reference, etc. info to the
        input sequence.

        Parameters
        ----------
        fp : input file of sequences

        Yields
        ------
        ``skbio.metadata``
        '''

    @abstractmethod
    def _parse_result(self):
        '''Parse the result from running an application.'''
